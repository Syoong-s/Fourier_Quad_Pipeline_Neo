#!/usr/bin/env python3
"""Configure and launch Fourier_Quad using Python 3.8+ and the standard library."""
import argparse
import contextlib
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import uuid

from schema import sections, path_keys, OUTPUT_PATHS, FILE_PATHS

HERE = Path(__file__).resolve().parent
SOURCE_TARGET = '/workspace/src_pipe'
CONFIG_TARGET = '/fq-config/pipeline.ini'
BEGIN = '# >>> FQ MANAGED PATHS'
END = '# <<< FQ MANAGED PATHS'
IMAGE = 'ghcr.io/syoong-s/fourier_quad_pipeline_neo:latest'
ENV_NAME = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*$')
PROTECTED = ('/bin', '/sbin', '/usr', '/lib', '/lib64', '/opt', '/etc', '/proc',
             '/sys', '/dev', '/fq-config', SOURCE_TARGET)


class ConfigError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise ConfigError(message)


def text(value, label):
    require(isinstance(value, str) and '\0' not in value and '\n' not in value
            and '\r' not in value, label + ' must be a single-line string')
    return value


def posint(value, label, minimum=1):
    require(type(value) is int and value >= minimum,
            '%s must be an integer >= %d' % (label, minimum))
    return value


def host_path(value, base, label):
    s = text(value, label)
    require(bool(s), label + ' must not be empty')
    require(not any(c in s for c in ':,'), label + ' must not contain colons or commas; use a path without bind delimiters')
    p = Path(s).expanduser()
    resolved = (p if p.is_absolute() else base / p).resolve()
    require(not any(c in str(resolved) for c in ':,\n\r'), label + ' resolves to a path with unsupported bind delimiters')
    return resolved


def inside(child, parent):
    return child == parent or child.startswith(parent.rstrip('/') + '/')


def container_path(value, label):
    s = text(value, label)
    require(s.startswith('/') and not s.startswith('//'), label + ' must be an absolute container path')
    require(not any(c in s for c in ':,'), label + ' must not contain colons or commas')
    require('..' not in PurePosixPath(s).parts, label + ' must not contain .. path components')
    return str(PurePosixPath(s))


def literal(value):
    # C++ INI accepts quoted strings and quote/backslash escapes, not JSON escapes.
    return '"' + value.replace('\\', '\\\\').replace('"', '\\"') + '"'


def json_no_duplicates(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'Duplicate JSON key: ' + key)
        result[key] = value
    return result


def environment(value, label):
    require(isinstance(value, dict), label + ' must be an object')
    result = {}
    for k, v in value.items():
        require(ENV_NAME.fullmatch(k) is not None, 'Invalid environment variable name: ' + k)
        require(not k.startswith(('SLURM_', 'PMI_', 'PMI2_', 'APPTAINER', 'SINGULARITY', 'FQ_')),
                'Reserved environment variable cannot be overridden: ' + k)
        result[k] = text(v, label + '.' + k)
    return result


def linked_paths(paths):
    """Keep C++ producer/consumer aliases consistent and reuse one mount."""
    p = {k: v for k, v in paths.items() if v is not None}
    producer = 'extcat.output_directory'
    consumer = 'lensing.source_cat'
    require(producer not in p or consumer not in p or p[producer] == p[consumer],
            'lensing.source_cat aliases extcat.output_directory; their values must agree')
    for producer, consumer in (
            ('astrocat.output_directory', 'lensing.astrometry_cat'),
            ('extcat.output_directory', 'lensing.source_cat')):
        if producer not in p and consumer in p:
            p[producer] = p[consumer]
        elif consumer not in p and producer in p:
            p[consumer] = p[producer]
    p.pop('lensing.source_cat', None)  # Write the authoritative path only.
    return p


def load_project(filename, create=False):
    filename = Path(filename).expanduser().resolve()
    require(filename.is_file(), 'JSON not found: %s; run init or copy project.example.json first' % filename)
    raw = filename.read_bytes()
    try:
        obj = json.loads(raw.decode('utf-8'), object_pairs_hook=json_no_duplicates)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ConfigError('Cannot parse JSON: ' + str(exc))
    require(isinstance(obj, dict), 'The JSON root must be an object')
    allowed = {'version', 'variant', 'source', 'image', 'ranks', 'build_jobs',
               'run_directory', 'mounts', 'paths', 'environment', 'hpc'}
    require(not set(obj) - allowed, 'Unknown JSON fields: ' + ', '.join(sorted(set(obj) - allowed)))
    require(type(obj.get('version')) is int and obj['version'] == 1, 'version must be 1')
    variant = obj.get('variant', 'standard')
    require(variant == 'standard', 'Only the standard variant is supported')
    catalogue = sections(variant)
    source = host_path(obj.get('source', '../cpp_Standard'), filename.parent, 'source')
    require(source.is_dir() and (source / 'Makefile').is_file(), 'source must be a C++ directory containing Makefile: ' + str(source))
    require((source / 'src/RuntimeConfig.cpp').is_file(), 'source is missing src/RuntimeConfig.cpp')
    header = source / 'include/RuntimeConfig.hpp'
    require(header.is_file(), 'source is missing include/RuntimeConfig.hpp')
    lens = re.search(r'struct\s+LensingRuntimeConfig\s*\{(.*?)\};', header.read_text(encoding='utf-8'), re.S)
    require(lens is not None, 'Cannot identify LensingRuntimeConfig in source; use a supported source tree')
    require(re.search(r'std::string\s+flat_path\s*;', lens.group(1)) is not None,
            'source is not a supported cpp_Standard tree')
    run_dir = host_path(obj.get('run_directory', './fq-run'), filename.parent, 'run_directory')
    require(not inside(str(run_dir), str(source)) and not inside(str(source), str(run_dir)),
            'run_directory must be separate from source; neither may contain the other')
    image = text(obj.get('image', IMAGE), 'image')
    require(bool(image) and not image.startswith('-') and not any(c.isspace() for c in image), 'Invalid image reference')
    ranks = posint(obj.get('ranks', 4), 'ranks', 2)
    jobs = posint(obj.get('build_jobs', 4), 'build_jobs')
    mounts = obj.get('mounts', [])
    require(isinstance(mounts, list), 'mounts must be an array')
    resolved = []
    targets = set()
    for i, m in enumerate(mounts):
        label = 'mounts[%d]' % i
        require(isinstance(m, dict) and not set(m) - {'host', 'container', 'mode', 'create'}, label + ' has invalid fields')
        src = host_path(m.get('host'), filename.parent, label + '.host')
        dst = container_path(m.get('container'), label + '.container')
        require(dst != '/' and not any(inside(dst, p) or inside(p, dst) for p in PROTECTED),
                'Mount target overlaps a reserved directory: ' + dst)
        require(dst not in targets, 'Duplicate mount target: ' + dst)
        targets.add(dst)
        mode = m.get('mode', 'ro')
        require(mode in ('ro', 'rw'), label + '.mode must be ro or rw')
        mkdir = m.get('create', False)
        require(type(mkdir) is bool and (not mkdir or mode == 'rw'), label + '.create is only supported for rw directories')
        require(not mkdir or not src.exists() or src.is_dir(), label + '.create cannot refer to a regular file')
        if create and mkdir:
            src.mkdir(parents=True, exist_ok=True)
        pending = mkdir and not src.exists()
        require(pending or src.is_file() or src.is_dir(), 'Mount source is missing or is not a regular file/directory: ' + str(src))
        if not pending:
            require(os.access(src, os.R_OK), 'Mount source is not readable: ' + str(src))
            require(mode != 'rw' or os.access(src, os.W_OK), 'rw mount source is not writable: ' + str(src))
        resolved.append({'host': src, 'container': dst, 'mode': mode,
                         'is_file': src.is_file(), 'pending': pending})
    paths = obj.get('paths', {})
    require(isinstance(paths, dict), 'paths must be an object')
    known_paths = path_keys(variant)
    require(not set(paths) - set(known_paths), 'Unknown path keys: ' + ', '.join(sorted(set(paths) - set(known_paths))))
    for key, value in paths.items():
        if value is not None:
            text(value, 'paths.' + key)
    paths = linked_paths(paths)
    for key, value in list(paths.items()):
        if not value:
            continue  # An intentional empty string retains its existing INI meaning.
        if known_paths[key] == 'relative_path' and not value.startswith('/'):
            require('..' not in PurePosixPath(value).parts and '\\' not in value,
                    'Relative paths must not contain .. or backslashes: ' + key)
        else:
            paths[key] = container_path(value, key)
            mount, mapped = resolve_mount(paths[key], resolved)
            require(mount is not None, key + ' is not covered by any mount: ' + paths[key])
            if key in OUTPUT_PATHS:
                require(mount['mode'] == 'rw', key + ' requires a rw mount')
                require(not mount['is_file'], key + ' requires a parent-directory mount so output files can be published or replaced')
            if key not in FILE_PATHS:
                require(not mount['is_file'], key + ' requires a directory, not a regular-file mount')
                require(not mapped.exists() or mapped.is_dir(), key + ' points to a file rather than a directory')
            else:
                require(not mapped.exists() or mapped.is_file(), key + ' requires a file rather than a directory')
    env = {'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'}
    env.update(environment(obj.get('environment', {}), 'environment'))
    hpc = obj.get('hpc', {})
    require(isinstance(hpc, dict) and not set(hpc) - {'sif', 'runtime', 'srun_args', 'environment'}, 'Invalid hpc fields')
    sif = host_path(hpc.get('sif', './cpppipeline.sif'), filename.parent, 'hpc.sif')
    runtime = text(hpc.get('runtime', 'apptainer'), 'hpc.runtime')
    require(bool(runtime) and not runtime.startswith('-'), 'Invalid hpc.runtime')
    srun_args = hpc.get('srun_args', [])
    require(isinstance(srun_args, list), 'hpc.srun_args must be an array of strings')
    # Only placement controls, never a command or an override of task-count / MPI / export.
    for arg in srun_args:
        text(arg, 'hpc.srun_args')
        require(re.fullmatch(r'--(?:cpu-bind|distribution|hint|mem-bind)=[^\s]+', arg) is not None,
                'srun_args accepts --cpu-bind=, --distribution=, --hint= and --mem-bind=; request resources in your own job script')
    henv = dict(env)
    henv.update(environment(hpc.get('environment', {}), 'hpc.environment'))
    return {'filename': filename, 'hash': hashlib.sha256(raw).hexdigest(), 'raw': obj,
            'variant': variant, 'source': source, 'run_dir': run_dir, 'image': image,
            'ranks': ranks, 'jobs': jobs, 'mounts': resolved, 'paths': paths,
            'env': env, 'hpc_env': henv, 'sif': sif, 'runtime': runtime, 'srun_args': srun_args,
            'catalogue': catalogue}


def resolve_mount(path, mounts):
    for mount in sorted(mounts, key=lambda x: len(x['container']), reverse=True):
        if path == mount['container']:
            return mount, mount['host']
        if not mount['is_file'] and inside(path, mount['container']):
            return mount, mount['host'] / path[len(mount['container']):].lstrip('/')
    return None, None


def path_block(p):
    lines = [BEGIN, '# Paths are managed by JSON paths. Commented entries retain compiled defaults.']
    for section, keys in p['catalogue'].items():
        lines += ['', '[' + section + ']']
        for key, (kind, desc) in keys.items():
            full = section + '.' + key
            if kind not in ('path', 'relative_path'):
                continue
            lines.append('# ' + desc)
            lines.append(key + ' = ' + literal(p['paths'][full]) if full in p['paths'] else '# ' + key + ' =')
    return '\n'.join(lines + [END])


def science_block(p):
    lines = ['# Scientific settings: uncomment and fill only the values to override.',
             '# Commented entries retain config/*.hpp defaults. An active empty value does NOT inherit.',
             '# No phases are selected here. Compiled defaults may enable the full pipeline.']
    for section, keys in p['catalogue'].items():
        nonpaths = [(k, kind, desc) for k, (kind, desc) in keys.items()
                    if kind not in ('path', 'relative_path')]
        if not nonpaths:
            continue
        lines += ['', '[' + section + ']']
        for key, kind, desc in nonpaths:
            hints = {'bool': 'true / false', 'datasets': 'e.g. gband:c4d_',
                     'columns': 'e.g. 1,2,3', 'tokens': 'comma-separated',
                     'positive': 'positive integer', 'nonnegative': 'non-negative integer', 'float': 'positive number'}
            hint = hints.get(kind, kind if '|' in kind else '')
            lines += ['# ' + desc + (' [' + hint + ']' if hint else ''), '# ' + key + ' =']
    return '\n'.join(lines) + '\n'


def render_ini(p, previous=None):
    block = path_block(p)
    if previous is None:
        return '# Fourier_Quad configuration: compiled defaults < this INI.\n\n' + block + '\n\n' + science_block(p)
    require(previous.count(BEGIN) == 1 and previous.count(END) == 1,
            'The existing INI has missing or duplicate managed-path markers; configure stopped without replacing it')
    before, rest = previous.split(BEGIN, 1)
    _, after = rest.split(END, 1)
    return before + block + after


def shell_array(name, values):
    return name + '=(' + ' '.join(shlex.quote(str(v)) for v in values) + ')\n'


def render_env(p, generation):
    values = {
        'FQ_ENV_VERSION': '1', 'FQ_GENERATION': generation,
        'FQ_JSON': str(p['filename']), 'FQ_JSON_SHA256': p['hash'],
        'FQ_SOURCE_HOST': str(p['source']), 'FQ_CONFIG_HOST': str(p['run_dir'] / 'pipeline.ini'),
        'FQ_SIF': str(p['sif']), 'FQ_RUNTIME': p['runtime'], 'FQ_BUILD_JOBS': str(p['jobs']),
    }
    result = '# Generated by fqctl.py configure. Edit JSON and regenerate. Source only trusted files.\n'
    result += ''.join(k + '=' + shlex.quote(v) + '\n' for k, v in values.items())
    for key, field in [('FQ_BIND_SOURCES', 'host'), ('FQ_BIND_TARGETS', 'container'), ('FQ_BIND_MODES', 'mode')]:
        result += shell_array(key, [m[field] for m in p['mounts']])
    result += shell_array('FQ_SRUN_ARGS', p['srun_args'])
    result += shell_array('FQ_CONTAINER_ENV', [k + '=' + v for k, v in sorted(p['hpc_env'].items())])
    return result


def atomic_write(path, content):
    fd, tmp = tempfile.mkstemp(prefix='.' + path.name, dir=str(path.parent))
    try:
        with os.fdopen(fd, 'w', encoding='utf-8', newline='\n') as f:
            f.write(content)
            f.flush()
            os.fsync(f.fileno())
        os.chmod(tmp, 0o600)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


@contextlib.contextmanager
def configure_lock(directory):
    import fcntl
    with (directory / '.configure.lock').open('a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        yield


def configure(p):
    p['run_dir'].mkdir(parents=True, exist_ok=True)
    ini = p['run_dir'] / 'pipeline.ini'
    env = p['run_dir'] / 'cpppipeline.env'
    with configure_lock(p['run_dir']):
        previous = ini.read_text(encoding='utf-8') if ini.exists() else None
        generated = render_ini(p, previous)
        generation = uuid.uuid4().hex
        generated = re.sub(r'^# FQ_GENERATION=.*\n', '', generated, flags=re.M)
        generated = '# FQ_GENERATION=' + generation + '\n' + generated
        # Validate before publishing. A user typo is not silently discarded.
        validate_ini(p, generated)
        if previous is not None:
            backup = ini.with_name(ini.name + '.bak.' + str(time.time_ns()))
            atomic_write(backup, previous)
            print('Previous INI saved to:', backup)
        atomic_write(ini, generated)
        atomic_write(env, render_env(p, generation))
    print('INI:', ini)
    print('HPC env:', env)
    print('Non-path defaults are commented placeholders. Review phases, datasets and science settings, then run check.')


def strip_comment(line):
    quote = None
    escaped = False
    for i, c in enumerate(line):
        if escaped:
            escaped = False
            continue
        if c == '\\' and quote:
            escaped = True
            continue
        if c in "\"'":
            if quote is None:
                quote = c
            elif quote == c:
                quote = None
        elif c in '#;' and quote is None:
            return line[:i]
    return line


def decode(value):
    v = strip_comment(value).strip()
    if not v or v[0] not in "\"'":
        return v
    require(len(v) >= 2 and v[-1] == v[0], 'Unmatched string quotes')
    out = ''
    i = 1
    while i < len(v) - 1:
        if v[i] == '\\':
            i += 1
            require(i < len(v) - 1, 'Incomplete escape at end of string')
            if v[i] not in ('\\', v[0]):
                out += '\\'
        out += v[i]
        i += 1
    return out


def split_values(value):
    value = strip_comment(value).strip()
    if not value:
        return []
    parts = []
    start = 0
    quote = None
    escaped = False
    for i, c in enumerate(value + ','):
        if escaped:
            escaped = False
            continue
        if c == '\\' and quote:
            escaped = True
            continue
        if c in "\"'":
            if quote is None:
                quote = c
            elif quote == c:
                quote = None
        elif c == ',' and quote is None:
            v = decode(value[start:i])
            require(bool(v), 'List entries must not be empty')
            parts.append(v)
            start = i + 1
    require(quote is None, 'Unmatched list quotes')
    return parts


def validate_ini(p, content):
    active = {}
    section = None
    for lineno, line in enumerate(content.splitlines(), 1):
        clean = strip_comment(line).strip()
        if not clean:
            continue
        try:
            if clean.startswith('['):
                require(clean.endswith(']'), 'Section header is missing ]')
                section = clean[1:-1].strip().lower()
                require(section in p['catalogue'], 'Unknown section: ' + section)
                continue
            require(section is not None and '=' in clean, 'Expected [section] and key=value')
            key, raw = clean.split('=', 1)
            key = key.strip().lower()
            require(key in p['catalogue'][section], 'Unknown setting: ' + section + '.' + key)
            full = section + '.' + key
            require(full not in active, 'Duplicate assignment: ' + full)
            kind = p['catalogue'][section][key][0]
            value = decode(raw)
            if kind == 'bool':
                require(value.lower() in ('true', 'false', '1', '0', 'yes', 'no', 'on', 'off'), 'Expected a boolean')
            elif kind in ('int', 'positive', 'nonnegative'):
                require(re.fullmatch(r'[+-]?[0-9]+', value) is not None, 'Expected an integer')
                n = int(value)
                require(kind != 'positive' or n > 0, 'Expected a positive integer')
                require(kind != 'nonnegative' or n >= 0, 'Expected a non-negative integer')
                limit = 2**64-1 if section == 'extcat' else (2**63-1 if key == 'process_stage' else 2**31-1)
                require(-2**31 <= n <= limit, 'Integer is out of range')
            elif kind == 'float':
                require(re.fullmatch(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?', value) is not None,
                        'Expected a decimal floating-point number')
                try:
                    n = float(value)
                except ValueError:
                    raise ConfigError('Expected a positive number')
                require(math.isfinite(n) and n > 0, 'Expected a finite positive number')
            elif '|' in kind:
                require(value.lower() in kind.split('|'), 'Allowed values: ' + kind)
            elif kind in ('tokens', 'datasets', 'columns'):
                vs = split_values(raw)
                if kind == 'datasets':
                    targets = set()
                    for item in vs:
                        require(item.count(':') == 1, 'datasets entries must use TARGET:PREFIX')
                        target, prefix = item.split(':')
                        require(target not in ('', '.', '..') and '/' not in target and '\\' not in target
                                and bool(prefix) and target not in targets, 'Invalid or duplicate dataset target/prefix')
                        targets.add(target)
                if kind == 'columns':
                    require(bool(vs) and all(re.fullmatch(r'\+?[0-9]+', v) and int(v) > 0 for v in vs),
                            'input_columns requires a non-empty list of positive integers')
            if full in path_keys(p['variant']):
                require(full in p['paths'] and value == p['paths'][full],
                        'Edit this path in JSON paths, then run configure: ' + full)
            active[full] = value
        except ConfigError as exc:
            raise ConfigError('INI line %d: %s' % (lineno, exc))
    for full, value in p['paths'].items():
        require(active.get(full) == value, 'INI path differs from JSON: ' + full + '; run configure again')
    true = lambda k: active.get(k, '').lower() in ('true', '1', 'yes', 'on')
    for phase, key in [('astrocat', 'astrocat.output_directory'), ('extcat', 'extcat.output_directory')]:
        if true('process.run_process_' + phase) and active.get(key):
            mount, _ = resolve_mount(active[key], p['mounts'])
            require(mount and mount['mode'] == 'rw', key + ' requires a rw mount when generating catalog tiles')
    if active.get('lensing.process_stage'):
        n = int(active['lensing.process_stage'])
        require(n % 23 != 0 or n % 19 == 0, 'process_stage with factor 23 must also have factor 19')
    for key, values in {'ccd_split':(1,2),'astrometry_trivial':(0,1),'astrometry_cat_type':(1,2),
                        'include_flat':(0,1),'include_mask':(0,1,2,3),'ext_cat':(0,1),'ext_psf':(0,1),
                        'deblending':(0,1),'psf_type':(1,2),'psf_ms':(0,1)}.items():
        full = 'lensing.' + key
        require(full not in active or int(active[full]) in values, full + ' has an invalid value')
    phase_keys = ['process.run_process_' + name for name in ('astrocat', 'extcat', 'init', 'main', 'rearr', 'fd')]
    require(not all(k in active and not true(k) for k in phase_keys), 'At least one processing phase must be enabled')
    if 'init.datasets' in active and not active['init.datasets']:
        require(not any(true('process.run_process_' + n) for n in ('init', 'main', 'rearr', 'fd')),
                'datasets must not be empty when a dataset phase is enabled')
    if true('process.run_process_fd'):
        require(active.get('lensing.ext_cat') != '0', 'process_fd requires ext_cat = 1')
    return active


def check_generated(p):
    ini, env = p['run_dir'] / 'pipeline.ini', p['run_dir'] / 'cpppipeline.env'
    require(ini.is_file() and env.is_file(), 'Generated files are missing; run configure first')
    it = ini.read_text(encoding='utf-8')
    et = env.read_text(encoding='utf-8')
    require(('FQ_JSON_SHA256=' + p['hash'] + '\n') in et, 'JSON has changed; run configure to update paths and env')
    match = re.search(r'^# FQ_GENERATION=([a-f0-9]{32})$', it, re.M)
    require(match is not None and ('FQ_GENERATION=' + match.group(1) + '\n') in et,
            'INI and env have different generation identifiers; run configure again')
    validate_ini(p, it)
    for m in p['mounts']:
        require(not m['pending'], 'Output directory not created: %s; run configure' % m['host'])
    print('Configuration check passed: JSON, explicit INI settings and declared mounts are consistent.')
    print('Commented settings retain compiled defaults. Undeclared paths, FITS contents and cluster compatibility are not checked.')


def docker_command(p, action):
    command = ['docker', 'run', '--rm', '--init', '--platform', 'linux/amd64',
               '--user', '%d:%d' % (os.getuid(), os.getgid()),
               '--workdir', SOURCE_TARGET, '--env', 'HOME=/tmp']
    for group in sorted(set(os.getgroups()) - {os.getgid()}):
        command += ['--group-add', str(group)]
    env = dict(p['env'])
    if os.getuid() == 0:
        env.update(OMPI_ALLOW_RUN_AS_ROOT='1', OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
    for key, val in sorted(env.items()):
        command += ['--env', key + '=' + val]
    mounts = [{'host':p['source'], 'container':SOURCE_TARGET, 'mode':'rw' if action == 'build' else 'ro'},
              {'host':p['run_dir'] / 'pipeline.ini', 'container':CONFIG_TARGET, 'mode':'ro'}] + p['mounts']
    for m in mounts:
        spec = 'type=bind,src=%s,dst=%s' % (m['host'], m['container'])
        if m['mode'] == 'ro':
            spec += ',readonly'
        command += ['--mount', spec]
    command.append(p['image'])
    if action == 'run':
        command += ['mpirun', '-np', str(p['ranks']), SOURCE_TARGET + '/Fourier_Quad_Pipe', '--config', CONFIG_TARGET]
    else:
        command += ['bash', '-c', 'set -euo pipefail; make clean; exec make -j"$1" "STACK_PREFIX=${CPP_STACK_PREFIX}"',
                    'fq-build', str(p['jobs'])]
    return command


def execute(command):
    require(shutil.which(command[0]) is not None, 'Command not found: ' + command[0])
    result = subprocess.run(command, check=False)
    return result.returncode if result.returncode >= 0 else 128 - result.returncode


def init_project(path):
    path = Path(path).expanduser().resolve()
    require(not path.exists(), 'Refusing to overwrite existing JSON: ' + str(path))
    obj = json.loads((HERE / 'project.example.json').read_text(encoding='utf-8'))
    obj['source'] = str((HERE.parent / 'cpp_Standard').resolve())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
        f.write('\n')
    print('Created', path)
    print('Fill in host directories, then run configure. Scientific settings are not selected automatically.')


def main(argv=None):
    parser = argparse.ArgumentParser(description='Configure and launch Fourier_Quad in Docker or Apptainer')
    sub = parser.add_subparsers(dest='command', required=True)
    for command, help_text in [('init','Create a JSON example'),('configure','Generate INI and HPC env from JSON'),
                               ('check','Check configuration without starting containers'),('pull','Pull the Docker toolchain image'),
                               ('build','Rebuild the C++ program in Docker'),('run','Run the MPI program with the generated INI')]:
        item = sub.add_parser(command, help=help_text)
        item.add_argument('--json', default=str(HERE / 'project.json'), help='JSON file; relative host paths are resolved against its directory')
    args = parser.parse_args(argv)
    try:
        if args.command == 'init':
            init_project(args.json)
            return 0
        p = load_project(args.json, create=args.command == 'configure')
        if args.command == 'configure':
            configure(p)
        elif args.command == 'check':
            check_generated(p)
        elif args.command == 'pull':
            return execute(['docker', 'pull', '--platform', 'linux/amd64', p['image']])
        else:
            check_generated(p)
            if args.command == 'run':
                require(os.access(p['source'] / 'Fourier_Quad_Pipe', os.X_OK), 'Executable missing; run build first')
            return execute(docker_command(p, args.command))
        return 0
    except (ConfigError, OSError, ValueError) as exc:
        print('Error: ' + str(exc), file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        return 130


if __name__ == '__main__':
    sys.exit(main())
