"""Offline regression tests; no Docker, SIF, Slurm or scientific data required."""
import importlib.util
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import fqctl as fq


class LauncherTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='fq tests ')
        self.root = Path(self.tmp.name)
        self.source = self.root / 'source code'
        (self.source / 'src').mkdir(parents=True)
        (self.source / 'include').mkdir()
        (self.source / 'Makefile').write_text('all:\n\ttrue\n')
        (self.source / 'src/RuntimeConfig.cpp').write_text('// fixture\n')
        self.binary = self.source / 'Fourier_Quad_Pipe'
        self.binary.write_text('#!/bin/sh\nexit 0\n')
        self.binary.chmod(0o755)
        (self.root / 'inputs').mkdir()
        (self.root / 'inputs/list.txt').write_text('fixture\n')
        self.cfg = self.root / 'project.json'
        self.obj = {
            'version':1, 'variant':'standard', 'source':'source code',
            'run_directory':'generated config', 'ranks':4, 'build_jobs':2,
            'mounts':[
                {'host':'inputs','container':'/data/input','mode':'ro'},
                {'host':'outputs','container':'/data/output','mode':'rw','create':True}],
            'paths':{'init.science_root':'/data/input', 'init.output_root':'/data/output'},
            'hpc':{'sif':'image.sif','runtime':'apptainer','srun_args':[]}}
        self.save()

    def tearDown(self):
        self.tmp.cleanup()

    def save(self):
        (self.source / 'include/RuntimeConfig.hpp').write_text(
            'struct LensingRuntimeConfig { std::string flat_path; };\n')
        self.cfg.write_text(json.dumps(self.obj), encoding='utf-8')

    def load(self, create=False):
        self.save()
        return fq.load_project(self.cfg, create=create)

    def configured(self):
        p = self.load(True)
        fq.configure(p)
        return p

    def test_relative_paths_are_json_relative(self):
        p = self.load()
        self.assertEqual(p['source'], self.source)
        self.assertEqual(p['run_dir'], self.root / 'generated config')

    def test_only_configure_creates_explicit_output(self):
        self.load()
        self.assertFalse((self.root / 'outputs').exists())
        self.load(True)
        self.assertTrue((self.root / 'outputs').is_dir())

    def test_missing_input_fails_without_creating_it(self):
        self.obj['mounts'][0]['host'] = 'missing'
        with self.assertRaises(fq.ConfigError): self.load(True)
        self.assertFalse((self.root / 'missing').exists())

    def test_all_nonpath_values_are_commented_blank(self):
        p = self.load()
        ini = fq.render_ini(p)
        active = fq.validate_ini(p, ini)
        self.assertEqual(set(active), set(p['paths']))
        for section, keys in p['catalogue'].items():
            for key, (kind, _) in keys.items():
                if kind not in ('path', 'relative_path'):
                    self.assertIn('# ' + key + ' =\n', ini)
        self.assertNotIn('run_process_init = true', ini)

    def test_reconfigure_preserves_science(self):
        p = self.configured()
        path = p['run_dir'] / 'pipeline.ini'
        original = path.read_text().replace('# run_process_init =\n', 'run_process_init = false\n# my note\n')
        path.write_text(original)
        self.obj['paths']['init.science_root'] = '/data/input/nested'
        p = self.load(True)
        fq.configure(p)
        now = path.read_text()
        self.assertIn('run_process_init = false\n# my note\n', now)
        self.assertIn('science_root = "/data/input/nested"', now)
        self.assertEqual(next(p['run_dir'].glob('pipeline.ini.bak.*')).read_text(), original)
        fq.check_generated(p)

    def test_reconfigure_without_markers_refuses(self):
        p = self.configured()
        f = p['run_dir'] / 'pipeline.ini'
        f.write_text('[init]\n# handwritten\n')
        with self.assertRaises(fq.ConfigError): fq.configure(p)
        self.assertEqual(f.read_text(), '[init]\n# handwritten\n')

    def test_reconfigure_removes_deleted_json_path(self):
        p = self.configured()
        self.obj['paths'].pop('init.science_root')
        p = self.load(True)
        fq.configure(p)
        self.assertIn('# science_root =\n', (p['run_dir'] / 'pipeline.ini').read_text())

    def test_catalog_links_use_single_mount(self):
        self.obj['paths'] = {'astrocat.output_directory':'/data/output/gaia',
                             'extcat.output_directory':'/data/output/source'}
        p = self.load()
        self.assertEqual(len(p['mounts']), 2)
        self.assertEqual(p['paths']['lensing.astrometry_cat'], '/data/output/gaia')
        self.assertNotIn('lensing.source_cat', p['paths'])

    def test_standard_alias_conflict_rejected(self):
        self.obj['paths'].update({'extcat.output_directory':'/data/output/a', 'lensing.source_cat':'/data/output/b'})
        with self.assertRaises(fq.ConfigError): self.load()

    def test_removed_variant_is_rejected(self):
        self.obj['variant'] = 'lite'
        with self.assertRaises(fq.ConfigError):
            self.load()

    def test_null_path_inherits_and_empty_path_clears(self):
        self.obj['paths'].update({'process.fd_expo_list':None, 'process.fd_output_base_directory':''})
        p = self.load()
        ini = fq.render_ini(p)
        self.assertIn('# fd_expo_list =\n', ini)
        self.assertIn('fd_output_base_directory = ""', ini)

    def test_single_file_mount(self):
        self.obj['mounts'].append({'host':'inputs/list.txt','container':'/lists/expo.list','mode':'ro'})
        self.obj['paths']['process.expo_list'] = '/lists/expo.list'
        p = self.load()
        self.assertTrue(p['mounts'][-1]['is_file'])

    def test_output_list_filename_is_file_not_directory(self):
        self.obj['paths']['process.rearranged_expo_list_filename'] = '/data/output/rearr.list'
        self.load()
        self.obj['paths']['process.rearranged_expo_list_filename'] = '/data/input/list.txt'
        with self.assertRaises(fq.ConfigError): self.load()

    def test_file_mount_cannot_contain_subpaths(self):
        self.obj['mounts'].append({'host':'inputs/list.txt','container':'/lists/expo.list','mode':'ro'})
        self.obj['paths']['process.expo_list'] = '/lists/expo.list/child'
        with self.assertRaises(fq.ConfigError): self.load()

    def test_nested_mount_longest_match(self):
        self.obj['mounts'].append({'host':'inputs','container':'/data/output/read-only','mode':'ro'})
        self.obj['paths']['init.output_root'] = '/data/output/read-only'
        with self.assertRaises(fq.ConfigError): self.load()

    def test_multiple_paths_under_one_parent_mount(self):
        self.obj['paths'].update({'lensing.flat_path':'/data/input/flat', 'lensing.psf_path':'/data/input/psf',
                                 'init.dq_root':'/data/input/dq'})
        self.assertEqual(len(self.load()['mounts']), 2)

    def test_all_path_interfaces_supported(self):
        self.obj['paths'] = {k: '/data/output/' + k.replace('.', '/')
                             for k in fq.path_keys('standard')}
        self.obj['paths']['lensing.source_cat'] = self.obj['paths']['extcat.output_directory']
        p = self.load()
        active = fq.validate_ini(p, fq.render_ini(p))
        self.assertEqual(set(active), set(p['paths']))

    def test_duplicate_and_protected_targets_fail(self):
        for target in ('/','/opt','/opt/x','/workspace','/fq-config/x','/data/input'):
            with self.subTest(target=target):
                self.obj['mounts'].append({'host':'inputs','container':target})
                with self.assertRaises(fq.ConfigError): self.load()
                self.obj['mounts'].pop()

    def test_unmapped_path_fails(self):
        self.obj['paths']['init.dq_root'] = '/other/dq'
        with self.assertRaises(fq.ConfigError): self.load()

    def test_unsafe_delimiters_fail(self):
        for path in ('/a:b', '/a,b', '/a\nb'):
            with self.subTest(path=path):
                self.obj['mounts'][0]['host'] = path
                with self.assertRaises(fq.ConfigError): self.load()

    def test_output_cannot_be_readonly(self):
        self.obj['paths']['init.output_root'] = '/data/input'
        with self.assertRaises(fq.ConfigError): self.load()

    def test_generation_outside_source(self):
        self.obj['run_directory'] = 'source code/generated'
        with self.assertRaises(fq.ConfigError): self.load()

    def test_bool_not_accepted_as_rank(self):
        self.obj['ranks'] = True
        with self.assertRaises(fq.ConfigError): self.load()

    def test_stale_json_fails(self):
        p = self.configured()
        self.obj['ranks'] = 8
        p = self.load()
        with self.assertRaises(fq.ConfigError): fq.check_generated(p)

    def test_mismatched_generation_fails(self):
        p = self.configured()
        f = p['run_dir'] / 'pipeline.ini'
        f.write_text(f.read_text().replace('# FQ_GENERATION=', '# removed='))
        with self.assertRaises(fq.ConfigError): fq.check_generated(p)

    def test_check_never_starts_subprocess(self):
        p = self.configured()
        with mock.patch.object(subprocess, 'run', side_effect=AssertionError('must not launch')):
            fq.check_generated(p)

    def test_run_passes_only_config_to_program(self):
        p = self.configured()
        cmd = fq.docker_command(p, 'run')
        i = cmd.index('/workspace/src_pipe/Fourier_Quad_Pipe')
        self.assertEqual(cmd[i+1:], ['--config', '/fq-config/pipeline.ini'])
        self.assertIn('mpirun', cmd)
        self.assertIn(str(p['source']), '\n'.join(cmd))
        self.assertIn('type=bind,src=%s,dst=/workspace/src_pipe,readonly' % p['source'], cmd)
        self.assertNotIn('FQ_PIPELINE_CONFIG', '\n'.join(cmd))

    def test_build_mount_is_writable(self):
        p = self.configured()
        cmd = fq.docker_command(p, 'build')
        self.assertIn('type=bind,src=%s,dst=/workspace/src_pipe' % p['source'], cmd)
        self.assertIn('make clean', cmd[-3])

    def test_env_is_shell_safe(self):
        self.obj['environment'] = {'EXAMPLE':'$(touch UNWANTED),quoted="x"=value'}
        p = self.configured()
        env = p['run_dir'] / 'cpppipeline.env'
        result = subprocess.run(['bash','-c','source "$1"; printf "%s" "${FQ_CONTAINER_ENV[@]}"','check',str(env)],
                                capture_output=True,text=True,check=True,cwd=self.root)
        self.assertIn('$(touch UNWANTED)', result.stdout)
        self.assertFalse((self.root / 'UNWANTED').exists())

    def test_unknown_json_or_ini_fields_fail(self):
        self.obj['rnak'] = 4
        with self.assertRaises(fq.ConfigError): self.load()
        del self.obj['rnak']
        p = self.load()
        with self.assertRaises(fq.ConfigError): fq.validate_ini(p, fq.render_ini(p)+'\n[process]\nrnu_process_init = true\n')

    def test_empty_scalar_fails_but_comment_passes(self):
        p = self.load()
        ini = fq.render_ini(p)
        fq.validate_ini(p, ini)
        for key in ('run_process_init', 'chipnx', 'pixel_size'):
            with self.subTest(key=key):
                with self.assertRaises(fq.ConfigError): fq.validate_ini(p, ini.replace('# '+key+' =\n', key+' =\n'))

    def test_catalog_output_permissions_when_enabled(self):
        self.obj['paths']['astrocat.output_directory'] = '/data/input/tiles'
        p = self.load()
        ini = fq.render_ini(p).replace('# run_process_astrocat =\n','run_process_astrocat = true\n')
        with self.assertRaises(fq.ConfigError): fq.validate_ini(p, ini)

    def test_srun_command_injection_rejected(self):
        self.obj['hpc']['srun_args'] = ['--ntasks=99']
        with self.assertRaises(fq.ConfigError): self.load()
        self.obj['hpc']['srun_args'] = ['--cpu-bind=cores']
        self.load()

    def test_hpc_check_requires_no_runtime_or_sif(self):
        p = self.configured()
        r = subprocess.run(['bash',str(HERE/'runner/run.sh'),'--check',str(p['run_dir']/'cpppipeline.env')],
                           capture_output=True,text=True,cwd='/tmp')
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertFalse(p['sif'].exists())

    def install_stubs(self):
        bindir = self.root / 'bin'
        bindir.mkdir()
        record = self.root / 'record.json'
        for name in ('apptainer', 'srun'):
            f = bindir / name
            f.write_text('#!'+sys.executable+'\nimport json,os,sys\n'
                         'json.dump({"argv":sys.argv,"env":dict(os.environ)},open(os.environ["FQ_TEST_RECORD"],"w"))\n'
                         'sys.exit(int(os.environ.get("FQ_TEST_EXIT","0")))\n')
            f.chmod(0o755)
        env = dict(os.environ, PATH=str(bindir)+os.pathsep+os.environ['PATH'],
                   FQ_TEST_RECORD=str(record), SLURM_JOB_ID='1234', SLURM_NTASKS='8')
        env.pop('SLURM_PROCID',None)
        (self.root/'image.sif').write_text('stub')
        return env, record

    def test_hpc_run_srun_single_launch_and_config(self):
        p = self.configured()
        env, record = self.install_stubs()
        r = subprocess.run(['bash',str(HERE/'runner/run.sh'),str(p['run_dir']/'cpppipeline.env')],
                           capture_output=True,text=True,env=env)
        self.assertEqual(r.returncode,0,r.stderr)
        argv = json.loads(record.read_text())['argv']
        self.assertFalse(any(a.startswith('--ntasks') for a in argv))
        self.assertIn('--mpi=pmi2',argv)
        self.assertEqual(argv[-1],'run')
        self.assertIn(str(HERE/'runner/runtime.sh'),argv)

    def test_per_rank_runtime_uses_rank_env_and_no_eval(self):
        self.obj['environment'] = {'EXAMPLE':'a,b=$(false)'}
        p = self.configured()
        env, record = self.install_stubs()
        env.update(PMI_RANK='5', SLURM_PROCID='5', APPTAINERENV_BAD='$(false)',
                   SINGULARITY_BINDPATH='/bad', APPTAINER_BIND='/bad')
        r = subprocess.run(['bash',str(HERE/'runner/runtime.sh'),str(p['run_dir']/'cpppipeline.env'),'run'],
                           capture_output=True,text=True,env=env)
        self.assertEqual(r.returncode,0,r.stderr)
        result = json.loads(record.read_text())
        self.assertEqual(result['argv'][-3:],['/workspace/src_pipe/Fourier_Quad_Pipe','--config','/fq-config/pipeline.ini'])
        self.assertIn('--no-eval',result['argv'])
        self.assertEqual(result['env']['APPTAINERENV_PMI_RANK'],'5')
        self.assertEqual(result['env']['APPTAINERENV_EXAMPLE'],'a,b=$(false)')
        self.assertNotIn('APPTAINERENV_BAD',result['env'])
        self.assertNotIn('APPTAINER_BIND',result['env'])

    def test_hpc_run_leaves_allocation_checks_to_srun(self):
        p = self.configured()
        env, record = self.install_stubs()
        for key in list(env):
            if key.startswith('SLURM_'):
                del env[key]
        r = subprocess.run(['bash',str(HERE/'runner/run.sh'),str(p['run_dir']/'cpppipeline.env')],
                           capture_output=True,text=True,env=env)
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(json.loads(record.read_text())['argv'][0], str(self.root/'bin/srun'))

    def test_hpc_runner_does_not_interpret_slurm_step_ids(self):
        p = self.configured()
        env, record = self.install_stubs()
        env.update(SLURM_PROCID='0', SLURM_STEP_ID='0')
        r = subprocess.run(['bash',str(HERE/'runner/run.sh'),str(p['run_dir']/'cpppipeline.env')],
                           capture_output=True,text=True,env=env)
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertTrue(record.exists())

    def test_build_does_not_require_slurm(self):
        p = self.configured()
        env, record = self.install_stubs()
        for key in list(env):
            if key.startswith('SLURM_'):
                del env[key]
        r = subprocess.run(['bash',str(HERE/'runner/build.sh'),str(p['run_dir']/'cpppipeline.env')],
                           capture_output=True,text=True,env=env)
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(json.loads(record.read_text())['argv'][0], str(self.root/'bin/apptainer'))

    def test_run_does_not_automatically_run_configuration_checks(self):
        p = self.configured()
        env, record = self.install_stubs()
        self.cfg.write_text('deliberately invalid JSON')
        r = subprocess.run(['bash',str(HERE/'runner/run.sh'),str(p['run_dir']/'cpppipeline.env')],
                           capture_output=True,text=True,env=env)
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertTrue(record.exists())
        check = subprocess.run(['bash',str(HERE/'runner/run.sh'),'--check',str(p['run_dir']/'cpppipeline.env')],
                               capture_output=True,text=True,env=env)
        self.assertNotEqual(check.returncode,0)

    def test_exit_code_preserved(self):
        p = self.configured()
        env, _ = self.install_stubs()
        env['FQ_TEST_EXIT']='17'
        r = subprocess.run(['bash',str(HERE/'runner/run.sh'),str(p['run_dir']/'cpppipeline.env')],capture_output=True,text=True,env=env)
        self.assertEqual(r.returncode,17)

    def test_hpc_build_does_not_force_pmi(self):
        p = self.configured()
        env, record = self.install_stubs()
        r = subprocess.run(['bash',str(HERE/'runner/build.sh'),str(p['run_dir']/'cpppipeline.env')],capture_output=True,text=True,env=env)
        self.assertEqual(r.returncode,0,r.stderr)
        result=json.loads(record.read_text())
        self.assertIn(str(p['source'])+':/workspace/src_pipe:rw',result['argv'])
        self.assertNotIn('APPTAINERENV_OMPI_MCA_ess',result['env'])

    def test_unsupported_source_signature_is_rejected(self):
        self.save()
        (self.source / 'include/RuntimeConfig.hpp').write_text('struct LensingRuntimeConfig { std::string source_cat; };')
        with self.assertRaises(fq.ConfigError): fq.load_project(self.cfg)

    def test_float_underscore_is_rejected(self):
        p = self.load()
        ini = fq.render_ini(p).replace('# pixel_size =\n', 'pixel_size = 1_000\n')
        with self.assertRaises(fq.ConfigError): fq.validate_ini(p, ini)

    def test_symlink_delimiter_rejected(self):
        target = self.root / 'bad:target'
        target.mkdir()
        (self.root / 'link').symlink_to(target, target_is_directory=True)
        self.obj['mounts'][0]['host'] = 'link'
        with self.assertRaises(fq.ConfigError): self.load()

    def test_all_phases_explicitly_disabled(self):
        p = self.load()
        ini = fq.render_ini(p)
        for name in ('astrocat', 'extcat', 'init', 'main', 'rearr', 'fd'):
            ini = ini.replace('# run_process_' + name + ' =\n', 'run_process_' + name + ' = false\n')
        with self.assertRaises(fq.ConfigError): fq.validate_ini(p, ini)

    def test_output_file_requires_parent_directory_mount(self):
        self.obj['mounts'].append({'host':'inputs/list.txt','container':'/lists/output.list','mode':'rw'})
        self.obj['paths']['process.rearranged_expo_list_filename'] = '/lists/output.list'
        with self.assertRaises(fq.ConfigError): self.load()

    def test_batch_environment_is_passed_without_policy_checks(self):
        p = self.configured()
        env, record = self.install_stubs()
        env.update(SLURM_PROCID='0', SLURM_STEP_ID='batch')
        for name in ('run.sh', 'build.sh'):
            with self.subTest(entry=name):
                r = subprocess.run(['bash',str(HERE/'runner'/name),str(p['run_dir']/'cpppipeline.env')],
                                   capture_output=True,text=True,env=env)
                self.assertEqual(r.returncode,0,r.stderr)


if __name__ == '__main__':
    unittest.main()
