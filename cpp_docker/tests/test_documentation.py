"""User documentation and template contracts; no scientific run is performed."""
from pathlib import Path
import re
import sys
import unittest
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / 'doc'
sys.path.insert(0, str(ROOT / 'cpp_docker'))
from schema import sections, path_keys


def headings(text):
    found = set()
    for heading in re.findall(r'^#{1,6}\s+(.+)$', text, re.M):
        heading = heading.lower().replace('`', '')
        heading = re.sub(r'[^\w\-\s]', '', heading)
        found.add(re.sub(r'\s', '-', heading))
    return found


def all_markdown_files():
    return sorted(list(ROOT.glob('*.md')) + list(DOC.glob('*.md')))


class DocumentationTests(unittest.TestCase):
    def test_only_readme_at_root_and_others_in_doc(self):
        root_mds = [p.name for p in ROOT.glob('*.md')]
        self.assertEqual(root_mds, ['README.md'])
        for name in ('CPP_GUIDE.md', 'DATA_PREPARATION.md', 'CPP_PIPELINE_PARAMETERS.md',
                     'DOCKER.md', 'HPC.md', 'PATHS.md', 'THIRD_PARTY_NOTICES.md'):
            self.assertTrue((DOC / name).is_file(), name)

    def test_no_mermaid_code_blocks_and_no_lite_references_in_docs(self):
        for p in all_markdown_files():
            text = p.read_text(encoding='utf-8')
            self.assertNotIn('```mermaid', text, str(p))
            self.assertNotIn('cpp_Lite', text, str(p))
            self.assertNotIn('Initialize.hpp', text, str(p))

    def test_all_runtime_keys_in_full_reference(self):
        ref = (DOC / 'CPP_PIPELINE_PARAMETERS.md').read_text(encoding='utf-8')
        for section, keys in sections().items():
            for key in keys:
                self.assertIn('`[%s].%s`' % (section, key), ref)

    def test_all_container_path_keys_documented(self):
        paths = (DOC / 'PATHS.md').read_text(encoding='utf-8')
        for key in path_keys():
            self.assertIn('`' + key + '`', paths)

    def test_native_templates_have_no_active_defaults(self):
        text = (ROOT / 'cpp_Standard' / 'pipeline.example.ini').read_text(encoding='utf-8')
        self.assertFalse(re.search(r'^\s*[^#;\[\s][^\n]*=', text, re.M))
        for section, keys in sections().items():
            for key in keys:
                self.assertIn('# ' + key + ' =\n', text)

    def test_user_documents_have_no_cjk(self):
        for p in all_markdown_files():
            self.assertFalse(re.search(r'[\u3400-\u9fff]', p.read_text(encoding='utf-8')), str(p))

    def test_root_navigation_prioritizes_cpp(self):
        readme = (ROOT / 'README.md').read_text(encoding='utf-8')
        for name in ('doc/CPP_GUIDE.md', 'doc/DATA_PREPARATION.md', 'doc/CPP_PIPELINE_PARAMETERS.md',
                     'doc/DOCKER.md', 'doc/HPC.md', 'doc/PATHS.md', 'LICENSE', 'doc/THIRD_PARTY_NOTICES.md'):
            self.assertIn('(' + name, readme)
        self.assertNotIn('Develop.md', readme)

    def test_local_document_links_and_anchors(self):
        for p in all_markdown_files():
            text = p.read_text(encoding='utf-8')
            for target in re.findall(r'\]\(([^\s)]+)\)', text):
                if '://' in target or target.startswith('mailto:'):
                    continue
                name, _, anchor = unquote(target).partition('#')
                q = (p.parent / name).resolve() if name else p
                if name.endswith('LICENSE') and not q.exists():
                    continue
                self.assertTrue(q.exists(), '%s -> %s' % (p.name, target))
                if anchor and q.suffix == '.md':
                    self.assertIn(anchor, headings(q.read_text(encoding='utf-8')),
                                  '%s -> %s' % (p.name, target))

    def test_data_preparation_contract_retained(self):
        data = (DOC / 'DATA_PREPARATION.md').read_text(encoding='utf-8')
        for required in ('Science images', 'Gaia catalog', 'External source catalog',
                         'DQ masks', 'gaia_p1_00.cat', 'gaia_m3_12.cat', 'gaia_p9.cat',
                         'astra_RA_123_124_Dec_m05_m04.dat',
                         'extern_RA_123_124_Dec_m05_m04.dat', 'CCDNUM',
                         'expand the upper and lower Dec', '2°, 4°,', 'header line'):
            self.assertIn(required, data)

    def test_full_header_parameter_inventory(self):
        files = sorted(ROOT.glob('cpp_Standard/config/*.hpp'))
        if not files:
            self.skipTest('A replacement-only package does not include unchanged configuration headers')
        ref = (DOC / 'CPP_PIPELINE_PARAMETERS.md').read_text(encoding='utf-8')
        for path in files:
            source = re.sub(r'//[^\n]*|/\*.*?\*/', '', path.read_text(), flags=re.S)
            names = re.findall(r'\b(?:constexpr|const)\s+[^;{}=]*?\b([A-Za-z_]\w*)\s*(?:\[\s*\])?\s*=', source)
            for name in names:
                self.assertRegex(ref, r'`(?:[A-Za-z_]\w*::)?' + re.escape(name) + r'`',
                                 '%s: missing parameter %s' % (path, name))

    def test_runner_has_no_scheduler_guard(self):
        for path in (ROOT / 'cpp_docker/runner').glob('*.sh'):
            source = path.read_text()
            self.assertNotIn('SLURM_JOB_ID', source)
            self.assertNotIn('SLURM_STEP_ID', source)
            self.assertNotIn('SLURM_NTASKS', source)
        run = (ROOT / 'cpp_docker/runner/run.sh').read_text()
        self.assertEqual(run.count('fq_check'), 1)
        self.assertIn('== --check', run)


if __name__ == '__main__':
    unittest.main()
