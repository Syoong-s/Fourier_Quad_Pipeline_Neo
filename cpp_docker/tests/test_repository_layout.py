"""Guard against stale build and CI paths after repository reorganizations."""

from pathlib import Path
import json
import re
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
CPP = ROOT / 'cpp_Standard'
CI = ROOT / '.github' / 'workflows'

sys.path.insert(0, str(ROOT / 'cpp_docker'))
from schema import sections


class RepositoryLayoutTests(unittest.TestCase):
    def test_supported_pipeline_is_present(self):
        for filename in ('Makefile', 'main.cpp', 'src/RuntimeConfig.cpp',
                         'include/RuntimeConfig.hpp', 'pipeline.example.ini',
                         'tests/RuntimeConfigTest.cpp'):
            self.assertTrue((CPP / filename).is_file(), filename)
        self.assertFalse((ROOT / 'cpp_Lite').exists(),
                         'The removed Lite tree must not be required by active builds')

    def test_every_makefile_source_exists(self):
        makefile = (CPP / 'Makefile').read_text(encoding='utf-8')
        match = re.search(r'^SRCS\s*=\s*(.*?)^OBJS\s*=', makefile, re.M | re.S)
        self.assertIsNotNone(match, 'Cannot identify Makefile SRCS list')
        sources = re.findall(r'[\w./-]+\.cpp', match.group(1))
        self.assertGreater(len(sources), 1, 'Empty or incomplete source list')
        self.assertEqual(len(sources), len(set(sources)), 'Duplicate Makefile source')
        for name in sources:
            with self.subTest(source=name):
                self.assertTrue((CPP / name).is_file(), name)

    def test_workflows_do_not_reference_removed_lite_tree(self):
        self.assertTrue(CI.is_dir())
        workflows = sorted(CI.glob('*.yml')) + sorted(CI.glob('*.yaml'))
        self.assertTrue(workflows, 'No GitHub workflows found')
        for workflow in workflows:
            with self.subTest(workflow=workflow.name):
                self.assertNotIn('cpp_Lite', workflow.read_text(encoding='utf-8'))
        config_workflow = (CI / 'configuration.yml').read_text(encoding='utf-8')
        for filename in ('cpp_Standard/src/RuntimeConfig.cpp',
                         'cpp_Standard/tests/RuntimeConfigTest.cpp',
                         'cpp_Standard/pipeline.example.ini'):
            self.assertIn(filename, config_workflow)

    def test_json_examples_use_existing_source_tree(self):
        for name in ('project.example.json', 'project.full.example.json'):
            path = ROOT / 'cpp_docker' / name
            config = json.loads(path.read_text(encoding='utf-8'))
            self.assertEqual(config.get('variant', 'standard'), 'standard')
            source = (path.parent / config['source']).resolve()
            self.assertEqual(source, CPP.resolve())
            self.assertTrue((source / 'Makefile').is_file())

    def test_launcher_supports_only_existing_variant(self):
        self.assertTrue(sections('standard'))
        with self.assertRaises(ValueError):
            sections('lite')


if __name__ == '__main__':
    unittest.main()
