import configparser
import ast
import os
import re
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class PackageBoundaryTests(unittest.TestCase):
    def test_core_import_and_path_search_without_computational_libraries(self):
        script = '''import sys
class Block:
    def find_spec(self, name, path=None, target=None):
        if name.split('.')[0] in {'torch', 'ase', 'numpy', 'scipy', 'matplotlib', 'fairchem', 'BUCToolkit'} or name.startswith('co2rr_network._vendor'):
            raise AssertionError('Unexpected core import: ' + name)
sys.meta_path.insert(0, Block())
from co2rr_network import build_network, minimax_path
g = build_network({'CO+CO -> COCO': -.4, 'COCO -> COCHO': -.2})
assert minimax_path(g, 'CO+CO', 'COCHO')['cost'] == -.2
'''
        env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), PYTHONDONTWRITEBYTECODE="1")
        result = subprocess.run([sys.executable, "-B", "-c", script], env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(result.returncode, 0, result.stderr.decode(errors="replace"))

    def test_no_computational_adapter_commands(self):
        from co2rr_network.cli import parser
        choices = next(a.choices for a in parser()._actions if hasattr(a, "choices") and isinstance(a.choices, dict))
        self.assertEqual(set(choices), {"network", "rules", "import-legacy", "path", "paths", "prune", "plot"})

    def test_no_mandatory_third_party_runtime_dependency(self):
        config = configparser.ConfigParser()
        config.read(str(ROOT / "setup.cfg"), encoding="utf-8")
        self.assertFalse(config.get("options", "install_requires", fallback="").strip())
        self.assertEqual(set(config["options.extras_require"]), {"plots", "kinetics"})

    def test_vendored_license_and_internal_import_prefix(self):
        vendor = ROOT / "src/co2rr_network/_vendor/fairchem"
        self.assertEqual((vendor / "LICENSE").read_bytes(), (vendor / "core/LICENSE.md").read_bytes())
        self.assertIn("MIT License", (vendor / "LICENSE").read_text(encoding="utf-8"))
        self.assertTrue((vendor / "core/models/equiformer_v2/LICENSE").is_file())
        for path in vendor.rglob("*.py"):
            self.assertIsNone(re.search(r"^\s*(?:from|import) fairchem\b", path.read_text(encoding="utf-8"), re.M), str(path))
        source = (vendor / "core/common/utils.py").read_text(encoding="utf-8")
        self.assertIn('module_name = "co2rr_network._vendor." + module_name', source)

    def test_legacy_qualified_model_name_resolves_to_internal_source(self):
        from types import SimpleNamespace
        source = (ROOT / "src/co2rr_network/_vendor/fairchem/core/common/registry.py").read_text(encoding="utf-8")
        tree = ast.parse(source)
        function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "_get_absolute_mapping")
        module = ast.parse("")
        module.body = [function]
        requested = []
        marker = object()
        def importer(name):
            requested.append(name)
            return SimpleNamespace(GemNetOC=marker)
        namespace = {"importlib": SimpleNamespace(import_module=importer)}
        exec(compile(module, "local_registry", "exec"), namespace)
        self.assertIs(namespace["_get_absolute_mapping"]("fairchem.core.models.gemnet_oc.GemNetOC"), marker)
        self.assertEqual(requested, ["co2rr_network._vendor.fairchem.core.models.gemnet_oc"])
