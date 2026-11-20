"""One-time reproducible extraction; never imported by the installed package.

Only reaction_pair's AST is executed. Original top-level jobs are not executed.
Dictionary assignment history recovers rules lost through legacy numeric key collisions.
"""
import ast
import copy
import hashlib
import json
import re
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = PACKAGE / "legacy" if (PACKAGE / "legacy").is_dir() else PACKAGE.parent


class HistoryDict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.history = []

    def __setitem__(self, key, value):
        self.history.append((key, value))
        super().__setitem__(key, value)

    def update(self, other):
        for k, v in other.items():
            self[k] = v


def extract(filename, profile):
    path = ROOT / filename
    source = path.read_text(encoding="utf-8")
    function = next(n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == "reaction_pair")
    module = ast.parse("")
    module.body = [function]
    namespace = {"copy": copy, "re": re, "dict": HistoryDict}
    exec(compile(module, str(path), "exec"), namespace)
    families = ["cc_coupling", "deoxygenation", "oxygen_hydrogenation", "carbon_hydrogenation"]
    if profile == "direct":
        families.append("desorption")
    families += ["rearrangement", "oxygen_fragment"]
    rules, seen, recovered = [], set(), 0
    for kind, mapping in zip(families, namespace["reaction_pair"]()):
        current = set(mapping.values())
        for key, (left, right) in mapping.history:
            signature = (tuple(left), tuple(right), kind)
            if signature in seen or left == right:
                continue
            seen.add(signature)
            lost = (left, right) not in current
            recovered += int(lost)
            rules.append({"source": list(left), "target": list(right), "kind": kind,
                          "metadata": {"legacy_file": filename, "recovered_key_collision": lost}})
    rules.sort(key=lambda r: (r["kind"], r["source"], r["target"]))
    data = {"schema_version": 1, "profile": profile,
            "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "recovered_key_collision_rules": recovered,
            "formulas": {}, "reactions": rules}
    destination = PACKAGE / "src/co2rr_network/resources" / (profile + "_rules.json")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {"profile": profile, "rules": len(rules), "recovered_key_collision_rules": recovered}


if __name__ == "__main__":
    print(extract("Judgment_reaction_pair_direct_des_2.py", "direct"))
    print(extract("Judgment_reaction_pair_out_des.py", "out"))
