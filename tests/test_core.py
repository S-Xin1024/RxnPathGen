import csv
import json
import math
import tempfile
import unittest
from pathlib import Path

from co2rr_network.analysis import best_path, path_statistics, prune_network
from co2rr_network.chemistry import KB_EV, balanced_step, carbon_valence, che_shift, composition
from co2rr_network.data import import_legacy_csv, load_energies
from co2rr_network.io import write_csv, write_json
from co2rr_network.network import build_network
from co2rr_network.rules import enumerate_states, load_rules
from co2rr_network.schema import Reaction, ReactionNetwork, State
from co2rr_network.workflow import run_network


def small_network():
    nodes = [State("CO+CO", ("CO", "CO"), 0), State("COCO", ("COCO",), .4),
             State("COCHO", ("COCHO",), .1), State("CHOCHO_d", ("CHOCHO_d",), -.2, True)]
    edges = [Reaction("a", "CO+CO", "COCO", "cc", .4, .6),
             Reaction("b", "COCO", "COCHO", "h", -.3, .2),
             Reaction("c", "COCHO", "CHOCHO_d", "h", -.3, .3),
             Reaction("d", "CO+CO", "COCHO", "other", .1, .8)]
    return ReactionNetwork(nodes, edges)


class ChemistryTests(unittest.TestCase):
    def test_formula_multidigit_and_legacy_fragments(self):
        self.assertEqual(composition("CH12O_CH3CH3_d"), {"C": 3, "H": 18, "O": 1})
        self.assertEqual(composition("COOHCHO"), {"C": 2, "H": 2, "O": 3})
        self.assertEqual(carbon_valence("CO"), 2)

    def test_invalid_formula(self):
        for name in ("", "CO?", "C0H", "CuCO"):
            with self.assertRaises(ValueError):
                composition(name)

    def test_balance(self):
        self.assertEqual(balanced_step(["COOH"], ["CO"]), {
            "electrons_consumed": 1, "protons_consumed": 1, "water_produced": 1})
        self.assertEqual(balanced_step(["CO", "CHO"], ["COCHO"])["electrons_consumed"], 0)
        with self.assertRaises(ValueError):
            balanced_step(["CO"], ["COCO"])

    def test_che_sign_and_ph_scale(self):
        counts = composition("CO")
        self.assertAlmostEqual(che_shift(counts, -.5), -1)
        self.assertEqual(che_shift(counts, 0, 14, scale="RHE"), 0)
        expected = 2 * KB_EV * 298.15 * math.log(10) * 14
        self.assertAlmostEqual(che_shift(counts, 0, 14, scale="SHE"), expected)
        with self.assertRaises(ValueError):
            che_shift(counts, temperature=0)

class PathTests(unittest.TestCase):
    def test_thermodynamic_and_kinetic_are_distinct(self):
        graph = small_network()
        result = best_path(graph, "CO+CO", "CHOCHO_d")
        self.assertAlmostEqual(result["cost"], .1)
        self.assertEqual(result["reactions"], ["d", "c"])
        result = best_path(graph, "CO+CO", "CHOCHO_d", "kinetic")
        self.assertAlmostEqual(result["cost"], .6)
        self.assertEqual(result["reactions"], ["a", "b", "c"])

    def test_sum_and_unreachable(self):
        graph = small_network()
        self.assertAlmostEqual(best_path(graph, "CO+CO", "CHOCHO_d", objective="sum")["cost"], .1)
        self.assertFalse(best_path(graph, "CHOCHO_d", "CO+CO")["reachable"])
        self.assertEqual(best_path(graph, "COCO", "COCO")["cost"], 0)

    def test_negative_steps_and_missing_barrier(self):
        graph = small_network()
        graph.reactions[0].barrier = None
        with self.assertRaises(ValueError):
            best_path(graph, "CO+CO", "COCO", "kinetic")
        self.assertGreaterEqual(best_path(graph, "COCO", "CHOCHO_d")["cost"], 0)

    def test_dag_counts_and_occurrence(self):
        result = path_statistics(small_network(), "CO+CO", "CHOCHO_d")
        self.assertTrue(result["exact"])
        self.assertEqual(result["total_paths"], 2)
        self.assertEqual(result["node_counts"]["COCO"], 1)
        self.assertEqual(result["node_counts"]["COCHO"], 2)

    def test_cycles_and_bounds(self):
        graph = small_network()
        graph.reactions.append(Reaction("back", "COCHO", "COCO", "back", .3, .5))
        full = path_statistics(graph, "CO+CO", "CHOCHO_d")
        self.assertEqual(full["total_paths"], 2)
        self.assertTrue(full["exact"])
        limited = path_statistics(graph, "CO+CO", "CHOCHO_d", max_paths=1)
        self.assertFalse(limited["exact"])
        depth = path_statistics(graph, "CO+CO", "CHOCHO_d", max_depth=1)
        self.assertFalse(depth["exact"])

    def test_prune_and_no_graph_mutation(self):
        graph = small_network()
        pruned = prune_network(graph, "CO+CO", ["CHOCHO_d"], .2)
        self.assertEqual({r.id for r in pruned.reactions}, {"d", "c"})
        self.assertEqual(len(graph.reactions), 4)

    def test_dag_large_integer_counts(self):
        # 2 parallel channels at each of 70 steps: no int64 overflow.
        nodes = [State("C" + "H" * i, ("C" + "H" * i,), 0) for i in range(71)]
        reactions = [Reaction("{}_{}".format(i, k), nodes[i].id, nodes[i+1].id, "channel", 0)
                     for i in range(70) for k in range(2)]
        self.assertEqual(path_statistics(ReactionNetwork(nodes, reactions), nodes[0].id, nodes[-1].id)["total_paths"], 2 ** 70)


class FileTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_network_round_trip_and_schema(self):
        graph = small_network()
        graph.save(self.root / "network.json")
        self.assertEqual(ReactionNetwork.load(self.root / "network.json"), graph)
        with self.assertRaises(ValueError):
            ReactionNetwork(graph.states, [Reaction("bad", "COCO", "COCHO", "h", 0)])

    def test_legacy_import_and_composite_name(self):
        source = self.root / "old.csv"
        source.write_text("Ga_O_CH4,-1.2\nGa_CO,-.4\n", encoding="utf-8")
        target = self.root / "new.csv"
        self.assertEqual(import_legacy_csv(source, target, facet="111"), 2)
        records, _ = load_energies(target, "Ga", "111")
        self.assertIn("O_CH4", records)

    def test_lowest_conformer_and_facet_isolation(self):
        path = self.root / "energies.csv"
        write_csv(path, [
            {"catalyst": "Cu", "facet": "100", "species": "CO", "free_energy": -1},
            {"catalyst": "Cu", "facet": "111", "species": "CO", "free_energy": -.5},
            {"catalyst": "Cu", "facet": "111", "species": "CO", "free_energy": -.8},
            {"catalyst": "Cu", "facet": "111", "species": "CO", "free_energy": -10, "converged": "false"}
        ], ["catalyst", "facet", "species", "free_energy", "converged"])
        records, skipped = load_energies(path, "Cu", "111")
        self.assertEqual(records["CO"]["free_energy"], -.8)
        self.assertEqual(skipped, 1)

    def fixture(self):
        rules = {"reactions": [
            {"source": ["CO", "CO"], "target": ["COCO"], "kind": "cc_coupling"},
            {"source": ["COCO"], "target": ["COCHO"], "kind": "carbon_hydrogenation"}]}
        write_json(self.root / "rules.json", rules)
        rows = [{"catalyst": "Cu", "species": name, "free_energy": g}
                for name, g in [("CO", -.1), ("COCO", .2), ("COCHO", 0)]]
        write_csv(self.root / "energies.csv", rows, ["catalyst", "species", "free_energy"])
        return {"energies": {name: g for name, g in [("CO", -.1), ("COCO", .2), ("COCHO", 0)]}, "rules": rules}

    def test_assembly_reference_and_che(self):
        kwargs = self.fixture()
        graph = build_network(**kwargs)
        self.assertAlmostEqual(graph.reactions[0].delta_g, .4)
        shifted = build_network(**kwargs, potential=-.5)
        self.assertAlmostEqual(shifted.reactions[0].delta_g, .4)  # C-C coupling is not PCET
        self.assertAlmostEqual(shifted.reactions[1].delta_g, -.7)

    def test_explicit_coadsorbed_energy(self):
        kwargs = self.fixture()
        kwargs["energies"]["CO+CO"] = -.4
        graph = build_network(**kwargs)
        self.assertAlmostEqual(graph.reactions[0].delta_g, .6)
        self.assertEqual(graph.states[0].metadata["energy_model"], "explicit_state")

    def test_missing_energy_is_explicit(self):
        kwargs = self.fixture()
        kwargs.pop("rules")
        with self.assertRaises(ValueError):
            build_network(**kwargs, missing_energy="error")
        graph = build_network(**kwargs, missing_energy="skip")
        self.assertGreater(len(graph.metadata["missing_state_energies"]), 0)

    def test_barrier_validation_filter_and_reverse(self):
        kwargs = self.fixture()
        kwargs["energies"]["CO+CO -> COCO"] = {"barrier": .5}
        graph = build_network(**kwargs, include_reverse=True)
        reverse = next(r for r in graph.reactions if r.metadata.get("reverse_of") == "r00000")
        self.assertAlmostEqual(reverse.barrier, .1)
        filtered = prune_network(graph, "CO+CO", ["COCHO"], threshold=.3)
        self.assertEqual(len(filtered.reactions), 0)
        kwargs["energies"]["CO+CO -> COCO"] = {"barrier": .1}
        with self.assertRaises(ValueError):
            build_network(**kwargs)

    def test_supplied_delta_g_agrees_with_state_energies(self):
        kwargs = self.fixture()
        kwargs["energies"]["CO+CO -> COCO"] = {"delta_g": .4, "barrier": .5}
        self.assertAlmostEqual(build_network(**kwargs).reactions[0].delta_g, .4)
        kwargs["energies"]["CO+CO -> COCO"]["delta_g"] = .1
        with self.assertRaises(ValueError):
            build_network(**kwargs)

    def test_network_report_export(self):
        self.fixture()
        config = {"energies": "energies.csv", "rules_path": "rules.json", "output": "out",
                  "network": {"catalyst": "Cu"}, "analysis": {"source": "CO+CO", "targets": ["COCHO"]}}
        write_json(self.root / "config.json", config)
        result = run_network(self.root / "config.json")
        self.assertEqual(result["states"], 3)
        self.assertTrue((self.root / "out/network.dot").is_file())


class CatalogTests(unittest.TestCase):
    def test_builtins_are_balanced_and_deduplicated(self):
        for name, n, r in [("direct", 629, 1102), ("out", 502, 957)]:
            rules = load_rules(name)
            self.assertEqual(len(enumerate_states(rules)), n)
            self.assertEqual(len(rules["reactions"]), r)
            self.assertTrue(all("balance" in rule for rule in rules["reactions"]))
        self.assertEqual(load_rules("out")["recovered_key_collision_rules"], 2)


if __name__ == "__main__":
    unittest.main()
