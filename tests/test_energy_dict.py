import unittest
from co2rr_network import ReactionNetwork, best_path, build_network, minimax_path


RULES = {"reactions": [
    {"source": ["CO", "CO"], "target": ["COCO"], "kind": "cc_coupling"},
    {"source": ["COCO"], "target": ["COCHO"], "kind": "hydrogenation"}]}


class EnergyDictionaryTests(unittest.TestCase):
    def test_state_formula_to_scalar(self):
        g = build_network({"CO": -.1, "COCO": .2, "COCHO": 0}, rules=RULES)
        self.assertAlmostEqual(g.reactions[0].delta_g, .4)
        self.assertAlmostEqual(g.reactions[1].delta_g, -.2)

    def test_reaction_formula_to_delta_g_scalar(self):
        g = build_network({"CO+CO -> COCO": -.4, "COCO -> COCHO": -.2})
        self.assertEqual(minimax_path(g, "CO+CO", "COCHO")["cost"], -.2)
        self.assertTrue(all(s.free_energy is None for s in g.states))

    def test_reaction_formula_to_barrier_scalar(self):
        g = build_network({"CO+CO -> COCO": .6, "COCO -> COCHO": .2}, energy_type="barrier")
        self.assertEqual(best_path(g, "CO+CO", "COCHO", mode="kinetic")["cost"], .6)
        self.assertTrue(all(r.delta_g is None for r in g.reactions))
        with self.assertRaises(ValueError):
            best_path(g, "CO+CO", "COCHO")

    def test_both_delta_g_and_barrier(self):
        g = build_network({"CO+CO → COCO": {"delta_g": .4, "barrier": .6}})
        self.assertEqual(g.reactions[0].delta_g, .4)
        self.assertEqual(g.reactions[0].barrier, .6)

    def test_target_formula_when_predecessor_is_unique(self):
        g = build_network({"COCO": {"delta_g": .4, "barrier": .6}}, rules=RULES)
        self.assertEqual(len(g.reactions), 1)
        self.assertEqual(g.reactions[0].source, "CO+CO")

    def test_ambiguous_target_formula_is_rejected(self):
        rules = {"reactions": RULES["reactions"] + [
            {"source": ["CO", "CHO"], "target": ["COCHO"], "kind": "cc_coupling"}]}
        with self.assertRaisesRegex(ValueError, "predecessors"):
            build_network({"COCHO": {"barrier": .6}}, rules=rules)

    def test_missing_quantities_stay_unknown(self):
        g = build_network({"CO+CO -> COCO": {"barrier": .5}})
        self.assertIsNone(g.reactions[0].delta_g)
        self.assertTrue(all(s.free_energy is None for s in g.states))
        self.assertIsInstance(g, ReactionNetwork)

    def test_empty_bad_and_nonfinite_inputs(self):
        for energy in ({}, "file.csv", {"CO": float("nan")}, {"COCO": {"energy_unknown": .1}},
                       {"CO+CO -> COCO": {"barrier": -.1}}):
            with self.assertRaises(ValueError):
                build_network(energy)

    def test_reverse_requires_reaction_energy(self):
        with self.assertRaises(ValueError):
            build_network({"CO+CO -> COCO": {"barrier": .5}}, include_reverse=True)
        g = build_network({"CO+CO -> COCO": {"delta_g": .4, "barrier": .5}}, include_reverse=True)
        self.assertAlmostEqual(g.reactions[1].barrier, .1)
