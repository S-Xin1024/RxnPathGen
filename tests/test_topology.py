import unittest
from co2rr_network import make_reaction_graph, recognize_reaction_pairs
from co2rr_network.rules import load_rules


class RecognitionTests(unittest.TestCase):
    def test_coupling_from_single_component_inventory(self):
        self.assertEqual(recognize_reaction_pairs({"CO"}, {"COCO"}), {"CO+CO -> COCO"})

    def test_list_and_set_have_same_inventory_semantics(self):
        expected = recognize_reaction_pairs({"CO", "COCO"}, {"COCO", "COCHO"})
        self.assertEqual(recognize_reaction_pairs(["COCO", "CO", "CO"], ["COCHO", "COCO"]), expected)
        self.assertIn("COCO -> COCHO", expected)

    def test_spectator_must_be_available_on_both_sides(self):
        pairs = recognize_reaction_pairs({"CO"}, {"CO", "CHO"})
        self.assertIn("CO+CO -> CO+CHO", pairs)
        self.assertNotIn("CO+CO -> CO+CHO", recognize_reaction_pairs({"CO"}, {"CHO"}))

    def test_empty_unknown_and_invalid_components(self):
        self.assertEqual(recognize_reaction_pairs(set(), {"COCO"}), set())
        self.assertEqual(recognize_reaction_pairs({"unknown"}, {"COCO"}), set())
        for values in ("CO", ["CO+CO"], [""], [1]):
            with self.assertRaises(ValueError):
                recognize_reaction_pairs(values, {"COCO"})

    def test_custom_rules_extend_supported_components(self):
        rules = {"reactions": [
            {"source": ["C2H4"], "target": ["C2H6"], "kind": "hydrogenation"}]}
        self.assertEqual(recognize_reaction_pairs({"C2H4"}, {"C2H6"}, rules=rules),
                         {"C2H4 -> C2H6"})

    def test_full_catalog_reconstructed_from_component_pools(self):
        for profile in ("direct", "out"):
            rules = load_rules(profile)["reactions"]
            source = {c for r in rules for c in r["source"]}
            target = {c for r in rules for c in r["target"]}
            expected = {"+".join(r["source"]) + " -> " + "+".join(r["target"]) for r in rules}
            self.assertEqual(recognize_reaction_pairs(source, target, profile=profile), expected)


class COOGraphTests(unittest.TestCase):
    def test_requested_single_tuple_return(self):
        result = make_reaction_graph({"CO+CO -> COCO", "COCO -> COCHO"})
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], [[0, 2], [2, 1]])

    def test_labels_decode_edges_and_preserve_stoichiometry(self):
        pairs = {"CO+CO -> COCO", "COCO -> COCHO"}
        indices, labels = make_reaction_graph(pairs, return_labels=True)
        decoded = {labels[a] + " -> " + labels[b] for a, b in zip(*indices)}
        self.assertEqual(decoded, pairs)
        self.assertEqual(labels, ("CO+CO", "COCHO", "COCO"))

    def test_order_is_deterministic_and_duplicates_normalized(self):
        a = ["COCO -> COCHO", "CO+CO -> COCO"]
        b = list(reversed(a)) + [" CO + CO → COCO "]
        self.assertEqual(make_reaction_graph(a, return_labels=True), make_reaction_graph(b, return_labels=True))

    def test_empty_and_self_loop(self):
        self.assertEqual(make_reaction_graph(set()), ([[], []],))
        self.assertEqual(make_reaction_graph(set(), return_labels=True), ([[], []], ()))
        self.assertEqual(make_reaction_graph({"CO -> CO"}), ([[0], [0]],))

    def test_direction_and_ordered_sites_preserved(self):
        indices, labels = make_reaction_graph({"CO+CHO -> COCHO", "CHO+CO -> COCHO"}, return_labels=True)
        self.assertEqual(len(labels), 3)
        self.assertEqual(len(indices[0]), 2)
        self.assertNotIn("COCHO -> CO+CHO", {labels[a] + " -> " + labels[b] for a, b in zip(*indices)})

    def test_bad_pair_input_is_rejected(self):
        for pairs in ("CO -> COH", {"CO"}, {"CO -> -> COH"}, {" -> CO"}, {"CO+ -> COH"}, {1}):
            with self.assertRaises(ValueError):
                make_reaction_graph(pairs)
