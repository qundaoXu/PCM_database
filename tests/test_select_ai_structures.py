"""Scientific and selection invariants for the AI structure collection."""
import unittest

import numpy as np
from ase import Atoms

from scripts import select_ai_structures as selection


class SelectionTests(unittest.TestCase):
    def test_capped_budget_preserves_rare_strata(self):
        counts = np.array([30001, 20001, 150, 9, 1])
        quotas = selection.allocate_budgets(counts, 1000)
        self.assertEqual(int(quotas.sum()), 1000)
        self.assertTrue(np.all((quotas >= 1) & (quotas <= counts)))
        self.assertEqual(quotas[-1], 1)
        np.testing.assert_array_equal(selection.allocate_budgets([1, 2], 10), [1, 2])

    def test_time_orders_and_missing_information(self):
        values = [None, 0, 1, 9.999, 10, 100, 1000, 10000]
        expected = ["none", "0ps", "1ps", "1ps", "10ps", "100ps", "1ns", "10ns"]
        self.assertEqual([selection.time_order_label(v) for v in values], expected)
        row = dict(file="GST225/example.dump", temperature_K=373, time_ps=None)
        normalized = selection.normalize_candidate(row)
        self.assertEqual(normalized["modeling_method"], "MLIP_MACE")
        self.assertEqual(normalized["temperature_K"], 300)
        self.assertEqual(normalized["reported_temperature_K"], 373)
        self.assertEqual(normalized["simulation_time_order"], "none")
        self.assertEqual(selection.header_time_fs(b'2\nTime=12000.00000000 pbc="T T T"\nGe 0 0 0\nTe 2 0 0\n'), 12000)
        self.assertIsNone(selection.header_time_fs(b'2\nLattice="10 0 0 0 10 0 0 0 10"\n'))

    def test_nonuniform_storage_is_sampled_by_physical_time(self):
        selected = selection.uniform_time_sample([0, 1, 2, 3, 4, 5], [0, 1, 2, 10, 90, 100], 3)
        self.assertEqual(selected, [0, 3, 5])
        self.assertEqual(len(set(selected)), 3)
        all_items = selection.uniform_time_sample([0, 1], [0, 100], 10)
        self.assertEqual(all_items, [0, 1])

    def test_reference_rounding_is_explicit_and_species_strict(self):
        refs = {"Sb:2;Te:1": {"material": "Sb2Te"}}
        ref, match, deviation = selection.reference_match("Sb:667;Te:333", refs, 0.001)
        self.assertEqual(ref["material"], "Sb2Te")
        self.assertEqual(match, "integer_stoichiometry_rounding")
        self.assertAlmostEqual(deviation, 0.0333333333333)
        self.assertIsNone(selection.reference_match("Sb:1;Te:1", refs, 0.001)[0])
        self.assertIsNone(selection.reference_match("Ge:2;Te:1", refs, 0.001)[0])

    def test_descriptor_is_periodic_and_atom_order_invariant(self):
        atoms = Atoms("Ge2Te2", positions=[[0.1, 1, 1], [11.9, 1, 1], [5, 5, 5], [7.6, 5, 5]],
                      cell=np.eye(3) * 12, pbc=True)
        a, g, _, _ = selection.structural_features(atoms)
        shifted = atoms[[3, 0, 2, 1]]
        shifted.positions += [5, 7, 12]
        b, h, _, _ = selection.structural_features(shifted)
        np.testing.assert_allclose(a, b, atol=1e-5)
        np.testing.assert_allclose(g, h, atol=1e-5)

    def test_farthest_point_selection_is_deterministic_and_covers_extremes(self):
        features = np.column_stack([np.linspace(0, 1, 101), np.zeros(101)]).astype(np.float32)
        selected, _ = selection.diversity_sample(features, 3)
        self.assertEqual(set(selected), {0, 50, 100})
        self.assertEqual(selection.diversity_sample(features, 3)[0], selected)

    def test_high_dimensional_diversity_retains_95_percent_variance(self):
        features = np.random.default_rng(17).normal(size=(300, 60)).astype(np.float32)
        selected, coverage = selection.diversity_sample(features, 10)
        self.assertEqual(len(set(selected)), 10)
        self.assertGreaterEqual(coverage["retained_variance"], 0.95)
        self.assertGreater(coverage["pca_dimensions"], 32)


if __name__ == "__main__":
    unittest.main()
