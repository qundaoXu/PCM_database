"""Verify imported catalog data, real references, database queries, and SOAP joins."""
import base64
import json
import unittest
from collections import Counter
from pathlib import Path

import numpy as np
from ase.db import connect
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]


class CatalogIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = json.loads((ROOT / "docs/data.json").read_text())
        cls.references = json.loads((ROOT / "docs/aimd_references.json").read_text())

    def test_catalog_counts_models_and_missing_AIMD_errors(self):
        self.assertEqual(len(self.catalog), 16944)
        self.assertEqual(Counter(r["generation_origin"] for r in self.catalog), {"AIMD": 2044, "AI": 14900})
        self.assertTrue(all(r["rdf_peak_error_pct"] is None for r in self.catalog if r["generation_origin"] == "AIMD"))
        self.assertTrue(all(isinstance(r["rdf_peak_error_pct"], float) for r in self.catalog if r["generation_origin"] == "AI"))
        ai = [r for r in self.catalog if r["generation_origin"] == "AI"]
        self.assertEqual(Counter(r["modeling_method"] for r in ai), {
            "MLIP_MACE": 4416, "MLIP_DeePMD": 833, "MLIP_NEP": 4014, "Diffusion": 4404, "Flow Matching": 1233})

    def test_reference_identity_and_real_POSCAR_downloads(self):
        self.assertEqual(len(self.references), 790)
        grouped = {}
        for r in self.references:
            grouped.setdefault(r["reference_id"], set()).add(r["reference_frame_0based"])
        rdf = json.loads((ROOT / "reorg/AIMD_reference/rdf_summary.json").read_text())
        self.assertEqual(len(grouped), 79)
        for trajectory in rdf["trajectories"]:
            frames = grouped[trajectory["reference_id"]]
            samples = trajectory["sampled_indices"]
            self.assertEqual(len(frames), 10)
            self.assertTrue(frames <= set(samples))
            self.assertTrue({samples[0], samples[-1], samples[len(samples) // 2]} <= frames)
        ids = {r["reference_id"] for r in self.references}
        for r in self.catalog + self.references:
            self.assertTrue((ROOT / "docs" / r["structure_url"]).is_file(), r["sid"])
            self.assertTrue(set(r["reference_ids"]) <= ids)
        for r in [self.catalog[0], self.catalog[2044], self.catalog[-1], self.references[0]]:
            atoms = read(ROOT / "docs" / r["structure_url"], format="vasp")
            self.assertEqual(len(atoms), r["natoms"])
            self.assertAlmostEqual(len(atoms) / atoms.get_volume(), r["density"], delta=1e-5)

    def test_joint_descriptors_use_all_catalog_and_reference_ids(self):
        d = json.loads((ROOT / "docs/descriptors.json").read_text())
        self.assertEqual(d["dim"], 252)
        self.assertEqual(set(d["sids"]), {r["sid"] for r in self.catalog + self.references})
        vectors = np.frombuffer(base64.b64decode(d["data"]), dtype="<f4").reshape(-1, 252)
        self.assertEqual(vectors.shape, (17734, 252))
        self.assertTrue(np.isfinite(vectors).all())
        np.testing.assert_allclose(np.linalg.norm(vectors, axis=1), 1, atol=1e-6)

    def test_database_queries_include_AI_and_separate_reference_role(self):
        db = connect(ROOT / "data/pcm.db")
        self.assertEqual(db.count(collection_role="catalog"), 16944)
        self.assertEqual(db.count(generation_origin="AI"), 14900)
        self.assertEqual(db.count(generation_origin="AIMD", collection_role="catalog"), 2044)
        self.assertEqual(db.count(collection_role="reference"), 790)
        ai = db.get(sid=self.catalog[2044]["sid"])
        self.assertAlmostEqual(ai.rdf_peak_error_pct, self.catalog[2044]["rdf_peak_error_pct"])
        self.assertIsNone(db.get(sid="s0000").get("rdf_peak_error_pct"))


if __name__ == "__main__":
    unittest.main()
