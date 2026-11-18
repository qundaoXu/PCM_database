"""Integrity checks for archive operations and periodic first-shell RDF."""
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from ase import Atoms
from ase.io import write

from scripts import organize_ai_sources as audit


class SourceAuditTests(unittest.TestCase):
    def test_periodic_rdf_and_finite_n_normalization(self):
        cell = np.eye(3) * 10
        positions = np.array([[0.1, 1, 1], [9.9, 1, 1], [5, 5, 5]])
        curve = audit.rdf(positions, cell)
        np.testing.assert_allclose(curve, audit.rdf(positions + [5, 7, 12], cell))
        self.assertAlmostEqual(float(np.sum(curve * 3 * audit.SHELLS / 1000)), 1)
        with self.assertRaises(ValueError):
            audit.rdf(positions, np.eye(3) * 6)

    def test_geometry_dedup_ignores_atom_order_but_preserves_species(self):
        atoms = Atoms("GeTe", positions=[[0, 0, 0], [2, 0, 0]], cell=np.eye(3) * 10)
        self.assertEqual(audit.geometry_id(atoms), audit.geometry_id(atoms[[1, 0]]))
        different = atoms.copy()
        different.set_chemical_symbols(["Te", "Ge"])
        self.assertNotEqual(audit.geometry_id(atoms), audit.geometry_id(different))

    def test_snapshot_timesteps_do_not_imply_a_chronological_trajectory(self):
        meta = audit.timeline_metadata([270, 627, 288])
        self.assertEqual(meta["timestep_order"], "non_monotonic_or_reset")
        self.assertEqual((meta["min_timestep"], meta["max_timestep"]), (270, 627))
        self.assertEqual(audit.timeline_metadata([0, 100, 200])["timestep_order"], "strictly_increasing")

    def test_xdatcar_rejects_incomplete_last_frame(self):
        atoms = Atoms("GeTe", scaled_positions=[[0, 0, 0], [0.2, 0, 0]], cell=np.eye(3) * 10)
        buffer = io.StringIO()
        write(buffer, [atoms, atoms], format="vasp-xdatcar")
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "XDATCAR"
            path.write_text(buffer.getvalue())
            self.assertEqual(len(list(audit.xdat_frames(path))), 2)
            path.write_text("\n".join(buffer.getvalue().splitlines()[:-1]) + "\n")
            with self.assertRaises(ValueError):
                list(audit.xdat_frames(path))

    def test_dump_stream_preserves_types_steps_and_rejects_truncation(self):
        header = ("ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n3\n"
                  "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
                  "ITEM: ATOMS id type x y z\n")
        coords = "3 3 3 0 0\n1 1 1 0 0\n2 2 2 0 0\n"
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "misleading_512atoms.dump"
            text = header.format(step=0) + coords + header.format(step=5000) + coords
            path.write_text(text)
            frames = list(audit.stream_lammps_frames(path, ["Ge", "Sb", "Te"]))
            self.assertEqual([a.info["timestep"] for a in frames], [0, 5000])
            self.assertEqual(frames[0].get_chemical_symbols(), ["Ge", "Sb", "Te"])
            self.assertEqual(len(frames[0]), 3)
            path.write_text(text.rsplit("\n", 2)[0] + "\n")
            with self.assertRaises(ValueError):
                list(audit.stream_lammps_frames(path, ["Ge", "Sb", "Te"]))

    def test_archive_duplicate_links_and_resume(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            ref, out = base / "AIMD_reference", base / "inventory"
            package = base / "package"
            ref.mkdir(); package.mkdir()
            original = ref / "XDATCAR-old"
            duplicate = package / "XDATCAR"
            original.write_bytes(b"original trajectory bytes\n")
            duplicate.write_bytes(original.read_bytes())
            target = ref / "GeTe" / "XDATCAR"
            sha = audit.digest(original)
            plan = dict(unique_trajectories=1, material_count=1, records=[dict(
                sha256=sha, archive_file=str(target.relative_to(base)),
                original_paths=[str(original.relative_to(base)), str(duplicate.relative_to(base))])])
            audit.save_json(ref / "archive_manifest.json", plan)
            with patch.multiple(audit, BASE=base, REF=ref, OUT=out):
                audit.organize_references(True)
                self.assertFalse(original.exists())
                self.assertTrue(duplicate.is_symlink())
                self.assertEqual(audit.digest(target), sha)
                audit.organize_references(True)
                self.assertEqual(audit.digest(duplicate), sha)
                target.write_bytes(b"changed bytes")
                with self.assertRaises(ValueError):
                    audit.organize_references(True)


if __name__ == "__main__":
    unittest.main()
