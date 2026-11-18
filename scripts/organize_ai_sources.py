#!/usr/bin/env python3
"""Consolidate local AIMD references and audit AI structure files (pcm_database).

No HPC mutation or website publication. Original package reference paths become
relative symlinks; archive_manifest.json records the original paths and SHA256.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from ase.io import iread, read
from scipy.signal import find_peaks
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "reorg"
REF = BASE / "AIMD_reference"
OUT = BASE / "inventory"
EDGES = np.arange(71) * 0.05
RADII = (EDGES[:-1] + EDGES[1:]) / 2
SHELLS = 4 * np.pi / 3 * np.diff(EDGES**3)


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def save_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def save_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v
                             for k, v in row.items()})


def composition(atoms):
    counts = Counter(atoms.get_chemical_symbols())
    divisor = math.gcd(*counts.values())
    return ";".join(f"{s}:{n // divisor}" for s, n in sorted(counts.items()))


def formula(key):
    return "".join(s + (n if n != "1" else "") for s, n in
                   (part.split(":") for part in key.split(";")))


def xdat_header(path):
    with path.open("rb") as handle:
        header = [handle.readline() for _ in range(7)]
        n = sum(map(int, header[6].split()))
        first = [handle.readline() for _ in range(n + 1)]
    # ASE handles VASP scale and element interpretation, without loading all frames.
    atoms = read(io.StringIO(b"".join(header + first).decode()), format="vasp-xdatcar", index=0)
    return atoms, header


def xdat_frames(path):
    atoms, _ = xdat_header(path)
    n = len(atoms)
    with path.open("rb") as handle:
        for _ in range(7):
            handle.readline()
        index = 0
        while marker := handle.readline():
            if not marker.strip():
                continue
            repeated_header = None
            if not re.fullmatch(rb"Direct\s+configuration\s*=\s*\d+\s*", marker):
                # ASE 3.29 writes a VASP header before every frame. Accept it
                # only after checking that this fixed-cell RDF protocol applies.
                repeated_header = [marker] + [handle.readline() for _ in range(6)]
                marker = handle.readline()
                if not re.fullmatch(rb"Direct\s+configuration\s*=\s*\d+\s*", marker):
                    raise ValueError(f"Unsupported XDATCAR marker: {path}, frame {index}: {marker!r}")
            lines = [handle.readline() for _ in range(n)]
            if any(not line for line in lines):
                raise ValueError(f"Truncated XDATCAR: {path}, frame {index}")
            if repeated_header:
                current = read(io.StringIO(b"".join(repeated_header + [marker] + lines).decode()),
                               format="vasp-xdatcar", index=0)
                if not np.array_equal(current.numbers, atoms.numbers) or not np.allclose(
                        current.cell.array, atoms.cell.array, atol=1e-10, rtol=0):
                    raise ValueError(f"Changing species or cell needs a separate RDF protocol: {path}")
            coords = np.fromstring(b"".join(lines).decode("ascii"), sep=" ")
            if coords.size != 3 * n or not np.isfinite(coords).all():
                raise ValueError(f"Invalid XDATCAR coordinates: {path}, frame {index}")
            yield index, coords.reshape(n, 3) @ atoms.cell.array
            index += 1


def rdf(positions, cell):
    if not np.allclose(cell, np.diag(np.diag(cell)), atol=1e-9, rtol=0):
        raise ValueError("This reference protocol requires an orthorhombic cell")
    box = np.diag(cell)
    if np.min(box) <= 2 * EDGES[-1]:
        raise ValueError("Cell too small for RDF cutoff")
    positions = positions % box
    pairs = cKDTree(positions, boxsize=box).query_pairs(EDGES[-1], output_type="ndarray")
    delta = positions[pairs[:, 1]] - positions[pairs[:, 0]]
    delta -= box * np.rint(delta / box)
    hist = np.histogram(np.linalg.norm(delta, axis=1), bins=EDGES)[0]
    n = len(positions)
    return hist / (n * (n - 1) / 2 * SHELLS / np.linalg.det(cell))


def peak(curve):
    smoothed = np.convolve(np.pad(curve, 2, mode="edge"), np.ones(5) / 5, mode="valid")
    window = (RADII >= 2) & (RADII <= 3.3)
    indices = np.flatnonzero(window)
    k = indices[np.argmax(smoothed[window])]
    if smoothed[k] <= 0 or k in (indices[0], indices[-1]):
        return None, "boundary_or_absent_peak"
    maxima, _ = find_peaks(smoothed, height=0.5 * smoothed[k], prominence=0.1 * smoothed[k], distance=3)
    maxima = maxima[window[maxima]]
    if len(maxima) != 1 or maxima[0] != k:
        return None, "ambiguous_first_shell_peak"
    den = smoothed[k - 1] - 2 * smoothed[k] + smoothed[k + 1]
    shift = 0.5 * (smoothed[k - 1] - smoothed[k + 1]) / den if den < 0 else 0
    return float(RADII[k] + np.clip(shift, -0.5, 0.5) * 0.05), "clear"


def reference_plan():
    indexed = {}
    with (BASE / "diffusion_rdf_dataset/references.csv").open() as handle:
        for row in csv.DictReader(handle):
            if row["included"] == "True":
                indexed[row["sha"]] = row
    candidates = sorted(p for directory in [REF, BASE / "2026_AI4S_uniform_400",
                                             BASE / "diffusion_rdf_dataset/dft_reference"]
                        for p in directory.rglob("*") if p.is_file() and not p.is_symlink()
                        and "XDATCAR" in p.name.upper())
    grouped = defaultdict(list)
    for path in candidates:
        grouped[digest(path)].append(path)
    records = []
    for sha, paths in grouped.items():
        original = indexed.get(sha)
        preferred = next((p for p in paths if p.name.startswith("XDATCAR_DFT")), paths[0])
        atoms, _ = xdat_header(preferred)
        key = composition(atoms)
        destination = REF / formula(key) / preferred.name
        if any(r["archive_file"] == str(destination.relative_to(BASE)) for r in records):
            destination = destination.with_name(destination.name + "__" + sha[:12])
        frames = sum(1 for _ in xdat_frames(preferred))
        rec = dict(reference_id="A_" + sha[:16], material=formula(key), composition=key,
                   natoms=len(atoms), frames=frames, volume_A3=atoms.get_volume(),
                   number_density_Atom_per_Angstrom3=len(atoms) / atoms.get_volume(),
                   temperature_K=int(original["temperature_K"]) if original and original["temperature_K"] else None,
                   phase=None, contributor=None, sha256=sha, bytes=preferred.stat().st_size,
                   archive_file=str(destination.relative_to(BASE)),
                   original_paths=[str(p.relative_to(BASE)) for p in paths],
                   source_reference_metadata=original)
        records.append(rec)
        print(f"reference {len(records)}/{len(grouped)} {rec['material']} {frames} frames", flush=True)
    return dict(schema=1, original_file_count=len(candidates), unique_trajectories=len(records),
                material_count=len({r["material"] for r in records}), records=records,
                source_policy="Original paths outside AIMD_reference retained as relative symlinks; no HPC changes")


def organize_references(apply):
    manifest = REF / "archive_manifest.json"
    if manifest.exists():
        plan = json.loads(manifest.read_text())
        known = {r["archive_file"] for r in plan["records"]}
        known.update(p for r in plan["records"] for p in r["original_paths"])
        unexpected = [str(p.relative_to(BASE)) for directory in
                      [REF, BASE / "2026_AI4S_uniform_400", BASE / "diffusion_rdf_dataset/dft_reference"]
                      for p in directory.rglob("*") if p.is_file() and not p.is_symlink()
                      and "XDATCAR" in p.name.upper() and str(p.relative_to(BASE)) not in known]
        if unexpected:
            raise ValueError(f"New references require an updated archive plan: {unexpected}")
    else:
        plan = reference_plan()
        save_json(OUT / "aimd_reorganization_plan.json", plan)
    if not apply:
        print("Dry run: use references --apply to consolidate files")
        return
    # Full plan is durable before any rename; rerunning resumes this same plan.
    save_json(manifest, plan)
    for record in plan["records"]:
        target = BASE / record["archive_file"]
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and digest(target) != record["sha256"]:
            raise ValueError(f"Destination differs: {target}")
        for original in record["original_paths"]:
            source = BASE / original
            if source == target:
                continue
            if source.is_symlink():
                if source.resolve() != target.resolve():
                    raise ValueError(f"Unexpected symlink: {source}")
                continue
            if source.exists():
                if digest(source) != record["sha256"]:
                    raise ValueError(f"Source differs: {source}")
                if not target.exists():
                    source.rename(target)
                else:
                    source.unlink()  # Verified identical local duplicate, archived bytes preserved.
            elif not target.exists():
                raise FileNotFoundError(source)
            if not source.is_relative_to(REF):
                source.symlink_to(os.path.relpath(target, source.parent))
        if digest(target) != record["sha256"]:
            raise ValueError(f"Archive verification failed: {target}")
    save_csv(REF / "references.csv", plan["records"])
    save_json(REF / "references.json", plan["records"])
    print(f"Archived {plan['unique_trajectories']} trajectories in {plan['material_count']} materials")


def build_reference_rdfs():
    plan = json.loads((REF / "archive_manifest.json").read_text())
    approved_300K = all(r.get("temperature_basis") == "user_assumed_300K" for r in plan["records"])
    materials = defaultdict(list)
    details = []
    for rec in plan["records"]:
        path = BASE / rec["archive_file"]
        atoms, _ = xdat_header(path)
        selected = set(np.linspace(max(0, rec["frames"] - 1000), rec["frames"] - 1,
                                   min(200, rec["frames"]), dtype=int).tolist())
        curves = [rdf(pos, atoms.cell.array) for i, pos in xdat_frames(path) if i in selected]
        if len(curves) != len(selected):
            raise ValueError(f"Reference frame count changed: {path}")
        curve = np.mean(curves, axis=0)
        value, status = peak(curve)
        first, _ = peak(np.mean(curves[:len(curves) // 2], axis=0))
        last, _ = peak(np.mean(curves[len(curves) // 2:], axis=0))
        drift = 100 * abs(first - last) / value if None not in (first, last, value) else None
        detail = dict(reference_id=rec["reference_id"], material=rec["material"],
                      archive_file=rec["archive_file"], sampled_indices=sorted(selected),
                      sampled_frames=len(curves), peak_A=value, peak_status=status,
                      first_half_peak_A=first, last_half_peak_A=last, drift_pct=drift,
                      mean_gr=curve.tolist(), temperature_K=rec["temperature_K"])
        details.append(detail)
        materials[rec["material"]].append(detail)
        print(f"RDF {len(details)}/{len(plan['records'])} {rec['material']}: {value} ({status})", flush=True)
    summaries = []
    for material, records in sorted(materials.items()):
        curve = np.mean([r["mean_gr"] for r in records], axis=0)
        value, status = peak(curve)
        flags = []
        if any(r["temperature_K"] is None for r in records):
            flags.append("temperature_unconfirmed")
        if any(r["peak_status"] != "clear" or r["drift_pct"] is None or r["drift_pct"] >= 5 for r in records):
            flags.append("individual_reference_peak_needs_review")
        summary = dict(material=material, unique_trajectories=len(records),
                       sampled_frames=sum(r["sampled_frames"] for r in records),
                       reference_peak_A=value, peak_status=status,
                       status="material_pooled_reference_300K" if approved_300K else "provisional_pooled_reference",
                       temperature_K=300 if approved_300K else None,
                       temperature_basis="user_assumed_300K" if approved_300K else None, review_flags=flags,
                       reference_ids=[r["reference_id"] for r in records], mean_gr=curve.tolist())
        summaries.append(summary)
        save_json(REF / material / "mean_rdf.json", dict(r_A=RADII.tolist(), **summary))
        save_csv(REF / material / "mean_rdf.csv", [dict(r_A=float(r), gr=float(g)) for r, g in zip(RADII, curve)])
    protocol = dict(bin_A=0.05, rdf_max_A=3.5, peak_search_A=[2.0, 3.3], smoothing_bins=5,
                    interpolation="three_point_parabola", samples_per_trajectory=200,
                    sampling="uniform_0based_indices_from_last_1000_stored_frames",
                    averaging="equal_weight_trajectory_means_after_file_SHA256_deduplication",
                    normalization="finite_N_total_number_RDF_N*(N-1)/2_pairs; periodic_MIC",
                    purpose=("user_approved_material_pooled_300K_reference" if approved_300K else
                             "provisional_material_pooled_baseline; condition matching still requires review"),
                    existing_diffusion_errors="retained in original structures.csv; not recomputed against pooled baseline")
    save_json(REF / "rdf_summary.json", dict(protocol=protocol, materials=summaries, trajectories=details))
    save_csv(REF / "rdf_summary.csv", [{k: v for k, v in s.items() if k != "mean_gr"} for s in summaries])


def geometry_id(atoms):
    xyz = np.rint(atoms.positions * 1e6).astype("<i8")
    labels = np.array(atoms.get_chemical_symbols(), dtype="U3")
    order = np.lexsort((xyz[:, 2], xyz[:, 1], xyz[:, 0], labels))
    h = hashlib.sha256(np.rint(atoms.cell.array.reshape(-1) * 1e6).astype("<i8").tobytes())
    h.update(labels[order].astype("S3").tobytes())
    h.update(xyz[order].tobytes())
    return h.hexdigest()


def stream_lammps_frames(path, specorder):
    """ASE parses one complete dump block at a time, bounding memory by N atoms."""
    with Path(path).open() as handle:
        index = 0
        while first := handle.readline():
            if not first.strip():
                continue
            header = [first] + [handle.readline() for _ in range(8)]
            if (header[0].strip() != "ITEM: TIMESTEP" or
                    header[2].strip() != "ITEM: NUMBER OF ATOMS" or
                    not header[4].startswith("ITEM: BOX BOUNDS") or
                    not header[8].startswith("ITEM: ATOMS")):
                raise ValueError(f"Invalid dump header: {path}, frame {index}")
            n = int(header[3])
            if n < 2:
                raise ValueError(f"Invalid atom count: {path}, frame {index}")
            lines = [handle.readline() for _ in range(n)]
            if any(not line or line.startswith("ITEM:") for line in lines):
                raise ValueError(f"Truncated dump: {path}, frame {index}")
            atoms = read(io.StringIO("".join(header + lines)), format="lammps-dump-text",
                         specorder=specorder)
            if len(atoms) != n:
                raise ValueError(f"Dump atom-count mismatch: {path}, frame {index}")
            yield atoms
            index += 1


def timeline_metadata(steps):
    differences = [b - a for a, b in zip(steps, steps[1:])]
    return dict(first_timestep=steps[0] if steps else None,
                last_timestep=steps[-1] if steps else None,
                min_timestep=min(steps) if steps else None,
                max_timestep=max(steps) if steps else None,
                timestep_differences=dict(Counter(differences)),
                timestep_order=("unknown" if len(steps) < 2 else "strictly_increasing"
                                if all(d > 0 for d in differences) else "non_monotonic_or_reset"))


def inventory_local():
    plan = json.loads((REF / "archive_manifest.json").read_text())
    reference_keys = {r["composition"] for r in plan["records"]}
    original_index = {}
    with (BASE / "diffusion_rdf_dataset/structures.csv").open() as handle:
        for row in csv.DictReader(handle):
            key = ("diffusion_rdf_dataset/" + row["output_file"], int(row["output_frame"]))
            if key in original_index:
                raise ValueError(f"Duplicate index locator: {key}")
            original_index[key] = row
    files = sorted(p for root in [BASE / "2026_AI4S_uniform_400", BASE / "diffusion_rdf_dataset/generated",
                                  BASE / "GST225"]
                   for p in root.rglob("*") if p.suffix.lower() in {".xyz", ".extxyz", ".dump"})
    remote_path = OUT / "hopper_counts.json"
    remote_files = json.loads(remote_path.read_text())["files"] if remote_path.exists() else []
    remote_by_sha = defaultdict(list)
    for remote in remote_files:
        remote_by_sha[remote.get("sha256")].append(remote)
    rows, file_rows, seen = [], [], {}
    for path in files:
        rel = str(path.relative_to(BASE))
        gst = rel.startswith("GST225/")
        specorder = ["Ge", "Sb", "Te"] if gst else ["Ge", "Te"]
        before = path.stat()
        file_sha = digest(path)
        matches = remote_by_sha.get(file_sha, [])
        known_aging = any("/1-GST225-aging/" in r["path"] for r in matches)
        declared = re.search(r"_(\d+)atoms", path.name)
        declared_natoms = int(declared[1]) if declared else None
        trajectory = (stream_lammps_frames(path, specorder) if path.suffix == ".dump"
                      else iread(path, index=":", format="extxyz"))
        start = len(rows)
        for i, atoms in enumerate(trajectory):
            if len(atoms) < 2 or not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.cell.array).all() or atoms.get_volume() <= 0:
                raise ValueError(f"Invalid structure: {rel}, frame {i}")
            key = composition(atoms)
            geom = geometry_id(atoms)
            source = original_index.get((rel, i))
            if rel.startswith("diffusion_rdf_dataset/"):
                if source is None or source["geometry"] != geom or int(source["n"]) != len(atoms) or source["composition"] != key:
                    raise ValueError(f"Original index disagrees: {rel}, frame {i}")
                volume = float(source["volume_A3"])
                if not np.isclose(volume, atoms.get_volume(), rtol=1e-8):
                    raise ValueError(f"Original volume disagrees: {rel}, frame {i}")
            uid = f"{rel}#{i}"
            timestep = atoms.info.get("timestep")
            if timestep is None:
                marker = re.search(r"timestep\s*=\s*(\d+)", atoms.info.get("comment", ""))
                timestep = int(marker[1]) if marker else None
            row = dict(structure_id="AI_" + geom[:20], file=rel, frame_index_0based=i,
                       method=source["generation_method"] if source else "MLIP-MD",
                       material=formula(key), composition=key, natoms=len(atoms),
                       volume_A3=atoms.get_volume(), number_density_Atom_per_Angstrom3=len(atoms) / atoms.get_volume(),
                       temperature_K=(373 if known_aging else int(re.search(r"_(\d+)K_", path.name)[1])
                                      if re.search(r"_(\d+)K_", path.name) else None),
                       phase="amorphous" if known_aging else None, contributor=None,
                       timestep=int(timestep) if timestep is not None else None,
                       timestep_ps=0.002 if known_aging else None,
                       time_ps=float(timestep) * 0.002 if known_aging and timestep is not None else None,
                       pbc=atoms.pbc.tolist(), geometry_sha256=geom,
                       duplicate_of=seen.get(geom), exact_composition_reference=key in reference_keys,
                       legacy_peak_error_pct=float(source["max_error_pct"]) if source else None,
                       legacy_reference_groups=source["reference_group_ids"] if source else None,
                       source_category=(source["source_category"] if source else
                                        "GST225_MLIP_MD" if gst else "uniform400_MLIP_MD"),
                       peak_error_against_pooled_reference_pct=None,
                       status="duplicate" if geom in seen else "geometry_valid_pending_unified_rdf")
            rows.append(row)
            seen.setdefault(geom, uid)
        group = rows[start:]
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError(f"Source changed during inventory: {path}")
        steps = [r["timestep"] for r in group if r["timestep"] is not None]
        file_rows.append(dict(file=rel, bytes=path.stat().st_size, sha256=file_sha,
                              frames=len(group), unique_geometries=sum(not r["duplicate_of"] for r in group),
                              materials=sorted({r["material"] for r in group}),
                              atom_counts=sorted({r["natoms"] for r in group}),
                              declared_filename_natoms=declared_natoms,
                              filename_atom_count_matches=(all(r["natoms"] == declared_natoms for r in group)
                                                           if declared_natoms is not None else None),
                              **timeline_metadata(steps),
                              source_category="GST225_MLIP_MD" if gst else "existing_local_dataset",
                              identical_hopper_paths=[r["path"] for r in matches],
                              methods=dict(Counter(r["method"] for r in group)),
                              exact_composition_reference=all(r["exact_composition_reference"] for r in group)))
        print(f"inventory {len(file_rows)}/{len(files)} {rel}: {len(group)} frames", flush=True)
        save_json(OUT / "ai_local_progress.json", dict(completed_files=len(file_rows), total_files=len(files),
                                                       latest_file=rel, parsed_frames=len(rows)))
    actual_locators = {(r["file"], r["frame_index_0based"]) for r in rows if r["file"].startswith("diffusion_rdf_dataset/")}
    if actual_locators != set(original_index):
        raise ValueError("Original index contains missing/unexpected frames")
    summary = dict(files=len(file_rows), complete_valid_frames=len(rows),
                   unique_geometries=len(seen), duplicate_frames=len(rows) - len(seen),
                   by_method=dict(Counter(r["method"] for r in rows)),
                   unique_by_method=dict(Counter(r["method"] for r in rows if not r["duplicate_of"])),
                   by_material=dict(Counter(r["material"] for r in rows)),
                   by_natoms=dict(Counter(r["natoms"] for r in rows)),
                   by_source_category=dict(Counter(r["source_category"] for r in rows)),
                   missing_exact_reference_frames=sum(not r["exact_composition_reference"] for r in rows),
                   unique_with_exact_reference=sum(r["exact_composition_reference"] and not r["duplicate_of"] for r in rows),
                   deduplication="species+raw_cartesian_coordinates+cell at 1e-6 A; atom-order invariant; no structural similarity clustering",
                   counting="stored complete frames; correlated MD frames are not statistically independent structures",
                   rdf_status="New pooled-reference error not yet calculated; legacy errors use original reference groups")
    save_json(OUT / "ai_local_files.json", file_rows)
    save_csv(OUT / "ai_local_files.csv", file_rows)
    save_json(OUT / "ai_local_frames.json", rows)
    save_csv(OUT / "ai_local_frames.csv", rows)
    save_json(OUT / "ai_local_summary.json", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["references", "rdf", "local"])
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    if args.command == "references":
        organize_references(args.apply)
    elif args.command == "rdf":
        build_reference_rdfs()
    else:
        inventory_local()


if __name__ == "__main__":
    main()
