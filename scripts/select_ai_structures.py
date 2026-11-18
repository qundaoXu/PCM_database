#!/usr/bin/env python3
"""Diversity selection, 300 K reference errors, and lossless structure export.

Run in pcm_database: python scripts/select_ai_structures.py build
Source structures remain untouched. Cached descriptors are tied to source SHA256
and parameters. Extracted frames retain their original bytes and provenance.
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
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from ase.io import read
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits

try:
    from .organize_ai_sources import BASE, REF, OUT, EDGES, RADII, SHELLS, digest, geometry_id, peak, save_csv, save_json
except ImportError:
    from organize_ai_sources import BASE, REF, OUT, EDGES, RADII, SHELLS, digest, geometry_id, peak, save_csv, save_json

DEST = BASE / "AI_selected"
CACHE = OUT / "diversity_cache"
SPECIES = ["Ge", "Sb", "Se", "Te"]
PAIR_TYPES = [(i, j) for i in range(4) for j in range(i, 4)]
PAIR_INDEX = np.empty((4, 4), dtype=int)
for index, (i, j) in enumerate(PAIR_TYPES):
    PAIR_INDEX[i, j] = PAIR_INDEX[j, i] = index
FEATURE_EDGES = np.linspace(1.5, 5.5, 41)
FEATURE_SHELLS = 4 * np.pi / 3 * np.diff(FEATURE_EDGES**3)
PARAMS = dict(version=1, cutoff_A=5.5, rdf_feature_bin_A=0.1, rdf_feature_min_A=1.5,
              species=SPECIES, pair_types=PAIR_TYPES, coordination_cutoff_A=3.5,
              coordination_overflow_bin=13, density_multiplier=100, coordination_multiplier=5,
              peak_bin_A=0.05, peak_search_A=[2, 3.3], smoothing_bins=5,
              distance_round_decimals=12,
              descriptor_dimension=511, pca_dimensions="adaptive", pca_variance_target=0.95, pca_scale_floor=0.05,
              random_seed=20261117)


def raw_frames(path):
    """Yield complete original blocks as bytes, without whole-trajectory loading."""
    dump = Path(path).suffix == ".dump"
    with Path(path).open("rb") as handle:
        index = 0
        while first := handle.readline():
            if not first.strip():
                continue
            if dump:
                header = [first] + [handle.readline() for _ in range(8)]
                if (header[0].strip() != b"ITEM: TIMESTEP" or header[2].strip() != b"ITEM: NUMBER OF ATOMS"
                        or not header[4].startswith(b"ITEM: BOX BOUNDS") or not header[8].startswith(b"ITEM: ATOMS")):
                    raise ValueError(f"Invalid dump header: {path} frame {index}")
                n = int(header[3])
            else:
                n = int(first)
                header = [first, handle.readline()]
            lines = [handle.readline() for _ in range(n)]
            if n < 2 or any(not line or line.startswith(b"ITEM:") for line in lines) or not header[-1]:
                raise ValueError(f"Incomplete structure: {path} frame {index}")
            yield index, b"".join(header + lines)
            index += 1


def parse_frame(raw, dump, gst):
    kwargs = {"specorder": ["Ge", "Sb", "Te"] if gst else ["Ge", "Te"]} if dump else {}
    return read(io.StringIO(raw.decode()), format="lammps-dump-text" if dump else "extxyz", **kwargs)


def modeling_method(row):
    path = row["file"]
    if path.startswith("GST225/"):
        return "MLIP_MACE"
    if path.startswith("2026_AI4S_uniform_400/gete/"):
        return "MLIP_DeePMD"
    if path.startswith("2026_AI4S_uniform_400/sb-te/"):
        return "MLIP_NEP"
    if row["method"] == "diffusion":
        return "Diffusion"
    if row["method"] == "flow_matching":
        return "Flow Matching"
    raise ValueError(f"Unknown modeling method: {row}")


def time_order_label(time_ps):
    if time_ps is None:
        return "none"
    if time_ps < 0 or not np.isfinite(time_ps):
        raise ValueError("Invalid simulation time")
    if time_ps == 0:
        return "0ps"
    exponent = int(np.floor(np.log10(time_ps) + 1e-12))
    return f"{10.0 ** (exponent - 3):g}ns" if exponent >= 3 else f"{10.0 ** exponent:g}ps"


def normalize_candidate(row):
    row = dict(row)
    row["modeling_method"] = modeling_method(row)
    row.setdefault("reported_temperature_K", row.get("temperature_K"))
    row["temperature_K"] = 300
    row["temperature_basis"] = "user_normalized_300K"
    row["simulation_time_ps"] = row.get("time_ps")
    row["simulation_time_order"] = time_order_label(row.get("time_ps"))
    row.setdefault("simulation_time_basis", "verified_HPC_timestep_ps" if row.get("time_ps") is not None else None)
    return row


def header_time_fs(raw):
    header = raw.split(b"\n", 2)[1]
    match = re.search(rb'(?:^|\s)Time=([^\s]+)', header)
    return float(match[1]) if match else None


def enrich_gpumd_times(rows):
    """Known NEP XYZ header schema matches GPUMD's Time-in-fs writer.

    Unit evidence: official src/measure/dump_observer.cu::output_line2;
    format-based provenance is explicit, not a guessed dt from frame spacing.
    """
    lookup = defaultdict(dict)
    for row in rows:
        if row["file"].startswith("2026_AI4S_uniform_400/sb-te/"):
            lookup[row["file"]][row["frame_index_0based"]] = row
    for rel, wanted in lookup.items():
        for frame, raw in raw_frames(BASE / rel):
            value = header_time_fs(raw)
            if value is None:
                continue
            row = wanted[frame]
            row.update(reported_simulation_time=value, reported_simulation_time_unit="fs",
                       simulation_time_ps=value / 1000,
                       simulation_time_order=time_order_label(value / 1000),
                       simulation_time_basis="GPUMD_Time_fs_inferred_from_matching_header_schema")


def dominant_peak(curve):
    value, status = peak(curve)
    if value is not None:
        return value, status
    if status != "ambiguous_first_shell_peak":
        return None, status
    # The requested main peak is the dominant first-shell maximum. Preserve the
    # ambiguity flag when another local peak exists; do not invent a missing peak.
    smooth = np.convolve(np.pad(curve, 2, mode="edge"), np.ones(5) / 5, mode="valid")
    indices = np.flatnonzero((RADII >= 2) & (RADII <= 3.3))
    k = indices[np.argmax(smooth[indices])]
    den = smooth[k - 1] - 2 * smooth[k] + smooth[k + 1]
    shift = 0.5 * (smooth[k - 1] - smooth[k + 1]) / den if den < 0 else 0
    return float(RADII[k] + np.clip(shift, -0.5, 0.5) * 0.05), "dominant_peak_multiple_local_maxima"


def structural_features(atoms):
    cell = atoms.cell.array
    if not atoms.pbc.all() or not np.isfinite(atoms.positions).all():
        raise ValueError("Nonperiodic or nonfinite structure")
    if not np.allclose(cell, np.diag(np.diag(cell)), atol=1e-8, rtol=0):
        raise ValueError("Descriptor protocol requires orthorhombic cells")
    box = np.diag(cell)
    if np.min(box) <= 2 * FEATURE_EDGES[-1]:
        raise ValueError("Cell too small for 5.5 A minimum-image descriptor")
    positions = atoms.positions % box
    pairs = cKDTree(positions, boxsize=box).query_pairs(FEATURE_EDGES[-1], output_type="ndarray")
    delta = positions[pairs[:, 1]] - positions[pairs[:, 0]]
    delta -= box * np.rint(delta / box)
    # Stabilize exact histogram-edge ties against floating-point MIC roundoff.
    distances = np.round(np.linalg.norm(delta, axis=1), 12)
    n, volume = len(atoms), atoms.get_volume()
    total = np.histogram(distances, bins=FEATURE_EDGES)[0] / (n * (n - 1) / 2 * FEATURE_SHELLS / volume)
    species = np.array([SPECIES.index(s) for s in atoms.get_chemical_symbols()])
    populations = np.bincount(species, minlength=4)
    types = PAIR_INDEX[species[pairs[:, 0]], species[pairs[:, 1]]]
    bins = np.searchsorted(FEATURE_EDGES, distances, side="right") - 1
    bins[distances == FEATURE_EDGES[-1]] = len(FEATURE_SHELLS) - 1
    valid = (bins >= 0) & (bins < len(FEATURE_SHELLS))
    hist = np.bincount(types[valid] * len(FEATURE_SHELLS) + bins[valid], minlength=400).reshape(10, 40)
    partial = np.zeros((10, 40))
    for k, (i, j) in enumerate(PAIR_TYPES):
        expected_pairs = populations[i] * (populations[i] - 1) / 2 if i == j else populations[i] * populations[j]
        if expected_pairs:
            partial[k] = hist[k] / (expected_pairs * FEATURE_SHELLS / volume)
    first = pairs[distances <= 3.5]
    degree = np.bincount(first.reshape(-1), minlength=n)
    degree = np.minimum(degree, 13)
    coordination = np.bincount(degree, minlength=14) / n
    species_coordination = np.zeros((4, 14))
    for i in range(4):
        if populations[i]:
            species_coordination[i] = np.bincount(degree[species == i], minlength=14) / populations[i]
    descriptor = np.concatenate([total, partial.ravel(), coordination * 5,
                                 species_coordination.ravel() * 5, [100 * n / volume]]).astype(np.float32)
    curve = np.histogram(distances, bins=EDGES)[0] / (n * (n - 1) / 2 * SHELLS / volume)
    value, status = dominant_peak(curve)
    assert descriptor.shape == (PARAMS["descriptor_dimension"],)
    return descriptor, curve.astype(np.float32), value, status


def cache_task(task):
    rel, source_sha, expected = task
    path = BASE / rel
    before = path.stat()
    if digest(path) != source_sha:
        raise ValueError(f"Source hash changed: {path}")
    tag = hashlib.sha256(json.dumps(dict(source_sha=source_sha, params=PARAMS), sort_keys=True).encode()).hexdigest()
    cache = CACHE / (tag + ".npz")
    # Version-1 feature arrays are identical: only the downstream PCA rule changed.
    previous_params = {k: v for k, v in PARAMS.items() if k != "pca_variance_target"}
    previous_params["pca_dimensions"] = 32
    previous_tag = hashlib.sha256(json.dumps(dict(source_sha=source_sha, params=previous_params), sort_keys=True).encode()).hexdigest()
    if not cache.exists() and (CACHE / (previous_tag + ".npz")).exists():
        cache = CACHE / (previous_tag + ".npz")
    if cache.exists():
        with np.load(cache, allow_pickle=False) as data:
            if data["indices"].tolist() == [r[0] for r in expected]:
                after = path.stat()
                if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                    raise ValueError(f"Source changed during cache validation: {path}")
                return rel, str(cache), len(expected), True
    by_frame = {frame: (global_index, geometry) for global_index, frame, geometry in expected}
    indices, descriptors, curves, peaks, statuses, raw_hashes = [], [], [], [], [], []
    for frame, raw in raw_frames(path):
        if frame not in by_frame:
            continue
        global_index, geom = by_frame[frame]
        atoms = parse_frame(raw, path.suffix == ".dump", rel.startswith("GST225/"))
        if geometry_id(atoms) != geom:
            raise ValueError(f"Inventory geometry differs: {path} frame {frame}")
        descriptor, curve, value, status = structural_features(atoms)
        indices.append(global_index); descriptors.append(descriptor); curves.append(curve)
        peaks.append(np.nan if value is None else value); statuses.append(status)
        raw_hashes.append(hashlib.sha256(raw).hexdigest())
    if indices != [r[0] for r in expected]:
        raise ValueError(f"Missing inventory frame: {path}")
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"Source changed while reading: {path}")
    cache.parent.mkdir(parents=True, exist_ok=True)
    temporary = cache.with_suffix(".npz.tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, indices=np.array(indices), descriptors=np.array(descriptors),
                            rdf=np.array(curves), peaks=np.array(peaks), statuses=np.array(statuses),
                            raw_sha256=np.array(raw_hashes))
    temporary.replace(cache)
    return rel, str(cache), len(indices), False


def allocate_budgets(counts, budget, square_root=True):
    """Capped allocation, never exceeding available items or the global budget."""
    counts = np.asarray(counts, dtype=int)
    if budget < np.count_nonzero(counts) or np.any(counts < 0):
        raise ValueError("Budget must cover one sample per nonempty stratum")
    budget = min(int(budget), int(counts.sum()))
    allocation = (counts > 0).astype(int)
    weights = np.sqrt(counts) if square_root else (counts > 0).astype(float)
    while int(allocation.sum()) < budget:
        active = allocation < counts
        remaining = budget - int(allocation.sum())
        shares = np.where(active, weights, 0)
        shares *= remaining / shares.sum()
        addition = np.minimum(counts - allocation, np.floor(shares).astype(int))
        if addition.sum() == 0:
            order = np.argsort(-shares, kind="stable")
            for k in order[:remaining]:
                if active[k]:
                    allocation[k] += 1
            continue
        allocation += addition
    return allocation


def diversity_sample(features, quota):
    n = len(features)
    if quota >= n:
        return list(range(n)), dict(retained_variance=1.0, pca_dimensions=0, max_distance=0.0)
    scale = np.maximum(np.std(features, axis=0), PARAMS["pca_scale_floor"])
    standardized = (features - np.mean(features, axis=0)) / scale
    nonconstant = np.std(standardized, axis=0) > 1e-7
    x = standardized[:, nonconstant]
    if not x.shape[1]:
        return [0], dict(retained_variance=1.0, pca_dimensions=0, max_distance=0.0)
    with threadpool_limits(limits=1):
        pca = PCA(n_components=PARAMS["pca_variance_target"], svd_solver="full")
        x = pca.fit_transform(x).astype(np.float32)
    dimensions = x.shape[1]
    center = np.mean(x, axis=0)
    selected = [int(np.argmin(np.sum((x - center)**2, axis=1)))]
    distances = np.full(n, np.inf, dtype=np.float32)
    for _ in range(quota):
        delta = x - x[selected[-1]]
        distances = np.minimum(distances, np.einsum("ij,ij->i", delta, delta))
        distances[selected] = -1
        if len(selected) >= quota or float(np.max(distances)) <= 1e-12:
            break
        selected.append(int(np.argmax(distances)))
    coverage = float(np.sqrt(max(0.0, float(np.max(distances)))))
    return selected, dict(retained_variance=float(np.sum(pca.explained_variance_ratio_)),
                          pca_dimensions=dimensions, max_distance=coverage)


def uniform_time_sample(indices, times, quota):
    if quota >= len(indices):
        return list(indices)
    order = np.argsort(times, kind="stable")
    ordered = np.asarray(indices)[order]
    time = np.asarray(times)[order]
    selected, previous = [], -1
    for j, target in enumerate(np.linspace(time[0], time[-1], quota)):
        lo, hi = previous + 1, len(time) - (quota - j)
        upper = int(np.clip(np.searchsorted(time, target), lo, hi))
        lower = max(lo, upper - 1)
        chosen = lower if abs(time[lower] - target) <= abs(time[upper] - target) else upper
        selected.append(int(ordered[chosen])); previous = chosen
    return selected


def time_sample(rows, indices, quota):
    bins = defaultdict(list)
    for i in indices:
        bins[rows[i]["simulation_time_order"]].append(i)
    names = sorted(bins, key=lambda k: min(rows[i]["simulation_time_ps"] for i in bins[k]))
    allocations = allocate_budgets([len(bins[k]) for k in names], quota, square_root=False)
    selected = []
    for name, q in zip(names, allocations):
        group = bins[name]
        selected += uniform_time_sample(group, [rows[i]["simulation_time_ps"] for i in group], int(q))
    return selected


def reference_match(composition, references, tolerance):
    fractions = {}
    for key, target in references.items():
        counts = {s: int(n) for s, n in (p.split(":") for p in key.split(";"))}
        total = sum(counts.values())
        fractions[key] = {s: n / total for s, n in counts.items()}
    if composition in references:
        return references[composition], "exact_composition", 0.0
    counts = {s: int(n) for s, n in (p.split(":") for p in composition.split(";"))}
    total = sum(counts.values()); fraction = {s: n / total for s, n in counts.items()}
    distances = [(max(abs(fraction[s] - f[s]) for s in fraction), key) for key, f in fractions.items() if set(f) == set(fraction)]
    if not distances:
        return None, "no_reference", None
    deviation, key = min(distances)
    if deviation > tolerance:
        return None, "no_reference", float(deviation * 100)
    return references[key], "integer_stoichiometry_rounding", float(deviation * 100)


def normalize_references():
    plan = json.loads((REF / "archive_manifest.json").read_text())
    for r in plan["records"]:
        r.setdefault("reported_temperature_K", r["temperature_K"])
        r.update(temperature_K=300, temperature_basis="user_assumed_300K")
    save_json(REF / "archive_manifest.json", plan)
    save_json(REF / "references.json", plan["records"])
    save_csv(REF / "references.csv", plan["records"])
    summary = json.loads((REF / "rdf_summary.json").read_text())
    summary["protocol"].update(temperature_K=300, temperature_basis="user_assumed_300K",
                               purpose="user_approved_material_pooled_300K_reference")
    for r in summary["trajectories"]:
        r.setdefault("reported_temperature_K", r["temperature_K"])
        r.update(temperature_K=300, temperature_basis="user_assumed_300K")
    result = {}
    by_material = {r["material"]: r["composition"] for r in plan["records"]}
    for r in summary["materials"]:
        r["review_flags"] = [f for f in r["review_flags"] if f != "temperature_unconfirmed"]
        r.update(temperature_K=300, temperature_basis="user_assumed_300K", status="material_pooled_reference_300K")
        r["composition"] = by_material[r["material"]]
        save_json(REF / r["material"] / "mean_rdf.json", dict(r_A=RADII.tolist(), **r))
        result[r["composition"]] = r
    save_json(REF / "rdf_summary.json", summary)
    save_csv(REF / "rdf_summary.csv", [{k: v for k, v in r.items() if k != "mean_gr"} for r in summary["materials"]])
    return result


def csv_metadata(path, rows):
    save_csv(path, [{k: "none" if v is None else v for k, v in r.items()} for r in rows])


def build(workers):
    policy = json.loads((OUT / "collection_policy.json").read_text())
    budget = int(policy["selection_policy"]["maximum_structures"])
    if budget > 15000:
        raise ValueError("Selection exceeds user cap")
    originals = json.loads((OUT / "ai_local_frames.json").read_text())
    rows = [normalize_candidate(r) for r in originals]
    enrich_gpumd_times(rows)
    files = json.loads((OUT / "ai_local_files.json").read_text())
    references = normalize_references()
    by_file = defaultdict(list)
    for i, row in enumerate(rows):
        if not row["duplicate_of"]:
            by_file[row["file"]].append((i, row["frame_index_0based"], row["geometry_sha256"]))
    source_lookup = {r["file"]: r for r in files}
    tasks = [(rel, source_lookup[rel]["sha256"], items) for rel, items in sorted(by_file.items())]
    caches = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(cache_task, task) for task in tasks]
        for future in as_completed(futures):
            rel, cache, n, reused = future.result()
            caches.append((rel, cache))
            progress = dict(stage="features", completed_files=len(caches), total_files=len(tasks), latest_file=rel)
            save_json(OUT / "selection_progress.json", progress)
            print(f"features {len(caches)}/{len(tasks)} {n} frames {rel} cache={reused}", flush=True)
    matrix = np.zeros((len(rows), PARAMS["descriptor_dimension"]), dtype=np.float32)
    curves = np.zeros((len(rows), len(RADII)), dtype=np.float32)
    measured = {}
    for rel, cache in caches:
        with np.load(cache, allow_pickle=False) as data:
            indices = data["indices"]
            matrix[indices] = data["descriptors"]; curves[indices] = data["rdf"]
            for i, value, status, raw_sha in zip(indices, data["peaks"], data["statuses"], data["raw_sha256"]):
                measured[int(i)] = (None if np.isnan(value) else float(value), str(status), str(raw_sha))
    eligible, matches, ineligible = defaultdict(list), {}, []
    for i, row in enumerate(rows):
        if row["duplicate_of"]:
            continue
        reference, match, deviation = reference_match(row["composition"], references,
                           policy["selection_policy"]["nearest_reference_max_atomic_fraction_deviation"])
        value, status, _ = measured[i]
        if reference is None or value is None:
            ineligible.append(dict(file=row["file"], frame=row["frame_index_0based"], reason=match if reference is None else status))
            continue
        matches[i] = reference, match, deviation
        key = (reference["material"], row["modeling_method"], row["natoms"])
        eligible[key].append(i)
    keys = sorted(eligible)
    allocations = allocate_budgets([len(eligible[k]) for k in keys], budget)
    selected, groups, selected_by = [], [], {}
    for j, (key, quota) in enumerate(zip(keys, allocations)):
        indices = eligible[key]
        known_times = [rows[i]["simulation_time_ps"] for i in indices]
        if all(t is not None for t in known_times) and max(known_times) - min(known_times) >= 1000:
            chosen = time_sample(rows, indices, int(quota))
            method, coverage = "physical_time_uniform_within_time_orders", None
        else:
            chosen_local, coverage = diversity_sample(matrix[indices], int(quota))
            chosen = [indices[k] for k in chosen_local]
            method = "diversity_farthest_point_scaled_PCA"
        selected += chosen
        for i in chosen:
            selected_by[i] = method
        groups.append(dict(material=key[0], modeling_method=key[1], natoms=key[2],
                           candidates=len(indices), allocated=int(quota), selected=len(chosen),
                           sampling_method=method, descriptor_coverage=coverage))
        save_json(OUT / "selection_progress.json", dict(stage="selection", completed_groups=j + 1,
                                                       total_groups=len(keys), selected=len(selected)))
        print(f"select {j+1}/{len(keys)} {key}: {len(indices)} -> {len(chosen)}", flush=True)
    if len(selected) > budget or len(set(selected)) != len(selected):
        raise ValueError("Invalid selection")
    selected_set = set(selected)
    for i, row in enumerate(rows):
        row["selected_for_ai_module"] = i in selected_set
        if i in measured:
            row["first_shell_peak_A"], row["first_shell_peak_status"], _ = measured[i]
        if i in matches:
            reference, match, deviation = matches[i]
            row.update(reference_material=reference["material"], reference_match=match,
                       reference_composition_deviation_percentage_points=deviation,
                       peak_error_against_pooled_reference_pct=100 * abs(measured[i][0] - reference["reference_peak_A"]) / reference["reference_peak_A"])
    save_json(OUT / "ai_local_frames.json", rows)
    csv_metadata(OUT / "ai_local_frames.csv", rows)
    DEST.mkdir(parents=True, exist_ok=True)
    metadata, exported_hashes, handles, positions = [], {}, {}, Counter()
    selected_files = defaultdict(dict)
    for i in selected:
        selected_files[rows[i]["file"]][rows[i]["frame_index_0based"]] = i
    try:
        for rel, wanted in sorted(selected_files.items()):
            path = BASE / rel
            for frame, raw in raw_frames(path):
                if frame not in wanted:
                    continue
                i = wanted[frame]; row = rows[i]
                reference, match, deviation = matches[i]
                value, status, raw_sha = measured[i]
                if hashlib.sha256(raw).hexdigest() != raw_sha:
                    raise ValueError(f"Source changed before export: {rel} frame {frame}")
                model = row["modeling_method"].replace(" ", "_")
                ext = ".dump" if path.suffix == ".dump" else ".extxyz"
                target_rel = Path("structures") / model / f"{row['material']}_N{row['natoms']}__{source_lookup[rel]['sha256'][:12]}{ext}"
                if target_rel not in handles:
                    (DEST / target_rel).parent.mkdir(parents=True, exist_ok=True)
                    handles[target_rel] = (DEST / target_rel).open("wb")
                handles[target_rel].write(raw)
                output_frame = positions[str(target_rel)]; positions[str(target_rel)] += 1
                metadata.append(dict(structure_id=row["structure_id"], material=row["material"],
                    composition=row["composition"], modeling_method=row["modeling_method"],
                    generation_method=row["method"], phase=row["phase"], temperature_K=300,
                    reported_temperature_K=row["reported_temperature_K"], temperature_basis="user_normalized_300K",
                    natoms=row["natoms"], number_density_Atom_per_Angstrom3=row["number_density_Atom_per_Angstrom3"],
                    volume_A3=row["volume_A3"], contributor=row["contributor"],
                    structure_file=str(target_rel), structure_frame_0based=output_frame,
                    structure_format="lammps-dump-text" if ext == ".dump" else "extxyz",
                    lammps_type_map={1: "Ge", 2: "Sb", 3: "Te"} if rel.startswith("GST225/") and ext == ".dump"
                                    else {1: "Ge", 2: "Te"} if ext == ".dump" else None,
                    source_file=rel, source_frame_0based=frame, source_file_sha256=source_lookup[rel]["sha256"],
                    raw_frame_sha256=raw_sha, geometry_sha256=row["geometry_sha256"],
                    source_hopper_paths=source_lookup[rel]["identical_hopper_paths"],
                    source_timestep=row["timestep"], timestep_ps=row["timestep_ps"],
                    simulation_time_ps=row["simulation_time_ps"], simulation_time_order=row["simulation_time_order"],
                    simulation_time_basis=row["simulation_time_basis"],
                    reported_simulation_time=row.get("reported_simulation_time"),
                    reported_simulation_time_unit=row.get("reported_simulation_time_unit"),
                    reference_material=reference["material"], reference_temperature_K=300,
                    reference_ids=reference["reference_ids"], reference_match=match,
                    reference_composition_deviation_percentage_points=deviation,
                    rdf_first_shell_peak_A=value, aimd_first_shell_peak_A=reference["reference_peak_A"],
                    peak_position_relative_error_pct=100 * abs(value - reference["reference_peak_A"]) / reference["reference_peak_A"],
                    peak_position_signed_error_pct=100 * (value - reference["reference_peak_A"]) / reference["reference_peak_A"],
                    rdf_peak_status=status, rdf_protocol="total_number_RDF_0.05A_5bin_smoothing_first_shell_2.0_to_3.3A_v1",
                    selection_method=selected_by[i], selection_stratum=[reference["material"], row["modeling_method"], row["natoms"]]))
    finally:
        for handle in handles.values():
            handle.close()
    for path in positions:
        exported_hashes[path] = digest(DEST / path)
    metadata.sort(key=lambda r: (r["modeling_method"], r["material"], r["natoms"], r["source_file"], r["source_frame_0based"]))
    save_json(DEST / "metadata.json", metadata)
    csv_metadata(DEST / "metadata.csv", metadata)
    manifest = dict(protocol=PARAMS, policy=policy, candidates=len(measured), eligible=sum(map(len, eligible.values())),
                    selected=len(metadata), cap=budget, ineligible=ineligible, strata=groups,
                    by_modeling_method=dict(Counter(r["modeling_method"] for r in metadata)),
                    by_material=dict(Counter(r["material"] for r in metadata)),
                    by_simulation_time_order=dict(Counter(r["simulation_time_order"] for r in metadata)),
                    peak_statuses=dict(Counter(r["rdf_peak_status"] for r in metadata)),
                    near_stoichiometry_reference_matches=sum(r["reference_match"] != "exact_composition" for r in metadata),
                    output_file_sha256=exported_hashes,
                    output_file_frames=dict(positions), source_file_sha256={r["file"]: r["sha256"] for r in files},
                    cache_files={rel: str(Path(path).relative_to(BASE)) for rel, path in caches}, temporal_independence_claim=False,
                    error_formula="100*abs(r_AI-r_AIMD)/r_AIMD", error_threshold_applied=False)
    save_json(DEST / "selection_manifest.json", manifest)
    save_csv(DEST / "selection_strata.csv", groups)
    inventory_summary = json.loads((OUT / "ai_local_summary.json").read_text())
    inventory_summary.update(by_modeling_method=dict(Counter(r["modeling_method"] for r in rows if not r["duplicate_of"])),
                             comparison_temperature_K=300, temperature_basis="user_normalized_300K",
                             selected_structures=len(metadata),
                             rdf_status="all_unique_candidates_evaluated_against_user_approved_300K_material_reference")
    save_json(OUT / "ai_local_summary.json", inventory_summary)
    (DEST / "README.md").write_text(
        f"# AI 生成结构：多样性筛选结果\n\n从 {len(measured):,} 个几何去重候选中选择 **{len(metadata):,}** 个结构（上限 {budget:,}）。\n\n"
        "`metadata.json` / `metadata.csv` 保存建模方式、源文件/帧号/哈希、300 K 对应参考、第一主峰及百分比误差。\n"
        "建模方式为 MLIP_MACE、MLIP_DeePMD、MLIP_NEP、Diffusion、Flow Matching。\n"
        "AIMD 及 MLIP 对应温度按用户约定统一为 300 K；原报告温度另保留，不表示重新测得温度。\n"
        "未知模拟时间在 JSON 的数值字段为 null，CSV 为 none；simulation_time_order 两种格式均为 none。\n"
        "已知时间以真实 ps 值记录，并按 10 的幂标注数量级，例如 1ps、10ps、100ps、1ns、10ns；0 时刻单列 0ps。\n\n"
        "Sb–Te XYZ 的 Time 字段按匹配 GPUMD 官方输出格式推断为 fs，再除以 1000 转 ps；\n"
        "simulation_time_basis 显式保存格式推断依据，不根据存储帧间隔反推 dt。\n"
        "单位来源：[GPUMD 官方 Time 输出代码](https://github.com/brucefan1983/GPUMD/blob/master/src/measure/dump_observer.cu)。\n\n"
        "按参考材料、建模方式和真实原子数分层，用候选数平方根分配名额并保留小类。\n"
        "511 维描述符包含元素分辨的 1.5–5.5 Å 径向分布、3.5 Å 配位数分布及数密度；\n"
        "分层缩放后 PCA 自适应保留至少 95% 方差，再做确定性的最远点采样。描述符侧重径向和配位多样性，不保证识别所有角度/拓扑差异。\n"
        "已知时间的老化轨迹先按时间数量级分层，每层按物理时间均匀选点，避免按非均匀存储帧索引选点。\n"
        "strata CSV 和 manifest 保存候选数、配额、实际采样数、PCA 方差覆盖及最远覆盖距离。\n\n"
        "峰位统一使用总数目 g(r)，0.05 Å 分箱、有限 N 原子对归一、5 点滑动平均，\n"
        "2.0–3.3 Å 第一壳主峰三点插值；误差为 100*abs(r_AI-r_AIMD)/r_AIMD。\n"
        "多个局部峰时记录占优主峰及歧义标记；未能定位内部主峰者不进入本次可评估集合。未按误差大小设筛选门槛。\n"
        "Sb667Te333 / Sb333Te667 按最近理想组分对应参考，最大原子分数偏差容差为 0.001（0.1 百分点），\n"
        "真实组分、参考组分和实际偏差全部保留。\n\n"
        "structures/ 内保留原始帧字节；LAMMPS dump 须使用 metadata 的 type map 用 ASE 读取。\n"
        "没有删除或裁剪原始轨迹，也没有发布网站。\n\n"
        "```bash\nconda activate pcm_database\npython scripts/select_ai_structures.py build --workers 2\npython scripts/select_ai_structures.py verify\n```\n")
    print(f"COMPLETE selected={len(metadata)} methods={manifest['by_modeling_method']}", flush=True)


def verify():
    manifest = json.loads((DEST / "selection_manifest.json").read_text())
    rows = json.loads((DEST / "metadata.json").read_text())
    if len(rows) > 15000 or len(rows) != manifest["selected"] or len({r["geometry_sha256"] for r in rows}) != len(rows):
        raise ValueError("Count or geometry uniqueness failed")
    for rel, sha in manifest["source_file_sha256"].items():
        if digest(BASE / rel) != sha:
            raise ValueError(f"Source hash changed since inventory: {rel}")
    for stratum in manifest["strata"]:
        coverage = stratum["descriptor_coverage"]
        if coverage and coverage["retained_variance"] < PARAMS["pca_variance_target"] - 1e-6:
            raise ValueError("Insufficient descriptor variance retained")
    lookup = defaultdict(dict)
    for r in rows:
        lookup[r["structure_file"]][r["structure_frame_0based"]] = r
        expected = 100 * abs(r["rdf_first_shell_peak_A"] - r["aimd_first_shell_peak_A"]) / r["aimd_first_shell_peak_A"]
        if (abs(expected - r["peak_position_relative_error_pct"]) > 1e-12 or r["temperature_K"] != 300
                or r["reference_temperature_K"] != 300):
            raise ValueError("Metadata error or temperature mismatch")
        if time_order_label(r["simulation_time_ps"]) != r["simulation_time_order"]:
            raise ValueError("Time order mismatch")
        if r["simulation_time_ps"] is not None and r["timestep_ps"] is not None and r["simulation_time_ps"] != r["source_timestep"] * r["timestep_ps"]:
            raise ValueError("Time conversion mismatch")
    count = 0
    for rel, expected_sha in manifest["output_file_sha256"].items():
        path = DEST / rel
        if digest(path) != expected_sha:
            raise ValueError(f"Export hash mismatch: {path}")
        n = 0
        for frame, raw in raw_frames(path):
            r = lookup[rel][frame]
            if hashlib.sha256(raw).hexdigest() != r["raw_frame_sha256"]:
                raise ValueError("Original frame bytes differ")
            atoms = parse_frame(raw, path.suffix == ".dump", r["modeling_method"] == "MLIP_MACE")
            if r.get("simulation_time_basis") == "GPUMD_Time_fs_inferred_from_matching_header_schema":
                if header_time_fs(raw) / 1000 != r["simulation_time_ps"]:
                    raise ValueError("XYZ header time conversion mismatch")
            if len(atoms) != r["natoms"] or geometry_id(atoms) != r["geometry_sha256"]:
                raise ValueError("Exported geometry differs")
            _, _, value, status = structural_features(atoms)
            if value is None or abs(value - r["rdf_first_shell_peak_A"]) > 1e-10 or status != r["rdf_peak_status"]:
                raise ValueError("Recomputed peak differs")
            count += 1; n += 1
        if n != manifest["output_file_frames"][rel]:
            raise ValueError("Exported file count differs")
        save_json(OUT / "selection_progress.json", dict(stage="verify", verified_frames=count, total_frames=len(rows), latest_file=rel))
        print(f"verify {count}/{len(rows)} {rel}", flush=True)
    with (DEST / "metadata.csv").open(encoding="utf-8-sig") as handle:
        csv_rows = list(csv.DictReader(handle))
        if len(csv_rows) != len(rows):
            raise ValueError("CSV count differs")
        for actual, expected in zip(csv_rows, rows):
            if actual["structure_id"] != expected["structure_id"] or actual["modeling_method"] != expected["modeling_method"]:
                raise ValueError("CSV metadata differs")
            if expected["simulation_time_ps"] is None and actual["simulation_time_ps"] != "none":
                raise ValueError("CSV missing time must be none")
    if count != len(rows):
        raise ValueError("Missing selected frames")
    result = dict(selected=count, geometry_unique=True, original_frame_bytes="passed",
                  recomputed_first_shell_peaks="passed", error_formula="passed", time_labels="passed",
                  temperature_K=300, csv_json_counts="passed", cap="passed")
    save_json(DEST / "verification.json", result)
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["build", "verify"])
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if args.command == "build":
        build(args.workers)
    else:
        verify()


if __name__ == "__main__":
    main()
