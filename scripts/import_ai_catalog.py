#!/usr/bin/env python3
"""Import selected AI structures and real AIMD comparison frames into the catalog.

Runs inside pcm_database; original trajectories and selected raw frames remain
unchanged. References are ten distributed frames per trajectory, not new
independent catalog models. SOAP uses every atom and existing site parameters.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import json
import re
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from ase.io import read, write
from dscribe.descriptors import SOAP

from annotate_functionality import classify
from organize_ai_sources import ROOT, BASE, digest, save_json, save_csv
from select_ai_structures import raw_frames, parse_frame
from build_aimd_reference_frames import build_reference_frames, REFERENCE_POLICY

RAW = ROOT / "data/raw"
META = ROOT / "data/metadata"
DOCS = ROOT / "docs"
CACHE = ROOT / "data/.cache/ai_catalog"
SOAP_PARAMS = dict(r_cut=6.0, n_max=8, l_max=6, periodic=True, average="outer")


# AI 生成结构的贡献者按体系（排序后的元素拼接）指定；AI 结构统一按非晶入库。
AI_CONTRIBUTORS = {
    "Ge": "Shaojie Yuan", "Se": "Shaojie Yuan", "GeSe": "Shaojie Yuan",
    "GeSbTe": "Qundao Xu",
    "SbTe": "Siqi Tang", "GeTe": "Siqi Tang", "Sb": "Siqi Tang", "Te": "Siqi Tang",
}
AI_DEFAULT_PHASE = "amorphous"


def ai_system(formula):
    return "".join(sorted(set(re.findall(r"[A-Z][a-z]?", formula))))


def ai_contributor(formula, contributor=None):
    return contributor or AI_CONTRIBUTORS.get(ai_system(formula))


def family(method):
    return "MLIP" if method.startswith("MLIP_") else method


def soap_vector(soap, atoms):
    dummy = atoms.copy()
    dummy.set_atomic_numbers(np.full(len(dummy), 6))
    dummy.set_pbc(True)
    vector = np.asarray(soap.create(dummy), dtype=np.float64).reshape(-1)
    norm = np.linalg.norm(vector)
    if not np.isfinite(vector).all() or norm <= 0 or len(vector) != 252:
        raise ValueError("Invalid SOAP vector")
    return (vector / norm).astype(np.float32)


def import_task(task):
    rel, expected_sha, entries = task
    path = BASE / "AI_selected" / rel
    if digest(path) != expected_sha:
        raise ValueError(f"Selected structure source changed: {path}")
    signature = hashlib.sha256(json.dumps(dict(sha=expected_sha, soap=SOAP_PARAMS,
                                entries=[(r["structure_id"], r["structure_frame_0based"]) for r in entries]), sort_keys=True).encode()).hexdigest()
    cache = CACHE / (signature + ".npz")
    if cache.exists():
        with np.load(cache, allow_pickle=False) as data:
            sids = data["sids"].tolist()
            if all((RAW / "AI_generated" / sid / "POSCAR").exists() for sid in sids):
                return str(cache), len(entries), True
    soap = SOAP(species=[6], **SOAP_PARAMS)
    wanted = {r["structure_frame_0based"]: r for r in entries}
    sids, vectors, hashes = [], [], []
    for i, block in raw_frames(path):
        if i not in wanted:
            continue
        record = wanted[i]
        atoms = parse_frame(block, path.suffix == ".dump", record["modeling_method"] == "MLIP_MACE")
        sid = "ai_" + record["structure_id"].removeprefix("AI_")
        vectors.append(soap_vector(soap, atoms))
        atoms.wrap()
        dest = RAW / "AI_generated" / sid / "POSCAR"
        dest.parent.mkdir(parents=True, exist_ok=True)
        write(dest, atoms, format="vasp", direct=False, sort=True, vasp5=True)
        sids.append(sid); hashes.append(digest(dest))
    if len(sids) != len(entries):
        raise ValueError(f"Missing selected frames: {path}")
    cache.parent.mkdir(parents=True, exist_ok=True)
    with cache.with_suffix(".npz.tmp").open("wb") as handle:
        np.savez_compressed(handle, sids=np.array(sids), vectors=np.array(vectors), hashes=np.array(hashes))
    cache.with_suffix(".npz.tmp").replace(cache)
    return str(cache), len(sids), False


def base_fields(sid, formula, phase, natoms, density, rel_path, contributor, source, source_frame):
    return dict(sid=sid, **{"化学式": formula, "相": phase or AI_DEFAULT_PHASE, "温度": "300K",
        "体系原子数": natoms, "密度_Atom_per_Angstrom3": density, "结构文件": rel_path,
        "贡献者": ai_contributor(formula, contributor), "HPC来源路径": source, "代表帧来源": source_frame,
        "功能性": classify(formula, is_mp=False)})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    selected = json.loads((BASE / "AI_selected/metadata.json").read_text())
    selection = json.loads((BASE / "AI_selected/selection_manifest.json").read_text())
    records = json.loads((META / "pcm_metadata.json").read_text())
    original = [r for r in records if r.get("生成方式", "AIMD") != "AI"]
    for i, row in enumerate(original):
        row.setdefault("sid", f"s{i:04d}")
        row.update({"生成方式": "AIMD", "AI类别": None, "建模方法": "AIMD",
                    "g(r)主峰相对误差_pct": None, "collection_role": "catalog"})
    by_file = defaultdict(list)
    for r in selected:
        by_file[r["structure_file"]].append(r)
    tasks = [(rel, selection["output_file_sha256"][rel], rows) for rel, rows in sorted(by_file.items())]
    caches = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(import_task, task) for task in tasks]
        for future in as_completed(futures):
            cache, n, reused = future.result(); caches.append(cache)
            print(f"import SOAP {len(caches)}/{len(tasks)}: {n} structures cached={reused}", flush=True)
            save_json(CACHE / "progress.json", dict(files=len(caches), total_files=len(tasks)))
    descriptors, structure_hashes = {}, {}
    for path in caches:
        with np.load(path, allow_pickle=False) as data:
            for sid, vector, sha in zip(data["sids"], data["vectors"], data["hashes"]):
                descriptors[str(sid)] = vector; structure_hashes[str(sid)] = str(sha)
    ai_records, rich = [], []
    for r in selected:
        sid = "ai_" + r["structure_id"].removeprefix("AI_")
        rel_path = f"AI_generated/{sid}/POSCAR"
        row = base_fields(sid, r["material"], r["phase"], r["natoms"],
                          r["number_density_Atom_per_Angstrom3"], rel_path, r["contributor"],
                          "; ".join(r["source_hopper_paths"]) or r["source_file"], "AI-generated")
        row.update({"生成方式": "AI", "AI类别": family(r["modeling_method"]), "建模方法": r["modeling_method"],
                    "g(r)主峰相对误差_pct": r["peak_position_relative_error_pct"],
                    "g(r)主峰_A": r["rdf_first_shell_peak_A"], "AIMD主峰_A": r["aimd_first_shell_peak_A"],
                    "AIMD参考材料": r["reference_material"], "AIMD参考轨迹": r["reference_ids"],
                    "模拟时间_ps": r["simulation_time_ps"], "模拟时间数量级": r["simulation_time_order"],
                    "原报告温度_K": r["reported_temperature_K"], "源文件": r["source_file"],
                    "源帧": r["source_frame_0based"], "几何SHA256": r["geometry_sha256"],
                    "POSCAR_SHA256": structure_hashes[sid], "collection_role": "catalog"})
        ai_records.append(row)
        rich.append(dict(r, sid=sid, phase=row["相"], contributor=row["贡献者"], catalog_structure_file=rel_path,
                         catalog_POSCAR_sha256=structure_hashes[sid],
                         catalog_export="periodic_positions_wrapped;species_sorted;Cartesian_A"))
    soap = SOAP(species=[6], **SOAP_PARAMS)
    reference_records = build_reference_frames(soap, descriptors, structure_hashes)
    for r in original:
        path = RAW / r["结构文件"]
        descriptors[r["sid"]] = soap_vector(soap, read(path, format="vasp"))
        structure_hashes[r["sid"]] = digest(path)
    all_records = original + ai_records
    save_json(META / "pcm_metadata.json", all_records); save_csv(META / "pcm_metadata.csv", all_records)
    save_json(META / "ai_generated_metadata.json", rich); save_csv(META / "ai_generated_metadata.csv", rich)
    save_json(META / "aimd_reference_frames.json", reference_records); save_csv(META / "aimd_reference_frames.csv", reference_records)
    sids = [r["sid"] for r in all_records + reference_records]
    array = np.stack([descriptors[sid] for sid in sids])
    if len(set(sids)) != len(sids) or not np.isfinite(array).all():
        raise ValueError("Descriptor IDs/values invalid")
    payload = dict(meta=dict(descriptor="SOAP", species_mode="element-agnostic", average="outer",
                    normalized="l2-per-structure", **{k: v for k, v in SOAP_PARAMS.items() if k != "average"},
                    atomic_centers="all atoms, including large cells; no center subsampling",
                    note="AI and actual AIMD reference frames share SOAP parameters and joint PCA"),
                   sids=sids, dim=252, data=base64.b64encode(array.astype("<f4").tobytes()).decode(),
                   structure_sha256=structure_hashes)
    save_json(DOCS / "descriptors.json", payload)
    save_json(META / "ai_catalog_import.json", dict(catalog=len(all_records), original_AIMD=len(original), AI=len(ai_records),
               reference_frames=len(reference_records), reference_frames_per_trajectory=10, reference_policy=REFERENCE_POLICY,
               by_modeling_method=dict(Counter(r["建模方法"] for r in ai_records)), soap=payload["meta"]))
    print(f"Catalog={len(all_records)}: AIMD={len(original)}, AI={len(ai_records)}; reference overlays={len(reference_records)}", flush=True)


if __name__ == "__main__":
    main()
