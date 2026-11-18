#!/usr/bin/env python3
"""Export ten comparison frames per AIMD trajectory without changing RDF baselines.

Sampling spans the existing RDF reference window and retains its former midpoint
representative and SID. Catalog structures and their SOAP vectors stay unchanged.
Run in pcm_database before build_site_data.py and build_asedb.py.
"""
import argparse
import base64
import json
from collections import Counter

import numpy as np
from ase.io import write
from dscribe.descriptors import SOAP

from annotate_functionality import classify
from build_descriptors import SOAP_PARAMS
from organize_ai_sources import ROOT, BASE, REF, digest, save_csv, save_json, xdat_header, xdat_frames

RAW = ROOT / "data/raw"
META = ROOT / "data/metadata"
DOCS = ROOT / "docs"
REFERENCE_POLICY = (
    "10 distributed frames per trajectory from the existing 200 RDF samples, "
    "including the prior midpoint; map overlay excluded from catalog counts"
)


def reference_indices(samples, count=10):
    if count < 10:
        raise ValueError("At least ten reference frames per trajectory are required")
    samples = sorted(set(samples))
    if len(samples) < count:
        raise ValueError(f"Only {len(samples)} distinct RDF samples available; requested {count}")
    positions = np.rint(np.linspace(0, len(samples) - 1, count)).astype(int)
    middle = len(samples) // 2
    if middle not in positions:
        positions[np.argmin(abs(positions - middle))] = middle
    return [samples[i] for i in sorted(positions)]


def build_reference_frames(soap, descriptors, structure_hashes, count=10):
    plan = json.loads((REF / "archive_manifest.json").read_text())
    rdf = json.loads((REF / "rdf_summary.json").read_text())
    lookup = {r["reference_id"]: r for r in rdf["trajectories"]}
    records = []
    for number, r in enumerate(plan["records"], 1):
        samples = lookup[r["reference_id"]]["sampled_indices"]
        indices = reference_indices(samples, count)
        wanted = set(indices)
        midpoint = samples[len(samples) // 2]
        path = BASE / r["archive_file"]
        atoms, _ = xdat_header(path)
        seen = set()
        for frame, positions in xdat_frames(path):
            if frame not in wanted:
                continue
            atoms.set_positions(positions)
            atoms.wrap()
            base_sid = "ref_" + r["reference_id"].removeprefix("A_")
            sid = base_sid if frame == midpoint else f"{base_sid}_f{frame:06d}"
            dest = RAW / "AIMD_reference" / sid / "POSCAR"
            dest.parent.mkdir(parents=True, exist_ok=True)
            write(dest, atoms, format="vasp", direct=False, sort=True, vasp5=True)
            dummy = atoms.copy()
            dummy.set_atomic_numbers(np.full(len(dummy), 6))
            dummy.set_pbc(True)
            vector = np.asarray(soap.create(dummy), dtype=np.float64).reshape(-1)
            norm = np.linalg.norm(vector)
            if len(vector) != 252 or not np.isfinite(vector).all() or norm <= 0:
                raise ValueError(f"Invalid reference SOAP: {sid}")
            descriptors[sid] = (vector / norm).astype(np.float32)
            structure_hashes[sid] = digest(dest)
            records.append({
                "sid": sid, "化学式": r["material"], "相": r["phase"] or "amorphous", "温度": "300K",
                "体系原子数": len(atoms), "密度_Atom_per_Angstrom3": len(atoms) / atoms.get_volume(),
                "结构文件": str(dest.relative_to(RAW)), "贡献者": r["contributor"],
                "HPC来源路径": r["archive_file"], "代表帧来源": "AIMD-reference-frame",
                "功能性": classify(r["material"], is_mp=False), "生成方式": "AIMD", "AI类别": None,
                "建模方法": "AIMD", "g(r)主峰相对误差_pct": None,
                "reference_id": r["reference_id"], "reference_material": r["material"],
                "reference_frame_0based": frame, "reference_trajectory_sha256": r["sha256"],
                "reference_frames_per_trajectory": count,
                "reference_role": "distributed_frames_from_RDF_samples;retains_previous_midpoint",
                "POSCAR_SHA256": structure_hashes[sid], "collection_role": "reference",
            })
            seen.add(frame)
            if seen == wanted:
                break
        if seen != wanted:
            raise ValueError(f"Missing reference frames in {path}: {wanted - seen}")
        print(f"AIMD comparison {number}/{len(plan['records'])}: {r['material']} × {count}", flush=True)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frames-per-trajectory", type=int, default=10)
    args = parser.parse_args()
    catalog = json.loads((META / "pcm_metadata.json").read_text())
    payload = json.loads((DOCS / "descriptors.json").read_text())
    if payload.get("dim") != 252 or any(payload["meta"].get(k) != v for k, v in SOAP_PARAMS.items()):
        raise ValueError("Existing catalog SOAP parameters differ; rebuild the catalog first")
    vectors = np.frombuffer(base64.b64decode(payload["data"]), dtype="<f4").reshape(-1, 252)
    index = {sid: i for i, sid in enumerate(payload["sids"])}
    descriptors, hashes = {}, {}
    for r in catalog:
        sid = r["sid"]
        sha = digest(RAW / r["结构文件"])
        if sha != payload.get("structure_sha256", {}).get(sid):
            raise ValueError(f"Catalog POSCAR changed; rebuild its descriptor first: {sid}")
        descriptors[sid] = vectors[index[sid]]
        hashes[sid] = sha
    soap = SOAP(species=[6], **SOAP_PARAMS)
    references = build_reference_frames(soap, descriptors, hashes, args.frames_per_trajectory)
    sids = [r["sid"] for r in catalog + references]
    if len(set(sids)) != len(sids):
        raise ValueError("Duplicate comparison frame IDs")
    payload.update(sids=sids, structure_sha256=hashes,
                   data=base64.b64encode(np.stack([descriptors[sid] for sid in sids]).astype("<f4").tobytes()).decode())
    payload["meta"]["note"] = "Catalog and ten actual AIMD frames per trajectory share SOAP parameters and joint PCA"
    save_json(META / "aimd_reference_frames.json", references)
    save_csv(META / "aimd_reference_frames.csv", references)
    save_json(DOCS / "descriptors.json", payload)
    summary_path = META / "ai_catalog_import.json"
    summary = json.loads(summary_path.read_text())
    summary.update(reference_frames=len(references), reference_frames_per_trajectory=args.frames_per_trajectory,
                   reference_policy=REFERENCE_POLICY)
    save_json(summary_path, summary)
    print(f"Catalog unchanged: {len(catalog)}; AIMD reference frames: {len(references)}; "
          f"trajectory counts: {dict(Counter(Counter(r['reference_id'] for r in references).values()))}")


if __name__ == "__main__":
    main()
