#!/usr/bin/env python3
"""由 pcm_metadata.json 生成前端站点数据（GitHub Pages 用）。

产出：
  docs/data.json          —— 记录数组（含 system/elements/sid）
  docs/structures/<sid>.poscar —— 每条记录的代表结构文件副本（供 3D 预览/下载）
"""
import json
import math
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
META = ROOT / "data" / "metadata" / "pcm_metadata.json"
ENERGY = ROOT / "data" / "metadata" / "energy.json"   # MLIP 预测能量(可选)
DOCS = ROOT / "docs"
OUT = DOCS / "data.json"
STRUCT_DIR = DOCS / "structures"
RAW = ROOT / "data" / "raw"


def site_record(r, index, energy):
    sid = r.get("sid", f"s{index:04d}")
    els = elements(r["化学式"])
    value = energy.get(sid)
    return {
        "sid": sid, "formula": r["化学式"], "system": "".join(els), "elements": els,
        "phase": r["相"] or "unknown", "temperature": r["温度"], "natoms": r["体系原子数"],
        "density": r["密度_Atom_per_Angstrom3"], "contributor": r.get("贡献者"),
        "path": r["结构文件"], "hpc_source": r.get("HPC来源路径", ""),
        "source_frame": r.get("代表帧来源", ""), "functionality": r.get("功能性", "Others"),
        "energy_per_atom": value if value is not None and math.isfinite(value) else None,
        "generation_origin": r.get("生成方式", "AIMD"), "generation_family": r.get("AI类别") or "AIMD",
        "modeling_method": r.get("建模方法", "AIMD"), "rdf_peak_error_pct": r.get("g(r)主峰相对误差_pct"),
        "rdf_peak_A": r.get("g(r)主峰_A"), "aimd_peak_A": r.get("AIMD主峰_A"),
        "reference_material": r.get("AIMD参考材料", r.get("reference_material")),
        "reference_ids": r.get("AIMD参考轨迹", []), "reference_id": r.get("reference_id"),
        "reference_frame_0based": r.get("reference_frame_0based"),
        "reference_frames_per_trajectory": r.get("reference_frames_per_trajectory"),
        "simulation_time_ps": r.get("模拟时间_ps"), "simulation_time_order": r.get("模拟时间数量级", "none"),
        "reported_temperature_K": r.get("原报告温度_K"), "source_file": r.get("源文件"),
        "source_frame_0based": r.get("源帧"), "collection_role": r.get("collection_role", "catalog"),
        "structure_url": f"structures/{sid}.poscar",
    }


def elements(formula: str):
    return sorted(set(re.findall(r"[A-Z][a-z]?", formula)))


def main():
    recs = json.loads(META.read_text(encoding="utf-8"))
    energy = json.loads(ENERGY.read_text())["values"] if ENERGY.is_file() else {}
    STRUCT_DIR.mkdir(parents=True, exist_ok=True)

    out = []
    for i, r in enumerate(recs):
        sid = r.get("sid", f"s{i:04d}")
        src = RAW / r["结构文件"]
        if not src.is_file():
            raise FileNotFoundError(src)
        shutil.copy2(src, STRUCT_DIR / f"{sid}.poscar")
        out.append(site_record(r, i, energy))
    OUT.write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    references = []
    ref_meta = META.parent / "aimd_reference_frames.json"
    if ref_meta.exists():
        for i, r in enumerate(json.loads(ref_meta.read_text())):
            shutil.copy2(RAW / r["结构文件"], STRUCT_DIR / f"{r['sid']}.poscar")
            references.append(site_record(r, i, {}))
    (DOCS / "aimd_references.json").write_text(json.dumps(references, ensure_ascii=False, separators=(",", ":")))
    n_struct = len(list(STRUCT_DIR.glob("*.poscar")))
    print(f"data.json: {len(out)} 条 ({OUT.stat().st_size//1024} KB)")
    print(f"structures/: {n_struct} 个结构文件")


if __name__ == "__main__":
    main()
