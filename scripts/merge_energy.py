#!/usr/bin/env python
"""把 MLIP 预测能量(energy.json)合并进各 metadata 产物。

写入:
  data/metadata/pcm_metadata.json  —— 加 "能量_eV_per_atom_MLIP" 字段(按 sid=索引)
  data/metadata/pcm_metadata.csv   —— 由上面的 json 重新生成(附带新列)
  docs/data.json                   —— 每条记录加 energy_per_atom(按 sid 匹配)

注意:该能量为 MLIP 预测(mace-omat-0-medium),反映相对能量,非 DFT 绝对值。
"""
import json
import math
from pathlib import Path

from organize_ai_sources import save_csv

ROOT = Path(__file__).resolve().parent.parent
ENERGY = ROOT / "data" / "metadata" / "energy.json"
META_JSON = ROOT / "data" / "metadata" / "pcm_metadata.json"
META_CSV = ROOT / "data" / "metadata" / "pcm_metadata.csv"
DATA = ROOT / "docs" / "data.json"
COL = "能量_eV_per_atom_MLIP"


def main():
    raw = json.loads(ENERGY.read_text())["values"]
    # 非有限值(inf/nan，结构异常导致)统一记为 null，避免污染 JSON / 颜色映射
    vals = {k: (v if (v is not None and math.isfinite(v)) else None) for k, v in raw.items()}
    dropped = [k for k, v in raw.items() if v is not None and not math.isfinite(v)]
    if dropped:
        print(f"coerced {len(dropped)} non-finite energies to null: {dropped}")

    # pcm_metadata.json：按记录自身 sid 匹配（AI 结构为 ai_*；旧 AIMD 记录无 sid 时回退 s{index:04d}）
    meta = json.loads(META_JSON.read_text(encoding="utf-8"))
    for i, r in enumerate(meta):
        r[COL] = vals.get(r.get("sid", f"s{i:04d}"))
    META_JSON.write_text(json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    # CSV 镜像 json 列（与其他 metadata 脚本同一写法）
    save_csv(META_CSV, meta)

    # docs/data.json (按 sid)
    data = json.loads(DATA.read_text(encoding="utf-8"))
    for r in data:
        r["energy_per_atom"] = vals.get(r["sid"])
    DATA.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")

    n = sum(1 for r in data if r.get("energy_per_atom") is not None)
    print(f"merged energy into metadata/csv/data.json: {n}/{len(data)} have values")


if __name__ == "__main__":
    main()
