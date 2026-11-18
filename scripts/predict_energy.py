#!/usr/bin/env python
"""用 MACE foundation model 对库内结构做单点能预测，记录 Energy per atom。

重要：这是 MLIP（机器学习势）预测能量，**不是 DFT 计算能量**。
它能较好反映结构间的相对能量，但不代表准确的绝对能量值。

环境：/opt/anaconda3/envs/mace
模型：data/foundation_MLIP_model/mace-omat-0-medium.model
输入：docs/structures/<sid>.poscar（代表帧，与 docs/data.json 的 sid 对应）
输出：data/metadata/energy.json
      {model, method:"MLIP", unit:"eV/atom", values:{sid: e_per_atom|null}}
      （增量写盘，可断点续跑；已算过的 sid 跳过）

用法（mace 环境）:
  /opt/anaconda3/envs/mace/bin/python scripts/predict_energy.py [--limit N] [--device mps|cpu|cuda] [--origin AI]

大批量（如 14,900 个 AI 结构）在 bwlab GPU 上跑：把本仓库的 scripts/predict_energy.py、docs/data.json、
docs/structures/、模型与 energy.json 按相同目录结构拷过去，CUDA_VISIBLE_DEVICES 固定一张空闲卡，
--device cuda 运行，结束后把 energy.json 拉回本地，再运行 build_site_data.py。
float32 单点能与 float64 相差 ~1e-6 eV/atom，与已有能量口径一致。
"""
import argparse
import json
import time
from pathlib import Path

from ase.io import read
from mace.calculators import MACECalculator

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "docs" / "data.json"
STRUCT_DIR = ROOT / "docs" / "structures"
MODEL = ROOT / "data" / "foundation_MLIP_model" / "mace-omat-0-medium.model"
OUT = ROOT / "data" / "metadata" / "energy.json"
MODEL_NAME = "mace-omat-0-medium"


def load_existing():
    if OUT.is_file():
        try:
            return json.loads(OUT.read_text())
        except Exception:
            pass
    return {"model": MODEL_NAME, "method": "MLIP", "unit": "eV/atom",
            "note": "MLIP-predicted energy (relative, not absolute DFT)", "values": {}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="只算前 N 个（0=全部）")
    ap.add_argument("--device", default="mps", choices=["mps", "cpu", "cuda"])
    ap.add_argument("--origin", default="all", choices=["all", "AI", "AIMD"], help="只算某一生成方式的结构")
    args = ap.parse_args()

    sids = [r["sid"] for r in json.loads(DATA.read_text())
            if args.origin == "all" or (r.get("generation_origin") or "AIMD") == args.origin]
    if args.limit:
        sids = sids[:args.limit]

    print(f"loading MACE model on {args.device} ...", flush=True)
    calc = MACECalculator(model_paths=[str(MODEL)], device=args.device, default_dtype="float32")

    store = load_existing()
    vals = store["values"]
    done0 = sum(1 for s in sids if s in vals)
    todo = [s for s in sids if s not in vals]
    print(f"total {len(sids)}, already {done0}, to compute {len(todo)}", flush=True)

    t0 = time.time()
    for k, sid in enumerate(todo, 1):
        p = STRUCT_DIR / f"{sid}.poscar"
        try:
            atoms = read(p, format="vasp")
            atoms.calc = calc
            e = float(atoms.get_potential_energy())
            vals[sid] = e / len(atoms)
        except Exception as exc:  # 记录失败但不中断
            vals[sid] = None
            print(f"  ! {sid} failed: {exc}", flush=True)
        if k % 20 == 0 or k == len(todo):
            OUT.write_text(json.dumps(store, ensure_ascii=False))
            rate = k / (time.time() - t0)
            eta = (len(todo) - k) / rate if rate else 0
            print(f"  {k}/{len(todo)}  {rate:.2f}/s  ETA {eta/60:.1f} min", flush=True)
    OUT.write_text(json.dumps(store, ensure_ascii=False))
    ok = sum(1 for v in vals.values() if v is not None)
    print(f"done: {ok}/{len(vals)} energies written to {OUT}", flush=True)


if __name__ == "__main__":
    main()
