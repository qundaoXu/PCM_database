#!/usr/bin/env python3
"""就地更新已入库 AI 结构的相与贡献者，无需重算 SOAP。

- AI 生成结构：未标注相 → amorphous；贡献者按体系取 import_ai_catalog.AI_CONTRIBUTORS。
- AIMD 参考帧（仅用于分布图对照）：未标注相 → amorphous。

更新 data/metadata 下的 pcm_metadata / ai_generated_metadata / aimd_reference_frames
（JSON + CSV）。之后运行 scripts/build_site_data.py 重新生成 docs/data.json。
"""
import json
from collections import Counter

from import_ai_catalog import AI_CONTRIBUTORS, AI_DEFAULT_PHASE, META, ai_system
from organize_ai_sources import save_csv, save_json

UNLABELED = (None, "", "unknown")


def relabel(rows, phase_key, contributor_key=None, formula_key=None):
    changed = Counter()
    for row in rows:
        if row.get(phase_key) in UNLABELED:
            row[phase_key] = AI_DEFAULT_PHASE
            changed["phase"] += 1
        if contributor_key:
            system = ai_system(row[formula_key])
            if system not in AI_CONTRIBUTORS:
                raise KeyError(f"No contributor mapping for AI system {system!r} ({row.get('sid')})")
            if row.get(contributor_key) != AI_CONTRIBUTORS[system]:
                row[contributor_key] = AI_CONTRIBUTORS[system]
                changed["contributor"] += 1
    return changed


def update(name, select, **keys):
    path = META / f"{name}.json"
    rows = json.loads(path.read_text(encoding="utf-8"))
    changed = relabel([r for r in rows if select(r)], **keys)
    save_json(path, rows)
    save_csv(META / f"{name}.csv", rows)
    print(f"{name}: {dict(changed) or 'no change'}")


def main():
    update("pcm_metadata", lambda r: r.get("生成方式") == "AI",
           phase_key="相", contributor_key="贡献者", formula_key="化学式")
    update("ai_generated_metadata", lambda r: True,
           phase_key="phase", contributor_key="contributor", formula_key="material")
    update("aimd_reference_frames", lambda r: True, phase_key="相")


if __name__ == "__main__":
    main()
