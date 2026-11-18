#!/usr/bin/env python3
"""Produce local Chinese collection report from verified inventories."""
import base64
import io
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from ase.io import read

from organize_ai_sources import BASE, REF, OUT, composition, formula, geometry_id, save_csv, save_json


def main():
    plan = json.loads((REF / "archive_manifest.json").read_text())
    rdfs = json.loads((REF / "rdf_summary.json").read_text())
    local = json.loads((OUT / "ai_local_summary.json").read_text())
    local_files = json.loads((OUT / "ai_local_files.json").read_text())
    frames = json.loads((OUT / "ai_local_frames.json").read_text())
    policy_path = OUT / "collection_policy.json"
    policy = json.loads(policy_path.read_text()) if policy_path.exists() else {"excluded_sources": []}
    approved_300K = policy.get("temperature_policy", {}).get("reference_temperature_K") == 300
    selection_path = BASE / "AI_selected/selection_manifest.json"
    selection = json.loads(selection_path.read_text()) if selection_path.exists() else None
    selection_text = ""
    if selection:
        selected_metadata = json.loads((BASE / "AI_selected/metadata.json").read_text())
        method_table = ["| 建模方式 | 筛选后结构 | 峰位误差中位数 / % | 峰位误差 P95 / % |", "|---|---:|---:|---:|"]
        for method, n in sorted(selection["by_modeling_method"].items()):
            errors = [r["peak_position_relative_error_pct"] for r in selected_metadata if r["modeling_method"] == method]
            method_table.append(f"| {method} | {n:,} | {np.median(errors):.4f} | {np.percentile(errors,95):.4f} |")
        selection_text = ("## 多样性筛选与正式 metadata\n\n"
            f"从 {selection['candidates']:,} 个几何去重候选选择 **{selection['selected']:,} 个结构**，上限 {selection['cap']:,}。Flow Matching 已纳入。\n\n"
            + "\n".join(method_table) + "\n\n"
            "按参考材料、建模方式与实际原子数分层，用候选数平方根分配预算并保留稀有尺寸。\n"
            "用元素分辨的径向分布、配位数分布及密度表示径向和配位多样性，缩放、PCA 后最远点采样。\n"
            "已知时间的老化轨迹按时间数量级分层，每层按物理时间均匀选点。此轮不以峰位误差大小筛选。\n\n"
            "AIMD 和对应评价温度按用户约定为 300 K，原报告温度另行保留。主峰误差逐结构重新计算。\n"
            f"{len(selection['ineligible']):,} 个候选因参考或内部主峰无法定位而未进入可评价集合；\n"
            f"已选结构峰位状态为 {json.dumps(selection['peak_statuses'],ensure_ascii=False)}。\n"
            f"模拟时间数量级分布为 {json.dumps(selection['by_simulation_time_order'],ensure_ascii=False)}；缺失者为 none。\n\n"
            "产物：`AI_selected/metadata.json` / `metadata.csv`、`structures/`、`selection_manifest.json`、\n"
            "`selection_strata.csv` 和 `verification.json`（完整验证后生成）。结构帧保留原始字节，来源哈希可追溯。\n"
            "每个 metadata 记录建模方式、模拟时间/数量级、300 K 参考、主峰位置及百分比误差。\n\n")
    local_by_sha = defaultdict(list)
    for r in local_files:
        local_by_sha[r["sha256"]].append(r["file"])
    remote = json.loads((OUT / "hopper_counts.json").read_text())["files"]
    samples = {r["path"]: r for r in json.loads((OUT / "hopper_first_frames.json").read_text())["files"]}
    context = json.loads((OUT / "hopper_context.json").read_text())
    by_sha = defaultdict(list)
    for r in remote:
        if r["format"] == "lammps-dump-text":
            by_sha[r["sha256"]].append(r)
    candidates = []
    for sha, copies in by_sha.items():
        row = min(copies, key=lambda r: ("/eval_staging_" in r["path"], r["path"]))
        atoms = read(io.StringIO(base64.b64decode(samples[row["path"]]["first_frame_base64"]).decode()),
                     format="lammps-dump-text", specorder=["Ge", "Sb", "Te"])
        if not np.isfinite(atoms.positions).all() or not atoms.pbc.all() or atoms.get_volume() <= 0:
            raise ValueError(f"Invalid representative frame: {row['path']}")
        key = composition(atoms)
        if key != "Ge:2;Sb:2;Te:5":
            raise ValueError(f"Unexpected remote composition: {row['path']}: {key}")
        category = ("aging_373K" if "/1-GST225-aging/" in row["path"] else
                    "liquid_benchmark" if "/canon_liquid_" in row["path"] else "amorphous_benchmark")
        selected = np.linspace(0, row["frames"] - 1, min(400, row["frames"]), dtype=int).tolist()
        candidates.append(dict(path=row["path"], source_host="hopper.nus.edu.sg", sha256=sha,
                               category=category, composition=key, material=formula(key),
                               natoms=len(atoms), stored_complete_frames=row["frames"],
                               first_timestep=row["first_timestep"], last_timestep=row["last_timestep"],
                               volume_first_frame_A3=atoms.get_volume(),
                               density_first_frame_Atom_per_Angstrom3=len(atoms) / atoms.get_volume(),
                               duplicate_file_paths=[r["path"] for r in copies if r is not row],
                               local_file_copies=local_by_sha.get(sha, []),
                               already_in_local_inventory=bool(local_by_sha.get(sha)),
                               status=("already_available_locally_do_not_count_twice" if local_by_sha.get(sha)
                                       else "candidate_not_bulk_downloaded_or_frame_deduplicated"),
                               count_validation="all_blocks_complete; first_frame_ASE_geometry_and_GST225_composition_checked",
                               max400_stored_index_plan=selected,
                               sampling_note=("Nonuniform original timestamps: select by physical time before collection"
                                              if "run_10ns/" in row["path"] else "Stored-index plan only; equilibrium/time criteria still to be chosen")))
    candidates.sort(key=lambda r: r["path"])
    save_json(OUT / "hopper_candidates.json", candidates)
    save_csv(OUT / "hopper_candidates.csv", candidates)
    categories = defaultdict(list)
    for r in candidates:
        categories[r["category"]].append(r)
    remote_summary = dict(unique_trajectory_files=len(candidates),
                          stored_complete_frames=sum(r["stored_complete_frames"] for r in candidates),
                          raw_trajectory_paths=sum(len(r) for r in by_sha.values()),
                          excluded_duplicate_file_paths=sum(len(r) - 1 for r in by_sha.values()),
                          by_category={k: dict(files=len(rows), frames=sum(r["stored_complete_frames"] for r in rows),
                                               atom_counts=sorted({r["natoms"] for r in rows}))
                                       for k, rows in categories.items()},
                          max400_stored_index_plan_frames=sum(len(r["max400_stored_index_plan"]) for r in candidates),
                          trajectories_already_local=sum(r["already_in_local_inventory"] for r in candidates),
                          frames_already_local=sum(r["stored_complete_frames"] for r in candidates if r["already_in_local_inventory"]),
                          additional_remote_candidate_frames=sum(r["stored_complete_frames"] for r in candidates if not r["already_in_local_inventory"]),
                          status="candidate_frames_not_unique_validated_structures",
                          note="File SHA256 deduplication only; shared initial MD frames may still repeat")
    save_json(OUT / "hopper_summary.json", remote_summary)

    deep_rows, deep_seen = [], {}
    for path, data in context["poscars_base64"].items():
        atoms = read(io.StringIO(base64.b64decode(data).decode()), format="vasp")
        geom = geometry_id(atoms)
        lineage = "1080atoms" if "/1080atoms/" in path else "216atoms/" + path.split("structure_")[1].split("/")[0]
        excluded = any(e["source_host"] == "hopper.nus.edu.sg" and path.startswith(e["root"] + "/")
                       for e in policy["excluded_sources"])
        deep_rows.append(dict(path=path, material=formula(composition(atoms)), natoms=len(atoms),
                              geometry_sha256=geom, duplicate_of=deep_seen.get(geom), model_lineage=lineage,
                              volume_A3=atoms.get_volume(), included=not excluded,
                              status="excluded_by_user" if excluded else "static_structure_generation_origin_unconfirmed"))
        deep_seen.setdefault(geom, path)
    save_json(OUT / "deeph_static_structures.json", deep_rows)
    save_csv(OUT / "deeph_static_structures.csv", deep_rows)

    reference_counts = Counter(r["material"] for r in plan["records"])
    reference_frames = Counter()
    for r in plan["records"]:
        reference_frames[r["material"]] += r["frames"]
    pooled = {r["material"]: r for r in rdfs["materials"]}
    ref_table = ["| 材料 | 去重轨迹 | 原始帧数 | 平均 g(r) 第一主峰 / Å |", "|---|---:|---:|---:|"]
    for material, n in sorted(reference_counts.items()):
        ref_table.append(f"| {material} | {n} | {reference_frames[material]:,} | {pooled[material]['reference_peak_A']:.6f} |")

    local_table = ["| 材料（实际组分） | MLIP-MD | Diffusion | Flow Matching |", "|---|---:|---:|---:|"]
    by_material = defaultdict(Counter)
    for r in frames:
        if not r["duplicate_of"]:
            by_material[r["material"]][r["method"]] += 1
    for mat, counts in sorted(by_material.items()):
        local_table.append(f"| {mat} | {counts['MLIP-MD']:,} | {counts['diffusion']:,} | {counts['flow_matching']:,} |")
    exact_requested = sum(r["method"] in {"MLIP-MD", "diffusion"} and not r["duplicate_of"]
                          and r["exact_composition_reference"] for r in frames)
    requested_count = sum(local["unique_by_method"].get(k, 0) for k in ["MLIP-MD", "diffusion"])
    gst_files = [r for r in local_files if r["file"].startswith("GST225/")]
    gst_frames = [r for r in frames if r["file"].startswith("GST225/")]
    gst_summary = dict(files=len(gst_files), complete_valid_frames=len(gst_frames),
                       unique_geometries=sum(not r["duplicate_of"] for r in gst_frames),
                       duplicate_frames=sum(bool(r["duplicate_of"]) for r in gst_frames),
                       by_actual_natoms=dict(Counter(r["natoms"] for r in gst_frames)),
                       filename_atom_count_mismatches=[r["file"] for r in gst_files if r["filename_atom_count_matches"] is False],
                       identical_hopper_files={r["file"]: r["identical_hopper_paths"] for r in gst_files if r["identical_hopper_paths"]},
                       max400_per_file_plan=[dict(file=r["file"], frame_indices_0based=np.linspace(
                           0, r["frames"] - 1, min(400, r["frames"]), dtype=int).tolist()) for r in gst_files],
                       sampling_status=("selected_in_AI_selected;original_frames_preserved" if selection else
                                        "count_plan_only_not_executed; nonuniform aging trajectory requires physical_time_selection"),
                       all_frames_unified_rdf_status="computed_against_300K_material_reference" if selection else "pending")
    save_json(OUT / "gst225_summary.json", gst_summary)
    save_json(OUT / "gst225_files.json", gst_files)
    save_csv(OUT / "gst225_files.csv", gst_files)
    gst_table = ["| 本地原文件名 | 实际原子数 | 完整帧 | 去重后 | timestep 最小–最大 |", "|---|---:|---:|---:|---|"]
    for r in gst_files:
        gst_table.append(f"| {Path(r['file']).name} | {', '.join(map(str,r['atom_counts']))} | {r['frames']:,} | {r['unique_geometries']:,} | {r['min_timestep']}–{r['max_timestep']} |")
    gst_table_text = "\n".join(gst_table)
    (BASE / "GST225/README.md").write_text("# GST225 本地结构清点\n\n" + gst_table_text + "\n\n"
        "原始文件名和字节保持不变。文件名 512/720 的实际原子数分别为 504/8001。\n"
        "720 文件与 Hopper 10 ns 老化轨迹的 SHA256 相同，只计一份来源。\n"
        "216 原子文件的 9 帧 timestep 不单调，作为快照集合登记。\n"
        "详情见 [清点报告](../design.md)，逐文件登记见 [CSV](../inventory/gst225_files.csv)。\n"
        "Hopper deeph 中此前 216/1080 静态结构已按用户要求排除。\n")
    reference_readme = """# AIMD 参考轨迹归档

物理文件统一位于 `<材料>/XDATCAR*`。原始文件名保留，材料名称由 ASE 读取真实元素及其约分比例确定。
GST225、GST467 分别使用 `Ge2Sb2Te5/`、`Ge4Sb6Te7/`。

81 个原路径对应 79 个文件 SHA256 唯一的轨迹，覆盖 24 种材料、152,985 个完整帧。
重复的 Sb、GeSe 文件只保留一份实体；原数据包路径改为相对符号链接，仍可运行原验证脚本。
`archive_manifest.json`、`references.json`、`references.csv` 保存来源路径、SHA256、组分、帧数及旧参考条件。
`XDATCAR-GeTe2` 的真实组分是 GeTe，已经放入 `GeTe/`，不能由文件名推断为 GeTe₂。

## 平均 g(r) 与第一峰

每条去重轨迹末尾 1000 个存储帧中均匀选择 200 帧（索引从 0 起），先平均帧的 g(r)，
再对同材料的各条轨迹均值等权平均。每条轨迹的采样索引、独立峰位、前后半段漂移均有记录。
共计算 15,800 个参考帧；不是对全部 152,985 帧平均。以 0.05 Å 分箱，周期最小镜像、有限 N
原子对数归一，r 最大 3.5 Å；在 2.0–3.3 Å 第一配位壳定位主峰，5 点滑动平均与三点抛物线插值。
该峰是总数目 g(r) 的第一壳主峰，不是衍射结构因子的主峰。

每个材料目录包含 `mean_rdf.json` / `mean_rdf.csv`；全局汇总为 `rdf_summary.json` / `rdf_summary.csv`。
全部 24 个材料平均峰及 79 条独立轨迹峰均能明确定位，独立参考前后半段峰位漂移均小于 5%。

按用户约定，全部 AIMD 参考统一按 300 K 对应评价，`temperature_basis=user_assumed_300K`。
原来 34 条报告 300 K、45 条未报告温度的信息保存在 reported_temperature_K 和原始来源记录中。
同材料各条轨迹等权合并，MLIP 原报告约 300–400 K 的数据统一与此 300 K 基准计算误差。
原 Ge–Se 数据包将 legacy 与已知 300 K 分组计算的历史误差保持原样，不能当作新 pooled 基准的误差。

""" + "\n".join(ref_table) + "\n"
    (REF / "README.md").write_text(reference_readme)

    remote_table = ["| Hopper 来源 | 去重轨迹文件 | 候选完整帧 | 每帧原子数 |", "|---|---:|---:|---:|"]
    names = {"amorphous_benchmark": "MLIP_benchmark：amorphous", "liquid_benchmark": "MLIP_benchmark：liquid", "aging_373K": "GST225：373 K 长时间老化"}
    for name, stats in remote_summary["by_category"].items():
        remote_table.append(f"| {names[name]} | {stats['files']} | {stats['frames']:,} | {', '.join(map(str,stats['atom_counts']))} |")
    local_table_text = "\n".join(local_table)
    remote_table_text = "\n".join(remote_table)
    report = f"""# AI 生成结构模块：归档、多样性筛选与 metadata

日期：2026-11-17（Asia/Singapore）。本轮完成本地 AIMD 归档、参考第一壳 g(r) 均值、
本地 AI 数据逐帧核验、Hopper 指定目录只读清点，以及本地结构的多样性筛选、建模方式 metadata 和 300 K 主峰误差计算。
前端模块留待下一阶段实现。

## 已完成的 AIMD 归档

`AIMD_reference/<材料>/` 已归档 **24 种材料、79 条 SHA256 去重轨迹、152,985 帧**。
全部原始字节按 SHA256 核验；81 个原路径中重复的 Sb、GeSe 各合并一份。
原两个数据包的参考文件位置保留相对链接，不改变结构字节、原索引或历史筛选结果。
`XDATCAR-GeTe2` 实际为 GeTe；`sb2te2` 实际为 SbTe，均按内容归档。
来源 JSON/CSV 与各材料的平均 g(r) 已生成，完整口径见 `AIMD_reference/README.md`。

每轨迹末尾 1000 个存储帧均匀选 200 帧，先算轨迹平均再按轨迹等权合并，共 15,800 个参考采样帧。
按用户约定，全部 AIMD 参考按 300 K 评价；原报告温度保存在 reported_temperature_K，不改变原始结构。
本轮复算精确复现 Ge–Se 原包的 44 条独立参考峰，24 个材料平均峰均成功定位。

## 本地 AI 数量（逐帧 ASE 核验）

| 类型 | 原始完整帧 | 去重后结构 |
|---|---:|---:|
| MLIP-MD：GeTe、Sb–Te、GST225 | {local['by_method']['MLIP-MD']:,} | {local['unique_by_method']['MLIP-MD']:,} |
| Diffusion：Ge–Se（含纯 Ge / Se） | {local['by_method']['diffusion']:,} | {local['unique_by_method']['diffusion']:,} |
| Flow Matching：Ge–Se，单独标注 | {local['by_method']['flow_matching']:,} | {local['unique_by_method']['flow_matching']:,} |
| 合计 | {local['complete_valid_frames']:,} | **{local['unique_geometries']:,}** |

23 个 MLIP 文件各 400 帧；34 个 Ge–Se 生成文件共 14,085 帧，与原 `structures.csv`
中的组分、原子数、晶胞体积和全部几何哈希逐帧一致。GeTe 有两条重复初始帧；
新增 GST225 文件另做全量核验，逐文件数量与重复项见下表及帧级索引。
去重对元素、笛卡尔坐标和晶胞以 1e-6 Å 量化后哈希，忽略原子排列顺序；未做平移、旋转或结构相似性聚类。
MD 存储帧数量不代表统计独立样本数量。

只计本次明确提出的 **MLIP-MD + Diffusion：{requested_count:,} 个去重结构**，已经达到 10,000+ 的数量目标。
其中有严格化学计量 AIMD 对应者 **{exact_requested:,} 个**，另有取整配比结构按下述容差匹配参考。
Flow Matching 的 1,350 个结构另外保留，不能写成 Diffusion。

800 帧暂无严格组分参考：`sb2te_1000` 为 Sb667Te333，`sbte2_1000` 为 Sb333Te667，
每个文件 400 帧。它们与理想 Sb₂Te、SbTe₂ 只有取整差异，但不能静默改组分或强配参考。
本轮按最大元素原子分数偏差 0.001（0.1 百分点）容差，匹配最近理想组分参考，
并在 metadata 中保留真实组分、参考材料、偏差和 `integer_stoichiometry_rounding` 标记。

原 Ge–Se 数据包已经按原参考组、第一峰误差 <5% 预筛选。历史误差保留在
`legacy_peak_error_pct` 中；本轮统一计算对新合并参考的误差，不把 <5% 设为新模块的准入门槛。

{local_table_text}

## 新增本地 GST225 多尺寸轨迹

按用户补充，将 `reorg/GST225/` 作为本地 GST225 轨迹来源。共 **{gst_summary['files']} 个文件、
{gst_summary['complete_valid_frames']:,} 个完整帧、{gst_summary['unique_geometries']:,} 个去重几何结构**。
实际组分由 ASE 读取，全部逐帧检查有限坐标、有效周期晶胞、原子数、几何哈希与重复项。
大体积 dump 按完整帧逐块交给 ASE 解析，避免一次将整条轨迹载入内存。

{gst_table_text}

两个文件名与内容不符，保留原文件并将真实值登记到 JSON/CSV：

- `GST225_512atoms.dump` 实际是 **504 原子**。
- `GST225_720atoms.dump` 实际是 **8,001 原子**，SHA256 与 Hopper 的
  `run_10ns/GST_aging_10ns.dump` 完全一致，作为同一来源的本地副本，不在跨来源数量中重复累计。

匹配 10 ns 轨迹的帧可由原模拟输入确认 373 K、dt=0.002 ps；其余文件的原报告温度和 dt 未提供，保留为空，
不根据文件名或帧数推算物理时长。原报告温度保留在 reported_temperature_K，评价温度统一为 300 K。
各文件原始 timestep 已保存。当前原始轨迹完整保留，筛选后的独立副本见 AI_selected。
`gst225_summary.json` 中最多 400 帧/文件的索引是数量规划，非均匀老化轨迹仍须按物理时间选点。
216 原子的 XYZ 共 9 帧，其 timestep 不单调，登记为快照集合而非时间连续的轨迹；不由首末帧推算时长。

## Hopper 的可采样数据

只读检查三个指定目录，没有删除或修改远端文件。排除训练/验证数据、DFT 标注轨迹、
软件示例以及同一轨迹的评测目录副本后，得到 **25 条去重轨迹文件、18,629 个完整候选帧**。
30 个原轨迹路径中有 5 个完全相同的文件副本。每条轨迹核对所有帧块完整性，并用 ASE
读取首帧，确认其真实组分为 GST225、周期晶胞有效。远端尚未逐帧做几何去重、平衡段选择或质量验收，
因此 18,629 是候选帧数，不能与本地去重后数量直接相加后宣称独立结构数。
其中 **{remote_summary['trajectories_already_local']} 条 / {remote_summary['frames_already_local']:,} 帧已存在本地 GST225 文件中**；
排除这些完全相同的文件后，远端还有 **{remote_summary['additional_remote_candidate_frames']:,} 个额外候选帧**。

{remote_table_text}

- `run_10ns/GST_aging_10ns.dump`：8,001 原子，2,108 帧，step 0–5,003,000，dt=0.002 ps，
  存储覆盖 0–10,006 ps。原输入已经使用非均匀采样；不能按恒定帧间隔推算时间。
- `run_10-12ns/GST_aging_10-12ns.dump`：8,001 原子，200 帧，step 5,005,000–6,000,000，
  每 5000 step 输出一帧，即 10 ps，存储覆盖 10,010–12,000 ps。两文件时间范围不重叠。
- liquid 目录包含 1200 K、3000 K；不能直接与当前未知/室温的非晶 GST225 参考计算同条件误差。
- `al_gen7/.../canon_liquid_3000K_seed0` 只有 29 帧；`al_gen6/.../seed1` 只有 471 帧，
  已按实际完整输出计数，没有将它们补算为常见的 601 帧。

每文件最多均匀选 400 个存储索引的数量上限草案合计 **{remote_summary['max400_stored_index_plan_frames']:,} 帧**，
索引已列在 `hopper_candidates.json`。这只是数量规划；长时间老化数据应另按物理时间选点，
并保留源 step、dt、实际时间及模型版本。尚未批量下载或执行抽样。

按用户明确要求，**忽略 Hopper `deeph` 中 4 个 216 原子结构、1 个 1080 原子结构及其处理副本**。
此前只读扫描记录保留用于审计，`deeph_static_structures.*` 中标记 `included=false`、`excluded_by_user`，
不纳入后续收集和计数。排除规则见 `inventory/collection_policy.json`，不影响新提供的本地 GST225 轨迹。

{selection_text}## 下一阶段的数据约定

1. 使用 AI_selected 的筛选结果和 metadata 实现第三个网页模块；当前尚未发布网页。
2. 保留全部原始数据、旧参考分组、原报告温度与历史误差供溯源。
3. 若继续收集 Hopper 液态高温候选，须单独处理；本轮未将远端未下载的液态结构计入筛选集合。

## 文件入口与复现

- `AIMD_reference/references.{{json,csv}}`：参考轨迹来源及物理属性。
- `AIMD_reference/rdf_summary.{{json,csv}}`、`AIMD_reference/<材料>/mean_rdf.{{json,csv}}`：初始参考均值与峰位。
- `inventory/ai_local_frames.{{json,csv}}`：每个本地 AI 结构的帧索引、方法、组分、哈希及重复项。
- `inventory/ai_local_files.{{json,csv}}`、`ai_local_summary.json`：文件哈希和汇总。
- `inventory/gst225_files.{{json,csv}}`、`gst225_summary.json`：新增 GST225 的真实尺寸、帧数、文件名差异和跨来源重复项。
- `inventory/hopper_discovery.json`：远端结构候选清单、排除的软件目录及扫描错误（本次无）。
- `inventory/hopper_counts.json`、`hopper_candidates.{{json,csv}}`、`hopper_summary.json`：远端实际帧数、文件去重及抽样草案。
- `inventory/deeph_static_structures.{{json,csv}}`：静态结构与处理副本来源。

```bash
conda activate pcm_database
python scripts/organize_ai_sources.py references --apply  # 可重复执行，校验已归档哈希
python scripts/organize_ai_sources.py rdf
python scripts/organize_ai_sources.py local
python scripts/report_ai_inventory.py
python scripts/select_ai_structures.py build --workers 2
python scripts/select_ai_structures.py verify
python -m unittest discover -s tests -v
```

Hopper 使用 SSH ControlMaster 会话复用。远端扫描使用 Python 标准库仅读文件，分析在本地 pcm_database 环境完成；
完整计数/扫描命令及控制 socket 默认值见 `scripts/inventory_hopper.py --help`。远端文件随后若继续增长，须重新清点。
"""
    (BASE / "design.md").write_text(report)
    (OUT / "README.md").write_text("# 数据清点产物\n\n完整数量、口径、限制与复现命令见 [../design.md](../design.md)。\n")
    print(f"Wrote {BASE / 'design.md'}; remote candidate frames={remote_summary['stored_complete_frames']}; exact local candidates={exact_requested}")


if __name__ == "__main__":
    main()
