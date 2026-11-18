# AI 生成结构模块：归档、多样性筛选与 metadata

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
| MLIP-MD：GeTe、Sb–Te、GST225 | 62,320 | 62,318 |
| Diffusion：Ge–Se（含纯 Ge / Se） | 12,735 | 12,735 |
| Flow Matching：Ge–Se，单独标注 | 1,350 | 1,350 |
| 合计 | 76,405 | **76,403** |

23 个 MLIP 文件各 400 帧；34 个 Ge–Se 生成文件共 14,085 帧，与原 `structures.csv`
中的组分、原子数、晶胞体积和全部几何哈希逐帧一致。GeTe 有两条重复初始帧；
新增 GST225 文件另做全量核验，逐文件数量与重复项见下表及帧级索引。
去重对元素、笛卡尔坐标和晶胞以 1e-6 Å 量化后哈希，忽略原子排列顺序；未做平移、旋转或结构相似性聚类。
MD 存储帧数量不代表统计独立样本数量。

只计本次明确提出的 **MLIP-MD + Diffusion：75,053 个去重结构**，已经达到 10,000+ 的数量目标。
其中有严格化学计量 AIMD 对应者 **74,253 个**，另有取整配比结构按下述容差匹配参考。
Flow Matching 的 1,350 个结构另外保留，不能写成 Diffusion。

800 帧暂无严格组分参考：`sb2te_1000` 为 Sb667Te333，`sbte2_1000` 为 Sb333Te667，
每个文件 400 帧。它们与理想 Sb₂Te、SbTe₂ 只有取整差异，但不能静默改组分或强配参考。
本轮按最大元素原子分数偏差 0.001（0.1 百分点）容差，匹配最近理想组分参考，
并在 metadata 中保留真实组分、参考材料、偏差和 `integer_stoichiometry_rounding` 标记。

原 Ge–Se 数据包已经按原参考组、第一峰误差 <5% 预筛选。历史误差保留在
`legacy_peak_error_pct` 中；本轮统一计算对新合并参考的误差，不把 <5% 设为新模块的准入门槛。

| 材料（实际组分） | MLIP-MD | Diffusion | Flow Matching |
|---|---:|---:|---:|
| Ge | 0 | 3,127 | 150 |
| Ge11Se9 | 0 | 134 | 0 |
| Ge2Sb2Te5 | 53,120 | 0 | 0 |
| Ge2Se | 0 | 715 | 150 |
| Ge3Se | 0 | 787 | 150 |
| Ge3Se17 | 0 | 194 | 0 |
| Ge3Se2 | 0 | 577 | 0 |
| Ge3Se7 | 0 | 912 | 150 |
| Ge4Se | 0 | 474 | 0 |
| Ge7Se3 | 0 | 955 | 150 |
| GeSe | 0 | 1,331 | 150 |
| GeSe19 | 0 | 67 | 0 |
| GeSe2 | 0 | 760 | 150 |
| GeSe3 | 0 | 808 | 150 |
| GeTe | 1,998 | 0 | 0 |
| Sb | 1,200 | 0 | 0 |
| Sb2Te | 800 | 0 | 0 |
| Sb2Te3 | 1,200 | 0 | 0 |
| Sb333Te667 | 400 | 0 | 0 |
| Sb667Te333 | 400 | 0 | 0 |
| SbTe | 1,200 | 0 | 0 |
| SbTe2 | 800 | 0 | 0 |
| Se | 0 | 1,894 | 150 |
| Te | 1,200 | 0 | 0 |

## 新增本地 GST225 多尺寸轨迹

按用户补充，将 `reorg/GST225/` 作为本地 GST225 轨迹来源。共 **5 个文件、
53,120 个完整帧、53,120 个去重几何结构**。
实际组分由 ASE 读取，全部逐帧检查有限坐标、有效周期晶胞、原子数、几何哈希与重复项。
大体积 dump 按完整帧逐块交给 ASE 解析，避免一次将整条轨迹载入内存。

| 本地原文件名 | 实际原子数 | 完整帧 | 去重后 | timestep 最小–最大 |
|---|---:|---:|---:|---|
| GST225_1080atoms.dump | 1080 | 30,001 | 30,001 | 0–3000000 |
| GST225_108atoms.dump | 108 | 1,001 | 1,001 | 0–100000 |
| GST225_216atoms.xyz | 216 | 9 | 9 | 207–627 |
| GST225_512atoms.dump | 504 | 20,001 | 20,001 | 5000–2005000 |
| GST225_720atoms.dump | 8001 | 2,108 | 2,108 | 0–5003000 |

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
其中 **1 条 / 2,108 帧已存在本地 GST225 文件中**；
排除这些完全相同的文件后，远端还有 **16,521 个额外候选帧**。

| Hopper 来源 | 去重轨迹文件 | 候选完整帧 | 每帧原子数 |
|---|---:|---:|---:|
| MLIP_benchmark：liquid | 18 | 14,316 | 306 |
| MLIP_benchmark：amorphous | 5 | 2,005 | 306 |
| GST225：373 K 长时间老化 | 2 | 2,308 | 8001 |

- `run_10ns/GST_aging_10ns.dump`：8,001 原子，2,108 帧，step 0–5,003,000，dt=0.002 ps，
  存储覆盖 0–10,006 ps。原输入已经使用非均匀采样；不能按恒定帧间隔推算时间。
- `run_10-12ns/GST_aging_10-12ns.dump`：8,001 原子，200 帧，step 5,005,000–6,000,000，
  每 5000 step 输出一帧，即 10 ps，存储覆盖 10,010–12,000 ps。两文件时间范围不重叠。
- liquid 目录包含 1200 K、3000 K；不能直接与当前未知/室温的非晶 GST225 参考计算同条件误差。
- `al_gen7/.../canon_liquid_3000K_seed0` 只有 29 帧；`al_gen6/.../seed1` 只有 471 帧，
  已按实际完整输出计数，没有将它们补算为常见的 601 帧。

每文件最多均匀选 400 个存储索引的数量上限草案合计 **9,429 帧**，
索引已列在 `hopper_candidates.json`。这只是数量规划；长时间老化数据应另按物理时间选点，
并保留源 step、dt、实际时间及模型版本。尚未批量下载或执行抽样。

按用户明确要求，**忽略 Hopper `deeph` 中 4 个 216 原子结构、1 个 1080 原子结构及其处理副本**。
此前只读扫描记录保留用于审计，`deeph_static_structures.*` 中标记 `included=false`、`excluded_by_user`，
不纳入后续收集和计数。排除规则见 `inventory/collection_policy.json`，不影响新提供的本地 GST225 轨迹。

## 多样性筛选与正式 metadata

从 76,403 个几何去重候选选择 **14,900 个结构**，上限 14,900。Flow Matching 已纳入。

| 建模方式 | 筛选后结构 | 峰位误差中位数 / % | 峰位误差 P95 / % |
|---|---:|---:|---:|
| Diffusion | 4,404 | 0.4762 | 1.9345 |
| Flow Matching | 1,233 | 0.5235 | 1.4746 |
| MLIP_DeePMD | 833 | 0.5385 | 1.0498 |
| MLIP_MACE | 4,416 | 0.8992 | 1.6827 |
| MLIP_NEP | 4,014 | 0.4711 | 1.2867 |

按参考材料、建模方式与实际原子数分层，用候选数平方根分配预算并保留稀有尺寸。
用元素分辨的径向分布、配位数分布及密度表示径向和配位多样性，缩放、PCA 后最远点采样。
已知时间的老化轨迹按时间数量级分层，每层按物理时间均匀选点。此轮不以峰位误差大小筛选。

AIMD 和对应评价温度按用户约定为 300 K，原报告温度另行保留。主峰误差逐结构重新计算。
0 个候选因参考或内部主峰无法定位而未进入可评价集合；
已选结构峰位状态为 {"clear": 14900}。
模拟时间数量级分布为 {"none": 10371, "0ps": 1, "0.01ps": 44, "0.1ps": 428, "1ps": 3090, "10ps": 695, "100ps": 121, "1ns": 120, "10ns": 16, "0.001ps": 14}；缺失者为 none。

产物：`AI_selected/metadata.json` / `metadata.csv`、`structures/`、`selection_manifest.json`、
`selection_strata.csv` 和 `verification.json`（完整验证后生成）。结构帧保留原始字节，来源哈希可追溯。
每个 metadata 记录建模方式、模拟时间/数量级、300 K 参考、主峰位置及百分比误差。

## 下一阶段的数据约定

1. 使用 AI_selected 的筛选结果和 metadata 实现第三个网页模块；当前尚未发布网页。
2. 保留全部原始数据、旧参考分组、原报告温度与历史误差供溯源。
3. 若继续收集 Hopper 液态高温候选，须单独处理；本轮未将远端未下载的液态结构计入筛选集合。

## 文件入口与复现

- `AIMD_reference/references.{json,csv}`：参考轨迹来源及物理属性。
- `AIMD_reference/rdf_summary.{json,csv}`、`AIMD_reference/<材料>/mean_rdf.{json,csv}`：初始参考均值与峰位。
- `inventory/ai_local_frames.{json,csv}`：每个本地 AI 结构的帧索引、方法、组分、哈希及重复项。
- `inventory/ai_local_files.{json,csv}`、`ai_local_summary.json`：文件哈希和汇总。
- `inventory/gst225_files.{json,csv}`、`gst225_summary.json`：新增 GST225 的真实尺寸、帧数、文件名差异和跨来源重复项。
- `inventory/hopper_discovery.json`：远端结构候选清单、排除的软件目录及扫描错误（本次无）。
- `inventory/hopper_counts.json`、`hopper_candidates.{json,csv}`、`hopper_summary.json`：远端实际帧数、文件去重及抽样草案。
- `inventory/deeph_static_structures.{json,csv}`：静态结构与处理副本来源。

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
