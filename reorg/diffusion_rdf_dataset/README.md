# 生成结构与DFT参考：第一峰位配对数据集

本地收集于2026-11-16。生成结构、DFT参考分别存放；按真实元素数的约分比例严格对应，允许原子数不同。

共检查 **14,335** 个去重候选，收集 **14,085** 个第一峰位误差严格小于5%的结构。

其中扩散生成 **12,735** 帧，Flow Matching **1,350** 帧；可按structures.csv中的generation_method列区分。纯Ge 3,277帧、纯Se 2,044帧，其余Ge–Se合金 8,764帧。

| 原子数 | 候选 | 收集 |
|---:|---:|---:|
| 100 | 100 | 100 |
| 300 | 14,213 | 13,963 |
| 1,000 | 21 | 21 |
| 30,000 | 1 | 1 |

| 组分 | 候选 | 收集 |
|---|---:|---:|
| Se | 2,102 | 2,044 |
| GeSe19 | 88 | 67 |
| Ge3Se17 | 214 | 194 |
| GeSe3 | 978 | 958 |
| Ge3Se7 | 1,082 | 1,062 |
| GeSe2 | 910 | 910 |
| GeSe | 1,501 | 1,481 |
| Ge11Se9 | 153 | 134 |
| Ge3Se2 | 588 | 577 |
| Ge2Se | 865 | 865 |
| Ge7Se3 | 1,116 | 1,105 |
| Ge3Se | 949 | 937 |
| Ge4Se | 487 | 474 |
| Ge | 3,302 | 3,277 |

## 查看与使用

- `generated/<组分>/atoms_<原子数>.extxyz`：通过验收的生成结构，按组分与尺寸合并；保留原始每帧字节。
- `dft_reference/<组分>/XDATCAR_<ID>`：独立的原始DFT参考轨迹。
- `structures.csv`：每行对应一个已收集结构；`output_frame`和源文件帧索引均从0开始，数据库来源的`locator`是attempt记录ID。
- `references.csv`：参考轨迹、实际组分、已知温度、使用帧号、独立峰位和参考组峰位。`included=False`的轨迹未用于基准、未复制。
- `collection_manifest.json`：机器读取的固定候选范围、源文件哈希、排除原因、参数和结果；一般无需人工浏览。

## 验收定义

总数目g(r)，0.05 Å bin，按晶胞体积和有限N原子对数归一；第一配位壳2.0–3.3 Å，统一5点滑动平均，主峰三点插值。计算邻近原子对，未创建大胞的N×N矩阵。边界峰、缺峰和明显多峰暂不收集。

`error_pct = 100 * abs(r_generated - r_DFT) / r_DFT < 5`。该指标是第一峰位误差，不代表整条g(r)曲线误差或全部结构性质的准确性。未按SOAP等额外压缩多样性。

参考用每条轨迹最后1000步内均匀抽出的200帧，完整2000帧原轨迹全部保留。相同组分的新300 K多seed轨迹按轨迹等权平均g(r)建立一个参考组；旧轨迹缺少充分条件记录时各自独立成组，不与新300 K轨迹混合。生成结构必须对其组分的全部参考组均通过5%；索引的`max_error_pct`保存最差误差，`error_reference_group_id`与`reference_peak_A`指产生该最差误差的组。参考组在筛选前确定，不逐帧挑更有利的参考。

参考稳定性使用后半段前100/后100个抽样帧的均值峰位比较，漂移≥5%或无法明确定位者不作为基准。旧轨迹温度未查明时留空并标记legacy，不自动标作300 K。

实际包含DFT参考44条；下列参考暂未纳入：
- `workspace_docs/supercomputer/data_increasement/Ge15Se285/seed_01/4-md/XDATCAR`：normal_end_marker_not_found_pending_review。

## 来源与数量口径

候选来自本地新旧生成模型、归档版本、消融与Flow Matching等分支，来源在索引中保留，不将这些结果全部当作当前主模型的表现。训练POSCAR、MLP弛豫轨迹、DFT优化中间步、骨架截取和未训练冒烟测试未混作独立生成样本。

原工作区逐帧初步去重候选93,733个，其中14,335个与DFT严格对应；无对应者不进入本包。纯Ge/Se也计入总量；参考轨迹帧数不计为生成数量。

筛选结果：
- `accepted`：14,085。
- `ambiguous_first_shell_peak`：57。
- `boundary_or_absent_peak`：47。
- `peak_error_ge5pct`：146。

## 验证与复现

依赖Python≥3.10、NumPy、SciPy。使用现有chaldiff环境即可：

```bash
/home/ysj/anaconda3/envs/chaldiff/bin/python -B collect_dataset.py verify
# 有原始工作区时，按固定清单在另一个新目录完整重建：
/home/ysj/anaconda3/envs/chaldiff/bin/python -B collect_dataset.py build \
  --source-root /home/ysj/amor_diff/amor_gen/graphite-main \
  --output /tmp/diffusion_rdf_dataset_rebuild --workers 4
```

verify无需原始工作区，逐帧重算峰位、核对严格组分对应和字节/几何哈希。build核对原始源文件哈希，使用固定候选范围和参考规则重建；拒绝覆盖已有生成数据。完整参考基准重算可使用verify --references。

收集后已核对14,085帧与索引一致，所有原始帧字节、参考文件哈希以及逐帧峰位复算通过。原工作区文件和public_release未改动，未上传。
