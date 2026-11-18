# AI 生成结构：多样性筛选结果

从 76,403 个几何去重候选中选择 **14,900** 个结构（上限 14,900）。

`metadata.json` / `metadata.csv` 保存建模方式、源文件/帧号/哈希、300 K 对应参考、第一主峰及百分比误差。
建模方式为 MLIP_MACE、MLIP_DeePMD、MLIP_NEP、Diffusion、Flow Matching。
AIMD 及 MLIP 对应温度按用户约定统一为 300 K；原报告温度另保留，不表示重新测得温度。
未知模拟时间在 JSON 的数值字段为 null，CSV 为 none；simulation_time_order 两种格式均为 none。
已知时间以真实 ps 值记录，并按 10 的幂标注数量级，例如 1ps、10ps、100ps、1ns、10ns；0 时刻单列 0ps。

Sb–Te XYZ 的 Time 字段按匹配 GPUMD 官方格式推断为 fs，除以 1000 转 ps，依据保存在 simulation_time_basis。
单位来源：[GPUMD 官方输出代码](https://github.com/brucefan1983/GPUMD/blob/master/src/measure/dump_observer.cu)。

按参考材料、建模方式和真实原子数分层，用候选数平方根分配名额并保留小类。
511 维描述符包含元素分辨的 1.5–5.5 Å 径向分布、3.5 Å 配位数分布及数密度；
分层缩放后 PCA 自适应保留至少 95% 方差，再做确定性的最远点采样。描述符侧重径向和配位多样性，不保证识别所有角度/拓扑差异。
已知时间的老化轨迹先按时间数量级分层，每层按物理时间均匀选点，避免按非均匀存储帧索引选点。
strata CSV 和 manifest 保存候选数、配额、实际采样数、PCA 方差覆盖及最远覆盖距离。

峰位统一使用总数目 g(r)，0.05 Å 分箱、有限 N 原子对归一、5 点滑动平均，
2.0–3.3 Å 第一壳主峰三点插值；误差为 100*abs(r_AI-r_AIMD)/r_AIMD。
多个局部峰时记录占优主峰及歧义标记；未能定位内部主峰者不进入本次可评估集合。未按误差大小设筛选门槛。
Sb667Te333 / Sb333Te667 按最近理想组分对应参考，最大原子分数偏差容差为 0.001（0.1 百分点），
真实组分、参考组分和实际偏差全部保留。

structures/ 内保留原始帧字节；LAMMPS dump 须使用 metadata 的 type map 用 ASE 读取。
没有删除或裁剪原始轨迹，也没有发布网站。

```bash
conda activate pcm_database
python scripts/select_ai_structures.py build --workers 2
python scripts/select_ai_structures.py verify
```
