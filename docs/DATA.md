# 网络数据接口

主接口为 `build_network(energies: dict, rules=None, ...)`，不要求文件路径、催化剂前缀或计算后端。JSON/CSV 只作为命令行的数据输入输出格式；相对路径以配置文件所在目录为基准。

## 统一能量字典

| 键 | 值 | 含义 |
|---|---|---|
| COCO | 0.2 或 {"free_energy":0.2} | 该物种/状态的形成自由能 |
| CO+CO -> COCO | 0.4 | 反应自由能变化 ΔG |
| CO+CO -> COCO | {"delta_g":0.4,"barrier":0.6} | 同时提供反应 ΔG 和活化能垒 |
| CO+CO -> COCO | {"barrier":0.6} | 只提供能垒，ΔG 保持未知 |

`energy_type="barrier"` 或 `"delta_g"` 可指定纯标量字典的含义。物种名携带逐反应能量时，只有规则中存在一个唯一前驱才接受，否则使用反应式键。只提供反应式字典时，拓扑直接从键构建；提供状态能量时使用给定规则或默认 CO2RR 规则。

仅提供边能量时，节点 free_energy 为 None；只给能垒时边 delta_g 为 None。不猜测绝对能量或将缺失数据填为零。状态能量和反应 ΔG 同时提供时必须一致。反向能垒的推导需要已知 ΔG。

反应 ΔG 与状态形成自由能以 U=0、pH=0 为输入约定；若传入非零 CHE 条件则一并修正。能垒作为所请求条件下的数值使用，不自动外推。

## 自由能

必需列：catalyst、species、free_energy。可选列：facet、structure_id、converged。

free_energy 是 U=0、pH=0 时相对于 CO2/H2/H2O 的形成自由能，单位 eV，已经完成裸表面参考和相应热力学修正。总能量、吸附能和形成自由能不可混填。

CSV 导入时，同一 catalyst/facet/species 的多个收敛构型取最低自由能；能量相同时按 structure_id 稳定选择。converged 支持 true/false、1/0、yes/no。缺失反应能量默认跳过，missing_energy=error 要求候选规则全部具有可用数据。缺失状态能量与跳过的反应记录在元数据中。

CO+CHO 表示有序双位点状态，具有该 ID 的显式能量行优先。缺失时将单吸附形成自由能相加，并标记 independent_adsorbates。单吸附形成自由能已各自减去裸表面，不能直接把两个包含完整表面的总能量相加。

aliases 只用于缺失能量的显式 fallback，已有原名能量优先；该配置不改变节点标签或断言两个构型相同。

## 能垒

字典中使用 barrier 字段，或使用 energy_type="barrier" 的标量值。已有 CSV 的导入边界使用 catalyst、source、target、barrier 列，可附 facet、potential、potential_scale、ph、temperature。

source、target 必须是完整状态 ID。barrier 是相对于 source 的正向活化自由能，不是 TS 总能量，也不是反应自由能。显式给出计算条件时必须与网络配置一致，程序不自动把能垒外推到其他电位。

输入由调用方统一为所选择的能量约定，程序不处理过渡态计算或能垒热力学修正。已知 ΔG 时要求 Ea≥max(0,ΔG)。同一反应式不能重复填入多个值；同一端点对的平行规则使用该反应式提供的能量。

## 自定义规则

```json
{
  "formulas": {"site1_CO": "CO", "site2_CO": "CO", "OCCO": "COCO"},
  "reactions": [
    {"source": ["site1_CO", "site2_CO"], "target": ["OCCO"], "kind": "cc_coupling"}
  ]
}
```

formulas 为带位点标签的物种提供 CHO 组成，仅用于元素守恒与 CHE 配平，不编码键连、电荷或机理。加号保留给状态组合。规则库是可扩展的有限候选网络。

## 网络对象

`recognize_reaction_pairs(comp_src, comp_dst)` 返回规则允许的 `Set[str]`，例如 {"CO+CO -> COCO"}。输入是候选单组分池，列表重复项与集合具有相同含义；反应中的重复组分和有序位点由规则保留。

`make_reaction_graph(pairs)` 默认返回 `(indices,)`。indices 为两行等长整数列表，逻辑形状是 (2,E)，每列给出一条有向边。可选择 return_labels=True 返回 `(indices,node_labels)`，其中 node_labels[i] 对应节点 i。节点与边均按字典序排列，归一化后的重复边合并，不添加权重或反向边。

COO 拓扑可独立使用，不要求创建 ReactionNetwork。ReactionNetwork 是使用标准库 dataclass 定义的带属性图容器；内部以字符串标记端点，路径算法按需要构建邻接表。

ReactionNetwork 包含 State 列表与 Reaction 列表。State 包含 ID、物种元组、可选形成自由能、终态标志和元数据；Reaction 包含唯一 ID、源、目标、类型、可选 ΔG、可选活化能垒和元数据。平行反应保留独立 ID。所有节点 ID 必须唯一，所有反应端点必须存在；绝对状态能量与 ΔG 均已知时必须一致。

对通用加权图算法，使用 minimax_path/maximin_path 的 weights 参数传入自定义映射或函数，无需改变节点自由能或把容量硬塞进物理能垒字段。
