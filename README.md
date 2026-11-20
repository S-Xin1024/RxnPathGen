# CO2RR Network

本库负责反应路径网络的构建、分析和瓶颈路径搜索。主接口接收一个能量字典：键为物种或反应式字符串，值为对应能量；同时提供 ΔG 和能垒时使用小字典。

```python
from co2rr_network import build_network, minimax_path

energies = {
    "CO+CO -> COCO": {"delta_g": 0.4, "barrier": 0.6},
    "COCO -> COCHO": {"delta_g": -0.3, "barrier": 0.2},
}
graph = build_network(energies)
path = minimax_path(graph, "CO+CO", "COCHO", weights="barrier")
```

也支持 `{"CO+CO -> COCO": 0.4}` 表示反应 ΔG，或配合 `energy_type="barrier"` 表示能垒。只给状态自由能时，传入 `{"CO": -0.1, "COCO": 0.2, ...}` 并使用反应规则建边。

目标物种只有唯一前驱时，可以直接用物种名携带 delta_g/barrier；存在多个前驱时用反应式区分。未提供的绝对状态能量或 ΔG 保留为 None。

核心功能使用 Python 标准库，无必装的第三方运行依赖：

- 由反应规则与能量表组装有向网络，保留有序位点、共吸附状态和反应通道。
- minimax：在全部路径中最小化最大边权。
- maximin：在全部路径中最大化最小边权。
- 热力学与能垒路径筛选、路径重建及瓶颈步骤定位。
- 可达性、有效子网络剪枝、DAG 精确路径计数和有环图的有界简单路径统计。
- JSON、CSV 和 Graphviz DOT 导出。

Markov/稳态分析和绘图为按需调用的分析模块，分别使用 NumPy/SciPy 和 Matplotlib。默认导入和网络构建、路径搜索不加载这些库。

## 使用

```powershell
python -m pip install . --no-deps
co2rr network examples/minimal/network.json
```

也可直接从源码运行：

```powershell
$env:PYTHONPATH = "src"
python -B -m co2rr_network network examples/minimal/network.json
```

结果保存到配置指定的目录：network.json、states.csv、reactions.csv、report.json、network.dot。命令行入口仅包含：

```text
co2rr network CONFIG
co2rr rules --profile direct
co2rr path NETWORK SOURCE TARGET
co2rr paths NETWORK SOURCE TARGET
co2rr prune NETWORK SOURCE TARGET --output FILE
co2rr import-legacy INPUT OUTPUT --facet 111
co2rr plot INPUT OUTPUT
```

## 反应识别与 COO 拓扑

```python
from co2rr_network import recognize_reaction_pairs, make_reaction_graph

pairs = recognize_reaction_pairs(
    comp_src={"CO", "COCO"},
    comp_dst={"COCO", "COCHO"},
)
# {"CO+CO -> COCO", "COCO -> COCHO"}

(indices,) = make_reaction_graph(pairs)
# [[0, 2], [2, 1]]，逻辑形状为 (2, E)

indices, node_labels = make_reaction_graph(pairs, return_labels=True)
# node_labels == ("CO+CO", "COCHO", "COCO")
```

comp_src/comp_dst 是单组分名称的候选池，接受 list 或 set，不表示某一次反应的化学计量表。因此 {"CO"} 可以按偶联规则识别 CO+CO，输出保留其计量数。只识别给定规则支持的反应；未知组分不会仅因元素配平而被判定为反应。可通过 rules 扩展规则。

COO 的第一行是源索引，第二行是目标索引。节点代表整条反应式一侧的完整状态，CO+CO 是一个节点。节点标签和边按字典序排列，确保 set 的迭代顺序不影响编号。空输入返回 ([[], []],)，不自动添加逆边。

该接口使用原生 Python 列表，不依赖图框架。现有 ReactionNetwork 是可选的带属性 dataclass 容器，能量、反应类型及元数据保存在其中；它不是 DGL、PyG 或 NetworkX 对象。

## 路径搜索接口

```python
from co2rr_network import ReactionNetwork, minimax_path, maximin_path, best_path

graph = ReactionNetwork.load("network.json")

# 原始带符号边权：不截断负数，不加统一平移量。
path = minimax_path(graph, "CO+CO", "CHOCHO_d", weights="delta_g")

# 通用最大瓶颈容量问题；weights 也可以是反应 ID 到权重的字典或函数。
capacity_path = maximin_path(graph, "CO+CO", "CHOCHO_d", weights=my_edge_weights)

# 保留化学筛选约定：代价 max(0, ΔG)，或使用已提供的活化能垒。
thermo = best_path(graph, "CO+CO", "CHOCHO_d", mode="thermodynamic")
kinetic = best_path(graph, "CO+CO", "CHOCHO_d", mode="kinetic")
```

minimax 与 maximin 是不同的图优化问题。大边权表示高反应能垒时应使用 minimax；maximin 适用于大权重代表高容量或高可靠性的场景。两种算法均返回路径、反应 ID、逐步权重和瓶颈边。无边的空路径不定义瓶颈权重，因此原始权重接口返回 cost=null。

```powershell
co2rr path network.json "CO+CO" CHOCHO_d --weight delta_g --objective minimax
co2rr paths network.json "CO+CO" CHOCHO_d
```

## 数据格式

Python 主接口使用 dict；保留命名列 CSV 的导入边界，供已有数据迁移：

```csv
catalyst,facet,species,free_energy,structure_id,converged
Cu,111,CO,-0.10,co_1,true
Cu,111,COCO,0.20,coco_1,true
```

`CO+CHO` 表示两个有序位点的完整状态，`COCHO` 表示一个 C2 中间体，`_d` 为脱附终态标签。显式共吸附能量优先；缺失时可使用单吸附形成自由能之和，近似方式会记录在节点中。不同催化剂和晶面不会混合。已存在的异构体能量不自动合并。

自由能定义、CHE 电位约定、能垒表及自定义规则见 [数据接口](docs/DATA.md) 和 [算法说明](docs/METHODS.md)。原始总能量到形成自由能的计算由调用方完成。

## CO2RR 规则

- direct：629 个状态，1102 条反应，包含脱附。
- out：502 个状态，957 条反应，不包含脱附。
- 自定义规则使用 source/target 数组和 kind，不依赖固定节点编号。

完整规则示例使用明确标为 synthetic/DEMO 的合成能量，仅演示网络构建：

```powershell
python examples/full_catalog/generate_demo.py
co2rr network examples/full_catalog/network.json
```

## 本地 FAIR-Chem 源码

用户提供的历史修改版已内联到 `src/co2rr_network/_vendor/fairchem/`。数值实现保留，只调整内部命名空间引用和动态模块发现。原 MIT 许可证位于 [fairchem/LICENSE](src/co2rr_network/_vendor/fairchem/LICENSE)，子模块许可证一并保留。

网络核心不依赖 FAIR-Chem。该源码不通过外部 fairchem 包替代，不自动安装其计算依赖；包内也不提供模型训练、结构优化、过渡态或计算作业调度入口。

## 开发

```powershell
$env:PYTHONPATH = "src"
python -B -m unittest discover -s tests -v
```

[迁移说明](docs/MIGRATION.md) · [来源与许可](NOTICE.md)
