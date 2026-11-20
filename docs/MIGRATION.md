# 网络代码迁移

| 原代码 | 当前实现 |
|---|---|
| Judgment_reaction_pair_direct_des_2.py | resources/direct_rules.json；有序状态和反应规则 |
| Judgment_reaction_pair_out_des.py | resources/out_rules.json；恢复字典键覆盖丢失的 2 条规则 |
| generate_graph 系列重复组图 | network.build_network；规则与能量输入配置化 |
| DijkstraAlgo.py | analysis.minimax_path、maximin_path、best_path |
| 路径频次脚本 | analysis.path_statistics；DAG 精确计数和有环图有界枚举 |
| 各版本筛边、子网 | analysis.prune_network |
| Markov.py 与稳态脚本 | 可选 kinetics 模块；使用活化能垒构建有权状态转移 |
| 反应自由能图与网络图 | plotting；JSON/CSV/DOT 输出及按需绘图 |
| C_valence_state_of_ads.py | chemistry.composition、carbon_valence |

模型训练、结构优化、过渡态、计算程序适配和通用作业串联已移出公开程序包。

用户提供的历史 FAIR-Chem 修改版保留在 src/co2rr_network/_vendor/fairchem，内部模块引用改为私有前缀，并附原主许可证和子模块许可证。网络模块不引用该实现。原工作区中的独立研究脚本保持原样。
