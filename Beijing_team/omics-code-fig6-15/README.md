# Fig6 / FigS2 的 15 张图：源码精简版

本地 Git 仓库仅整理目标压缩包 `output/pdf/fig6_48pt_editable/volcano_balanced_v2/Fig6_JM_FigS2_15_vector_PDFs_Arial48_balanced_v2.zip` 对应的制图代码，以及为这些图提供差异分析、GSEA 和 GO/KEGG 结果的 bulk 前置代码。**不复制 PDF、测序数据、结果矩阵或其他图。**

## 15 张图与代码

| 图 | 数量 | 制图源码 |
| --- | ---: | --- |
| Fig6 J、FigS2 A 火山图 | 2 | `figures/volcano_balanced_v2/draw_volcano_balanced_v2.py` |
| Fig6 K、FigS2 B 的 GO BP / KEGG 单图 | 8 | `figures/draw_go_kegg_48pt.py` |
| Fig6 L/M、FigS2 C/D/E 的 GSEA 单图 | 5 | `figures/draw_gsea_panels.py` |
| 按目标文件名组合 ZIP | 15 | `figures/volcano_balanced_v2/package_balanced_v2.py` |

目标 ZIP 中两张火山图的成员名保留旧名称，但 PDF 内容来自 `balanced_v2` 新版。`TARGET_PANELS.csv` 记录 15 个 ZIP 成员的 SHA-256；`SOURCE_MANIFEST.csv` 记录收录源码的原始路径和 SHA-256。

## 前置 bulk 代码

| 步骤 | 源码目录 | 输出给后续的内容 |
| --- | --- | --- |
| PDC/PDO 标签修正、QC、配对差异、通路分析 | `bulk/pdc_pdo/`，按 01–04 顺序 | PDC/PDO 火山数据、DE、GSEA、GO/KEGG |
| Ulixertinib/Ctrl 样本表、featureCounts、配对差异、通路分析 | `bulk/ulixertinib/`，05b → featureCounts → 05c → 05f | 药物对照 DE、GSEA、GO/KEGG |
| 两组共用的全库 GSEA、RDS 转换、选定通路清单与 GO/KEGG top 10 | `bulk/shared/` | GSEA manifest 和 8 张富集图的 top 10 表 |
| 从已保存的 bulk 结果提取这 15 张图使用的 5 条 GSEA 曲线和药物火山 CSV | `figures/prepare_plot_inputs.py` | 6 个绘图输入 CSV |

原始 `05d_Ulixertinib_pathway_GSEA_GO_KEGG.py` 的运行日志记录了输入缺失异常；已完成的结果对应 `05f_Ulixertinib_GSEA_GO_KEGG.py`，本版只收入 05f。此前的大仓库中的 D/Q 配对分析、去卷积、PPT、美化试验稿及旧火山图脚本均不在此版。

## 输入与使用

四个制图脚本读取的 10 个直接输入文件列在 `REQUIRED_INPUTS.csv`，当前均位于上级 PCP 工作区；它们没有复制进仓库。`figures/prepare_plot_inputs.py` 可从保存的 bulk 结果重建其中 6 个 CSV，其余 4 个由上表的 bulk 脚本生成。该脚本还需原工作区已有的 `output/pdf/fig6_v3/data/gsea_exact_inputs.json`；未找到这个 JSON 的生成源码。

图脚本原样保留了 `D:\projects\PCP` 路径，bulk 原脚本还保留了 `/home/susu/...` 或 `E:\PCP` 路径。若在另一台机器运行，先修改路径并准备相应输入、Arial 字体和依赖。`05b_Ulixertinib_audit.py` 会先写样本表，随后读取原位置缺失的 `vendor_gene_count.xls` 做额外核对；运行前需提供该文件或调整其路径。当前工作区也未找到更早的 fastp/STAR 脚本及原始 BAM/GTF，因此 bulk 的原始测序起点无法仅凭本仓库复跑。

在本仓库目录检查源码完整性：

```powershell
python tools/verify_sources.py
```

如需在原 PCP 工作区更新源码清单：

```powershell
python tools/collect_sources.py
```

重收集会覆盖仓库中的原样源码；`figures/prepare_plot_inputs.py` 是从原 `fig6_v3/rebuild_figures.py` 提取的聚焦版本，不受该操作覆盖。Python 包名见 `requirements-python.txt`；R 脚本还需 `clusterProfiler`、`fgsea`、`msigdbr`、`org.Hs.eg.db`、`enrichplot`、`ggplot2`、`dplyr` 和 `tidyr`。
