# 剑鱼论文的逐步分析

[English](README.md)

这 32 个步骤重复论文（`docs/swordfish_bycatch_manuscript.md`）及其补充材料（`docs/supplementary_tables.md`）的
分析，从 IATTC 公开延绳钓文件开始，到每一幅图、每一张表和每一项模拟，文末列出的四项除外。各步骤按方法部分的顺序
排列。每一步是一个从上到下阅读和运行的短脚本，每条注释先写英文，再写中文。

四个文件存放多个步骤共用的部分：`common.R` 把新表和新图与存档结果比较并进行自助抽样，`model_functions.R` 拟合
各船队模型，`simulation_functions.R` 设定模拟，`plot_functions.R` 存放绘图设定。`run_all.R` 按顺序运行 32 个步骤。
所有步骤都是 R 脚本。

## 需要安装的软件

- R 4.6.0，以及 sdmTMB 1.2.0、TMB 1.9.25、fmesher 0.8.0、ggplot2 4.0.4、ggsci 5.2.0、sf 1.1.1、patchwork 1.3.2、
  rnaturalearth 1.2.0、rnaturalearthdata 1.0.0、jsonlite 2.0.0、cluster 2.1.8.2、generics 0.1.4、scales 1.4.0、
  ragg 1.5.2 和 png 0.1.9。R 会提示 sdmTMB 是在 R 4.6.1 下编译的，这不影响结果。

## 设置

把这些脚本放在项目文件夹的 `analysis_walkthrough/` 文件夹中。各步骤从 `data/raw/` 读取 IATTC 公开延绳钓文件
（https://www.iattc.org/getmedia/b8f0bdbb-595d-4c16-9965-cbedcf122aaa/PublicLLTunaBillfish.zip），从 `config/` 读取设置，核对时读取 `results/` 和 `tmp/` 中的存档结果。

在项目文件夹（包含 `data/`、`config/` 和 `results/` 的文件夹）中打开终端，在运行 R 之前先设定下面这些变量，使每个
R 进程只用一个核心，并使 R 能写入中文文件夹名。`run_all.R` 会自行设定这些变量。

```sh
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
```

## 步骤

按顺序运行各步骤。每一步只写入 `analysis_walkthrough/output/`；`data/`、`results/`、`tmp/` 和 `config/` 中的文件只
读取。所有步骤都读取 `config/directed_fishery.json`。

### 渔业、观测与预测区域（第 2.1 节）

1. `01_check_public_data.R` 读取保存的 IATTC 压缩包和 `PublicLLTunaBillfishNum.csv`，检查两者是同一文件。不写出文件。
2. `02_draw_study_design.R` 不读取文件，写出 `Fig1_study_design`（图 1）。
3. `03_prepare_records.R` 读取公开文件，把六支船队在 H2 内的记录及各网格的位置和面积写入 `records.rds`。
4. `04_define_areas.R` 读取 `records.rds`，写出各网格有努力量的年数（`support.csv`）、有剑鱼正尾数的网格（`footprint.csv`）、31 个拟合（`tasks.csv`）、各评估分区的网格（`area_structure.csv`，表 S1）和模型输入 `input.rds`。
5. `05_map_catch_cells.R` 读取 `footprint.csv`，写出 `Fig2_coverage`（图 2）和 `coastal_grid_check.csv`。
6. `06_count_catch_cells.R` 读取公开文件，写出各时段有剑鱼的网格（`period_coverage.csv`），以及金枪鱼延绳钓船队新增网格中有正尾数的年份（`added_cell_positive_years.csv`）（第 3.1 节、表 S1 注释）。

### 与西班牙重叠区域内的年度指数（第 2.2 节）

7. `07_fit_baseline_models.R` 读取 `input.rds`，以全部合格记录拟合各船队的基线模型，写出 `fits_baseline.csv`、`indices_baseline.csv`，并在 `fits/` 中写出各拟合的诊断、按网格的残差和各网格的相对 CPUE。
8. `08_spain_area_indices.R` 读取 `input.rds` 和 `indices_baseline.csv`，写出各评估分区的西班牙指数（`spain_area_indices.csv`，图 S1）。
9. `09_compare_annual_indices.R` 读取 `indices_baseline.csv`，写出与西班牙的秩相关（`pair_comparisons.csv`，表 S1）和年度指数的精度（`index_precision.csv`，表 S3）。
10. `10_count_annual_overlap.R` 读取 `input.rds` 和 `indices_baseline.csv`，写出逐年两支船队都有记录的网格（`pair_overlap.csv`）以及在这些年份上中国与西班牙的相关（`china_spain_overlap_correlations.csv`）（补充材料 S1）。
11. `11_plot_annual_indices.R` 读取 `indices_baseline.csv`、`support.csv` 和 `spain_area_indices.csv`，写出 `Fig3_annual`（图 3）和 `FigS1_spain_areas`（图 S1）。
12. `12_check_residuals.R` 读取第 07 步在 `fits/` 中的文件，输出表 S2，并写出 `spatial_residual_cell_years.csv`、`residual_pair_correlations.csv`（补充材料 S4）、`FigS2_fleet_residuals`（图 S2）和 `FigS3_spatial_residuals`（图 S3）。
13. `13_retrospective.R` 读取 `input.rds` 和 `indices_baseline.csv`，写出 `retrospective_fits.csv`、`retrospective_indices.csv`、`retrospective_rho.csv`（表 S4）和 `FigS4_retrospective`（图 S4）。
14. `14_quarterly_change.R` 读取 `input.rds` 和 `fits_baseline.csv`，写出相同网格和季度内中国与西班牙的变化比（`quarterly_matched_cells.csv`、`quarterly_fit_checks.csv`、`quarterly_shared_change.csv`、`quarterly_between_seasons.csv`；图 4b，补充材料 S9）。
15. `15_compare_cells_fished_by_both.R` 读取第 07 步的相对 CPUE、`support.csv`、`pair_overlap.csv` 和 `quarterly_shared_change.csv`，写出 `shared_cell_cpue_by_year.csv`、`shared_cell_direction_by_year.csv`（表 S11）、`shared_cell_cpc_spain_ratios.csv`、`Fig4_local_contrasts`（图 4）、`FigS11_shared_cell_cpue`（图 S11）、`FigS12_local_contrasts_2021_2023`（图 S12）和 `shared_cell_cpue_atlas.pdf`（补充地图）。

### 年度指数对渔获物种组成的敏感性（第 2.3 节）

16. `16_describe_catch_composition.R` 读取 `input.rds`，写出主成分与聚类的选择（`composition_selection.csv`，表 S5）、`china_cluster_swordfish.csv`、载荷、聚类组成和各聚类的年度占比（排除和包含剑鱼），以及 `FigS5_fleet_pca`、`FigS6_fleet_clusters`、`FigS7_fleet_pca_inclusive` 和 `FigS8_fleet_clusters_inclusive`（图 S5 至 S8）。
17. `17_fit_composition_models.R` 读取 `input.rds`、`fits_baseline.csv` 和 `composition_selection.csv`，拟合其余 25 个模型，并为全部 31 个拟合写出 `fits.csv`、`indices.csv`、`model_diagnostics.csv`（表 S2 和 S6）、`model_parameters.csv`（表 S2）和 `residual_quantiles.csv`。
18. `18_summarize_composition.R` 读取 `indices.csv` 和 `fits.csv`，写出 `composition_adjustment_summary.csv`（表 3）和 `composition_changes.csv`。
19. `19_plot_composition_indices.R` 读取 `indices.csv`，写出 `Fig5_composition`（图 5）和 `FigS9_composition_inclusive`（图 S9）。

### 空间覆盖与模拟的区域丰度（第 2.4 节）

20. `20_simulation_inputs.R` 读取 `input.rds` 和公开文件，写出取样设计（`coverage_sampling_designs.csv`，表 2）和模拟的渔获水平（`simulated_catch_levels.csv`，表 S7）。
21. `21_simulate_coverage.R` 读取 `input.rds` 和公开文件，模拟情景 A1 至 A3 在三种设计下的渔获（900 个拟合），写出 `coverage_simulation.csv`、`coverage_annual_intervals.csv` 和 `coverage_uncorrected_check.csv`。
22. `22_summarize_coverage.R` 读取第 21 步的表，写出 `coverage_summary.csv`（表 S9）、`coverage_annual_paired_summary.csv`（图 6a）和 `coverage_paired_errors.csv`。
23. `23_coverage_change_intervals.R` 读取 `coverage_simulation.csv` 和 `coverage_paired_errors.csv`，写出变化中位数及其自助法区间（`coverage_change_intervals.csv`，图 6b）。
24. `24_simulate_equal_records.R` 读取 `input.rds` 和公开文件，模拟记录数相同的两种设计（600 个拟合），写出 `equal_records_samples.csv`、`equal_records_fits.csv` 和 `equal_records_annual.csv`。
25. `25_summarize_equal_records.R` 读取第 21 和 24 步的表，写出 `equal_records_balance.csv`、`equal_records_pairs.csv` 和 `equal_records_summary.csv`（表 S8）。
26. `26_simulate_decline_location.R` 读取 `input.rds`、公开文件和 `equal_records_samples.csv`，模拟下降发生在核心作业区以外、以内或两个区域（A2、L1、L2；900 个拟合），写出 `decline_location_fits.csv`、`decline_location_annual.csv` 和 `decline_location_truth.csv`。
27. `27_summarize_decline_location.R` 读取第 21、24 和 26 步的表，写出 `decline_location_errors.csv` 和 `decline_location_summary.csv`（补充材料 S11）。
28. `28_plot_decline_location.R` 读取 `decline_location_summary.csv`，写出 `FigS14_decline_location`（图 S14）。
29. `29_simulate_catchability_trend.R` 读取 `input.rds` 和公开文件，模拟情景 Q1 至 Q6（1,200 个拟合），写出 `catchability_trend_fits.csv`、`catchability_trend_annual.csv` 和 `catchability_trend_count_checks.csv`。
30. `30_summarize_catchability_trend.R` 读取第 21、23 和 29 步的表，写出 `catchability_trend_summary.csv`（表 S12、图 6b）。
31. `31_compare_all_and_complete.R` 读取第 21、24 和 29 步的表，写出全部重复上和拟合全部被接受的重复上的中位数（`coverage_all_vs_complete.csv`，表 S10 及表 S12 最后一列）。
32. `32_plot_coverage_experiment.R` 读取第 21、22、23 和 30 步的表，写出 `Fig6_coverage_experiment`（图 6）和 `FigS13_single_replicate`（图 S13）。

表 1 列出 `simulation_functions.R` 中定义的情景，没有单独的表文件。每幅图都写出 PDF 和 600 dpi 的 PNG。

模拟步骤 21、24、26 和 29 从配置文件读取重复次数（100）。快速运行时可指定较小的次数，例如 `--replicates 10`。各次
重复默认在 `--workers 8` 个进程上并行运行。每次重复沿用正式分析的随机种子，因此结果与进程数无关。每次重复保存为
`output/replicates/` 下的一个文件，中断后可接着运行；代码、配置或输入数据改动后需删除该文件夹。

## 命令

按顺序运行全部 32 个步骤，每步在新的 R 进程中按上述设置运行：

```sh
Rscript analysis_walkthrough/run_all.R --workers 8
```

`--workers` 传给模拟步骤 21、24、26 和 29。运行在第一个失败的步骤处停止，并指出该步骤。逐步运行时：

```sh
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
Rscript analysis_walkthrough/01_check_public_data.R
Rscript analysis_walkthrough/02_draw_study_design.R
Rscript analysis_walkthrough/03_prepare_records.R
Rscript analysis_walkthrough/04_define_areas.R
Rscript analysis_walkthrough/05_map_catch_cells.R
Rscript analysis_walkthrough/06_count_catch_cells.R
Rscript analysis_walkthrough/07_fit_baseline_models.R
Rscript analysis_walkthrough/08_spain_area_indices.R
Rscript analysis_walkthrough/09_compare_annual_indices.R
Rscript analysis_walkthrough/10_count_annual_overlap.R
Rscript analysis_walkthrough/11_plot_annual_indices.R
Rscript analysis_walkthrough/12_check_residuals.R
Rscript analysis_walkthrough/13_retrospective.R --workers 5
Rscript analysis_walkthrough/14_quarterly_change.R
Rscript analysis_walkthrough/15_compare_cells_fished_by_both.R
Rscript analysis_walkthrough/16_describe_catch_composition.R
Rscript analysis_walkthrough/17_fit_composition_models.R --workers 4
Rscript analysis_walkthrough/18_summarize_composition.R
Rscript analysis_walkthrough/19_plot_composition_indices.R
Rscript analysis_walkthrough/20_simulation_inputs.R
Rscript analysis_walkthrough/21_simulate_coverage.R --workers 8
Rscript analysis_walkthrough/22_summarize_coverage.R
Rscript analysis_walkthrough/23_coverage_change_intervals.R
Rscript analysis_walkthrough/24_simulate_equal_records.R --workers 8
Rscript analysis_walkthrough/25_summarize_equal_records.R
Rscript analysis_walkthrough/26_simulate_decline_location.R --workers 8
Rscript analysis_walkthrough/27_summarize_decline_location.R
Rscript analysis_walkthrough/28_plot_decline_location.R
Rscript analysis_walkthrough/29_simulate_catchability_trend.R --workers 8
Rscript analysis_walkthrough/30_summarize_catchability_trend.R
Rscript analysis_walkthrough/31_compare_all_and_complete.R
Rscript analysis_walkthrough/32_plot_coverage_experiment.R
```

## 核对

每一步最后读取其新表以及 `results/directed_fishery/` 中同名的存档表（第 03、04、21 和 29 步另读取
`tmp/directed_fishery/` 中存档的模型输入和核对表），输出"一致"，或不一致的列及最大差值。数值的相对差小于约 1e-8
时视为一致，这可吸收 R 写表时把数字保留 15 位有效数字造成的差别。对每幅 PNG 图，该步报告文件与存档图完全相同、像素
相同，或不同像素所占比例。重复次数较少的运行与存档结果中相同的重复比较，其汇总预期会不同。

生产脚本以 NumPy 默认随机数生成器进行第 22、25、27 和 30 步的自助抽样。`common.R` 中的 `numpy_rng()` 在基础 R 中
重复该生成器，因此相同种子给出相同的抽样和相同的区间。

每一步都已按上述方式运行和核对：每张表都与存档表一致（最大绝对差为 2.5e-12，出现在 `model_diagnostics.csv`），
每幅 PNG 图都与存档文件完全相同。补充地图（`shared_cell_cpue_atlas.pdf`）未比较，
因为 PDF 记录了生成时间。表 S2、S5 和 S6 没有单独的存档表：第 12 步由参数和诊断文件输出表 S2，第 16 步把表 S5 与
`results/directed_fishery/composition_figure_receipt*.json` 中保存的候选方案比较，表 S6 即 `model_diagnostics.csv`
的 `delta_AIC_vs_common` 列。

## 未包括的内容

- 捕捞策略模拟和图 S10：该存档模拟使用 sdmTMB 1.0.0（结果见 `results/sdmtmb_full/` 和 `results/paper_figures/main_simulation_contrasts.csv`）。
- 图册和 Word 文件：只汇编已有的图和文字，不做计算。
- 重新下载公开文件：IATTC 每次发布都会替换该文件，因此第 01 步核对 2026-11-06 下载的保存副本。
- 去掉空间随机场的重新拟合：没有存档拟合需要这样做，因此若有拟合未被接受，第 07 和 17 步给出说明后停止。
