# =====================================================================================================
# Run all 32 steps in order
# 按顺序运行全部 32 个步骤
#
# Runs steps 01 to 32, each in a new Rscript process with the settings of the README, so that no step depends
# on objects left by another. The simulation steps 21, 24, 26 and 29 receive the number of worker processes;
# steps 13 and 17 use their own defaults (5 and 4). The run stops at the first step that fails and names it.
# Each step prints its own comparison with the archived results.
# 依次运行第 01 至 32 步，每步在新的 Rscript 进程中按 README 的设置运行，因此任何步骤都不依赖其他步骤留下的对象。
# 模拟步骤 21、24、26 和 29 接收工作进程数；第 13 和 17 步使用各自的默认值（5 和 4）。运行在第一个失败的步骤处停止，
# 并指出该步骤。每个步骤输出各自与存档结果的比较。
#
# Reads / 读取:  the files read by each step / 各步骤读取的文件
# Writes / 写出: the files written by each step, all in analysis_walkthrough/output/ / 各步骤写出的文件，均在
#                analysis_walkthrough/output/ 中
#
# Run from the project folder (the folder that holds data/, config/ and results/).
# 在项目文件夹（包含 data/、config/ 和 results/ 的文件夹）中运行。
#   Rscript analysis_walkthrough/run_all.R --workers 8
# =====================================================================================================

if (!file.exists("config/directed_fishery.json")) stop("Run this from the project folder / 请在项目文件夹中运行")
argv <- commandArgs(trailingOnly=TRUE)
workers <- if ("--workers" %in% argv) as.integer(argv[match("--workers", argv) + 1L]) else 8L
stopifnot(!is.na(workers), workers >= 1)

# The settings of the README: one core per process and UTF-8 for the Chinese folder names.
# README 中的设置：每个进程一个核心，并以 UTF-8 读取中文文件夹名。
Sys.setenv(LANG="en_US.UTF-8", LC_ALL="en_US.UTF-8", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1")

# The 32 steps in order; the list must match the step files in this folder.
# 按顺序排列的 32 个步骤；该清单必须与本文件夹中的步骤文件一致。
steps <- c("01_check_public_data", "02_draw_study_design", "03_prepare_records", "04_define_areas",
           "05_map_catch_cells", "06_count_catch_cells", "07_fit_baseline_models", "08_spain_area_indices",
           "09_compare_annual_indices", "10_count_annual_overlap", "11_plot_annual_indices", "12_check_residuals",
           "13_retrospective", "14_quarterly_change", "15_compare_cells_fished_by_both",
           "16_describe_catch_composition", "17_fit_composition_models", "18_summarize_composition",
           "19_plot_composition_indices", "20_simulation_inputs", "21_simulate_coverage", "22_summarize_coverage",
           "23_coverage_change_intervals", "24_simulate_equal_records", "25_summarize_equal_records",
           "26_simulate_decline_location", "27_summarize_decline_location", "28_plot_decline_location",
           "29_simulate_catchability_trend", "30_summarize_catchability_trend", "31_compare_all_and_complete",
           "32_plot_coverage_experiment")
simulation <- c("21_simulate_coverage", "24_simulate_equal_records", "26_simulate_decline_location",
                "29_simulate_catchability_trend")
files <- sub("\\.R$", "", list.files("analysis_walkthrough", pattern="^[0-9]{2}_.*\\.R$"))
if (!identical(sort(files, method="radix"), steps)) stop("The step list does not match the step files / 步骤清单与步骤文件不一致")

# Run each step in a new process and stop at the first failure.
# 在新进程中运行每个步骤，并在第一次失败时停止。
rscript <- file.path(R.home("bin"), "Rscript")
for (step in steps) {
  args <- c(file.path("analysis_walkthrough", paste0(step, ".R")), if (step %in% simulation) c("--workers", workers))
  cat(sprintf("\n==== %s ====\n", paste(c("Rscript", args), collapse=" ")))
  status <- system2(rscript, args)
  if (status != 0) {
    cat(sprintf("\nStep %s failed (exit status %d); the later steps were not run.\n", substr(step, 1, 2), status))
    cat(sprintf("第 %s 步失败（退出状态 %d）；之后的步骤未运行。\n", substr(step, 1, 2), status))
    quit(save="no", status=1)
  }
}
cat("\nAll 32 steps finished / 全部 32 个步骤已完成\n")
