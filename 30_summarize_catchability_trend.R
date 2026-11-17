# =====================================================================================================
# Step 30. Estimated change with an increasing catchability of the tuna-targeting longline fleets (Table S12)
# 第 30 步  金枪鱼延绳钓船队可捕性增加时的估计变化（表 S12）
#
# Combines the fits at r = 0.01 and 0.02 (step 29) with the fits at constant catchability and the Spain-only
# fits of step 21, whose generating model is the r = 0 case. Over the replicates in which all seven fits of a
# scenario were accepted (Spain; + core area and + outside core area at 0, 1 and 2%) it gives the median true
# and estimated change, the median MAE of the log index and the median error of the estimated change. The
# intervals at r > 0 come from 10,000 paired bootstrap samples of the replicates (NumPy generator with seed
# cfg$seed + 5100, repeated by numpy_rng() in common.R); the rows at constant catchability take their intervals
# from step 23. These are simulation results, not estimates from the fishery data.
# 把 r = 0.01 和 0.02 的拟合（第 29 步）与第 21 步可捕性恒定的拟合及仅西班牙的拟合合并，后者的生成模型即 r = 0 的
# 情形。在某情景七个拟合（西班牙；+ 核心作业区和 + 核心作业区以外在 0、1 和 2% 下）均被接受的重复上，给出真实与估计
# 变化的中位数、对数指数 MAE 的中位数和估计变化误差的中位数。r > 0 的区间来自对重复进行的 10,000 次配对自助抽样
# （NumPy 生成器，种子 cfg$seed + 5100，由 common.R 中的 numpy_rng() 重复）；可捕性恒定的各行采用第 23 步的区间。
# 这些是模拟结果，不是渔业数据的估计。
#
# Reads / 读取:  analysis_walkthrough/output/catchability_trend_fits.csv, catchability_trend_annual.csv,
#                catchability_trend_count_checks.csv (step 29), coverage_simulation.csv,
#                coverage_annual_intervals.csv (step 21), coverage_change_intervals.csv (step 23),
#                config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/catchability_trend_summary.csv (Table S12 and Fig. 6b / 表 S12 和图 6b)
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/30_summarize_catchability_trend.R
# =====================================================================================================

source("analysis_walkthrough/common.R")
cat("Step 30: estimated change with increasing catchability (Table S12) / 第 30 步：可捕性增加时的估计变化（表 S12）\n")
cfg <- jsonlite::fromJSON("config/directed_fishery.json")
scenarios <- c("stable", "periphery_decline", "periphery_fluctuation")
rates <- c(0, 0.01, 0.02)
designs <- c("Spain", "Shared", "Extended")
keys <- c("rep_id", "scenario", "design", "annual_q_increase")

# Read the fits of steps 21 and 29; the r = 0 rows come from step 21.
# 读取第 21 和 29 步的拟合；r = 0 的各行来自第 21 步。
main <- read.csv(file.path(out, "coverage_simulation.csv")); main$annual_q_increase <- 0
main_annual <- read.csv(file.path(out, "coverage_annual_intervals.csv")); main_annual$annual_q_increase <- 0
main_intervals <- read.csv(file.path(out, "coverage_change_intervals.csv"))
trend <- read.csv(file.path(out, "catchability_trend_fits.csv"))
trend_annual <- read.csv(file.path(out, "catchability_trend_annual.csv"))
n_rep <- length(unique(trend$rep_id))
stopifnot(nrow(trend) == n_rep*length(scenarios)*2*2, !anyDuplicated(row_key(trend, keys)))

# The generating model differs from step 21 only in the catchability of the five tuna-targeting longline fleets.
# 生成模型与第 21 步的差别只在五支金枪鱼延绳钓船队的可捕性。
fit_keys <- c("rep_id", "scenario", "design")
for (r in c(0.01, 0.02)) {
  one <- trend[trend$annual_q_increase == r, ]
  at <- match(row_key(one, fit_keys), row_key(main, fit_keys))
  stopifnot(!anyNA(at), identical(as.numeric(one$truth), as.numeric(main$truth[at])),
            identical(as.numeric(one$records), as.numeric(main$records[at])))
}
counts <- read.csv(file.path(out, "catchability_trend_count_checks.csv"))
stopifnot(nrow(counts) == n_rep*length(scenarios)*2, all(counts$esp_counts_identical))

# MAE of the log index of each accepted fit, and the estimated and true change in percent.
# 每个被接受拟合的对数指数 MAE，以及以百分比表示的估计与真实变化。
columns <- c(keys, "status", "estimate", "truth")
fits <- rbind(main[columns], trend[columns])
annual <- rbind(main_annual[c(keys, "year", "estimated_index", "true_index")],
                trend_annual[c(keys, "year", "estimated_index", "true_index")])
metrics <- fits[fits$status == "accepted", c(keys, "estimate", "truth")]
log_error <- abs(log(annual$estimated_index/annual$true_index))
metrics$annual_abs_log_error <- tapply(log_error, row_key(annual, keys), mean)[row_key(metrics, keys)]
stopifnot(!anyNA(metrics$annual_abs_log_error))
metrics$change_pct <- 100*(metrics$estimate - 1)
metrics$true_change_pct <- 100*(metrics$truth - 1)
metrics$abs_change_error_pp <- abs(metrics$change_pct - metrics$true_change_pct)
metrics$abs_log_period_error <- abs(log(metrics$estimate/metrics$truth))

# Seven fits per replicate and scenario: Spain, and + core area and + outside core area at 0, 1 and 2%.
# 每次重复和情景七个拟合：西班牙，以及 + 核心作业区和 + 核心作业区以外在 0、1 和 2% 下。
required <- c(paste("Spain", 0), paste(rep(c("Shared", "Extended"), each=3), rates))
rng <- numpy_rng(cfg$seed + 5100)
rows <- list()
for (scenario in scenarios) {
  one <- metrics[metrics$scenario == scenario, ]
  have <- tapply(paste(one$design, one$annual_q_increase), one$rep_id, function(x) all(required %in% x))
  ids <- sort(as.numeric(names(have)[have]))
  cat(sprintf("  %s: %d replicates with all seven fits accepted / %d 次重复的七个拟合均被接受\n", scenario, length(ids), length(ids)))
  draws <- bootstrap_draws(rng, length(ids))
  for (r in rates) {
    for (design in designs) {
      source_rate <- if (design == "Spain") 0 else r
      g <- one[one$design == design & one$annual_q_increase == source_rate, ]
      g <- g[match(ids, g$rep_id), ]
      values <- g$change_pct
      if (source_rate == 0) {
        # Rows at constant catchability: the same replicates and medians as step 23, and its intervals.
        # 可捕性恒定的各行：重复和中位数与第 23 步相同，并采用其区间。
        ref_row <- main_intervals[main_intervals$scenario == scenario & main_intervals$design == design, ]
        stopifnot(nrow(ref_row) == 1, ref_row$pairs == length(ids), abs(ref_row$change - median(values)) < 1e-9)
        interval <- c(ref_row$lower, ref_row$upper)
      } else interval <- bootstrap_interval(values, draws)
      rows[[length(rows) + 1L]] <- data.frame(scenario=scenario, annual_q_increase=r, design=design,
        source=if (source_rate == 0) "main run" else "catchability-trend run", replicates=length(ids),
        true_change_pct=median(g$true_change_pct), change_pct=median(values), change_lower=interval[1],
        change_upper=interval[2], median_abs_change_error_pp=median(g$abs_change_error_pp),
        median_abs_log_period_error=median(g$abs_log_period_error), annual_abs_log_error=median(g$annual_abs_log_error))
    }
  }
}
summary_table <- do.call(rbind, rows)
write.csv(summary_table, file.path(out, "catchability_trend_summary.csv"), row.names=FALSE)
print(summary_table[c("scenario", "annual_q_increase", "design", "change_pct", "change_lower", "change_upper",
                      "annual_abs_log_error")], row.names=FALSE)

# ---- Compare with the archived results / 与存档结果比较 ----
if (n_rep < cfg$coverage_replicates) {
  cat(sprintf("  Only %d replicates were run, so the summary differs from the archived run. / 只运行了 %d 次重复，因此汇总与存档结果不同。\n",
              n_rep, n_rep))
}
compare_table("catchability_trend_summary.csv")
