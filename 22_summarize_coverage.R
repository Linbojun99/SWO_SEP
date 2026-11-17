# =====================================================================================================
# Step 22. Errors of the three nested sampling designs (Section 2.4, Table S9)
# 第 22 步  三种嵌套取样设计的误差（第 2.4 节，表 S9）
#
# Uses the replicates in which all three fits were accepted. For each scenario and design it gives the median
# error of the estimated change (absolute difference between the logarithms of the estimated and true change
# ratios), the median estimated and true change, and the median MAE of the log index over the eight years
# (each series scaled to its early mean). For the decline and recovery scenario it also gives the mean index
# in 2015 and 2017 divided by the early mean. The 95% intervals of the medians come from 10,000 bootstrap
# samples of whole replicates within each scenario, so that the designs remain paired (NumPy generator with
# seed cfg$seed + 3200, repeated by numpy_rng() in common.R; scenarios in the order A2, A1, A3). These are
# simulation results, not estimates from the fishery data.
# 使用三个拟合均被接受的重复。对每个情景和设计，给出估计变化误差（估计与真实变化比值的对数之差的绝对值）的中位数、
# 估计与真实变化的中位数，以及八个年份对数指数 MAE 的中位数（每条序列都除以其早期均值）。对下降后回升情景，另给出
# 2015 和 2017 年指数均值除以早期均值。中位数的 95% 区间来自在每个情景内对整次重复进行的 10,000 次自助抽样，使各
# 设计保持配对（NumPy 生成器，种子 cfg$seed + 3200，由 common.R 中的 numpy_rng() 重复；情景顺序为 A2、A1、A3）。
# 这些是模拟结果，不是渔业数据的估计。
#
# Reads / 读取:  analysis_walkthrough/output/coverage_simulation.csv, coverage_annual_intervals.csv (step 21),
#                config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/coverage_summary.csv (Table S9 / 表 S9),
#                coverage_annual_paired_summary.csv (Fig. 6a / 图 6a), coverage_paired_errors.csv
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/22_summarize_coverage.R
# =====================================================================================================

source("analysis_walkthrough/common.R")
cat("Step 22: errors of the three nested sampling designs (Table S9) / 第 22 步：三种嵌套取样设计的误差（表 S9）\n")
cfg <- jsonlite::fromJSON("config/directed_fishery.json")

# Read the fits of step 21: three designs for each scenario and replicate, all with the same true change.
# 读取第 21 步的拟合：每个情景和重复有三种设计，真实变化相同。
sim <- read.csv(file.path(out, "coverage_simulation.csv"))
n_rep <- length(unique(sim$rep_id))
stopifnot(nrow(sim) == n_rep*9, setequal(sim$rep_id, seq_len(n_rep)),
          all(tapply(sim$truth, row_key(sim, c("rep_id", "scenario")), function(x) length(unique(x))) == 1))
fit_keys <- c("scenario", "rep_id", "design")
stopifnot(setequal(sim$design, c("Spain", "Shared", "Extended")), !anyDuplicated(row_key(sim, fit_keys)))

# Keep the replicates in which all three fits were accepted, and calculate the error of the estimated change.
# 保留三个拟合均被接受的重复，并计算估计变化的误差。
replicate_key <- row_key(sim, c("rep_id", "scenario"))
ok <- ave(sim$status == "accepted", replicate_key, FUN=all)
paired <- sim[ok, ]
paired$abs_error <- abs(log(paired$estimate/paired$truth))

# Median annual estimate, interquartile range and median truth by year over the same replicates (Fig. 6a).
# 在相同重复上按年份计算年度估计的中位数、四分位距和真实值的中位数（图 6a）。
annual <- read.csv(file.path(out, "coverage_annual_intervals.csv"))
stopifnot(setequal(annual$year, seq(2009, 2023, 2)), all(table(row_key(annual, fit_keys)) == 8))
complete <- annual[row_key(annual, fit_keys) %in% row_key(paired, fit_keys), ]
stopifnot(nrow(complete) == 8*nrow(paired))
band_keys <- c("scenario", "design", "year")
bands <- do.call(rbind, lapply(split(complete, row_key(complete, band_keys)), function(g)
  data.frame(g[1, band_keys], pairs=nrow(g), median_estimate=median(g$estimated_index),
             q25=unname(quantile(g$estimated_index, .25, type=7)), q75=unname(quantile(g$estimated_index, .75, type=7)),
             median_truth=median(g$true_index))))
bands <- sort_rows(bands, band_keys)
write.csv(bands, file.path(out, "coverage_annual_paired_summary.csv"), row.names=FALSE)

# MAE of the log index over the eight years, and the mean index in 2015 and 2017 relative to the early mean.
# 八个年份对数指数的 MAE，以及 2015 和 2017 年指数均值相对于早期均值的比值。
annual$log_error <- abs(log(annual$estimated_index/annual$true_index))
fit_key <- row_key(annual, fit_keys)
middle <- annual$year %in% c(2015, 2017)
paired_key <- row_key(paired, fit_keys)
paired$annual_abs_log_error <- tapply(annual$log_error, fit_key, mean)[paired_key]
paired$trough_ratio <- tapply(annual$estimated_index[middle], fit_key[middle], mean)[paired_key]
paired$true_trough_ratio <- tapply(annual$true_index[middle], fit_key[middle], mean)[paired_key]
stopifnot(!anyNA(paired$annual_abs_log_error))
stats <- do.call(rbind, lapply(split(paired, row_key(paired, c("scenario", "design"))), function(g)
  data.frame(scenario=g$scenario[1], design=g$design[1], pairs=nrow(g), median_error=median(g$abs_error),
             median_trend=median(g$estimate), median_truth=median(g$truth),
             annual_abs_log_error=median(g$annual_abs_log_error))))
stats <- sort_rows(stats, c("scenario", "design"))

# 10,000 bootstrap samples of whole replicates within each scenario, shared by the three designs.
# 在每个情景内对整次重复进行 10,000 次自助抽样，三种设计共用。
rng <- numpy_rng(cfg$seed + 3200)
rows <- list()
for (scenario in c("periphery_decline", "stable", "periphery_fluctuation")) {
  one <- paired[paired$scenario == scenario, ]
  ids <- sort(unique(one$rep_id))
  draws <- bootstrap_draws(rng, length(ids))
  for (design in sorted(one$design)) {
    group <- one[one$design == design, ]
    group <- group[match(ids, group$rep_id), ]
    row <- data.frame(scenario=scenario, design=design)
    measures <- c("annual_abs_log_error", if (scenario == "periphery_fluctuation") "trough_ratio")
    for (measure in measures) {
      interval <- bootstrap_interval(group[[measure]], draws)
      row[[paste0(measure, "_lower")]] <- interval[1]; row[[paste0(measure, "_upper")]] <- interval[2]
    }
    if (scenario == "periphery_fluctuation") {
      row$trough_ratio <- median(group$trough_ratio); row$true_trough_ratio <- median(group$true_trough_ratio)
    } else {
      row[c("trough_ratio_lower", "trough_ratio_upper", "trough_ratio", "true_trough_ratio")] <- NA_real_
    }
    rows[[length(rows) + 1L]] <- row
  }
}
intervals <- do.call(rbind, rows)
stopifnot(nrow(intervals) == nrow(stats))
stats <- cbind(stats, intervals[match(row_key(stats, c("scenario", "design")), row_key(intervals, c("scenario", "design"))),
                                setdiff(names(intervals), c("scenario", "design"))])
stats <- stats[c("scenario", "design", "pairs", "median_error", "median_trend", "median_truth", "annual_abs_log_error",
                 "annual_abs_log_error_lower", "annual_abs_log_error_upper", "trough_ratio", "trough_ratio_lower",
                 "trough_ratio_upper", "true_trough_ratio")]
rownames(stats) <- NULL
write.csv(stats, file.path(out, "coverage_summary.csv"), row.names=FALSE)
print(stats[c("scenario", "design", "pairs", "median_error", "annual_abs_log_error")], row.names=FALSE)

# Errors of each complete replicate side by side, and the reduction in error from adding records outside the core area.
# 并列各完整重复的误差，以及加入核心作业区以外记录后误差的减少量。
piv <- sort_rows(unique(paired[c("scenario", "rep_id")]), c("scenario", "rep_id"))
side <- function(design, column, scenarios=unique(paired$scenario)) {
  x <- paired[paired$design == design & paired$scenario %in% scenarios, ]
  x[[column]][match(row_key(piv, c("scenario", "rep_id")), row_key(x, c("scenario", "rep_id")))]
}
designs <- c("Extended", "Shared", "Spain")
for (d in designs) piv[[d]] <- side(d, "abs_error")
for (mode in c("Spain", "Shared")) piv[[paste0("gain_over_", mode)]] <- piv[[mode]] - piv$Extended
for (d in designs) piv[[paste0("annual_", d)]] <- side(d, "annual_abs_log_error")
for (d in designs) piv[[paste0("trough_", d)]] <- side(d, "trough_ratio", "periphery_fluctuation")
write.csv(piv, file.path(out, "coverage_paired_errors.csv"), row.names=FALSE)

# ---- Compare with the archived results / 与存档结果比较 ----
if (n_rep < cfg$coverage_replicates) {
  cat(sprintf("  Only %d replicates were run, so these tables differ from the archived run of %d. / 只运行了 %d 次重复，因此这些表与存档的 %d 次重复不同。\n",
              n_rep, cfg$coverage_replicates, n_rep, cfg$coverage_replicates))
}
compare_table("coverage_summary.csv", keys=c("scenario", "design"))
compare_table("coverage_annual_paired_summary.csv", keys=c("scenario", "design", "year"))
compare_table("coverage_paired_errors.csv", keys=c("scenario", "rep_id"))
