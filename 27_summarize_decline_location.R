# =====================================================================================================
# Step 27. Location of the decline and the regional index (Section S11, Fig. S14)
# 第 27 步  下降位置与区域指数（补充材料 S11，图 S14）
#
# Checks that the true regional abundance is the same in each year of a replicate in the three cases (A2, L1,
# L2) and that the A2 fits equal those of steps 21 and 24. For each accepted fit it calculates the MAE of the
# log index and the error of the estimated change. Panel (a) of Fig. S14 uses the replicates in which all nine
# fits were accepted; panel (b) uses those in which the six fits of the two designs with the same number of
# records were accepted. Medians have 95% intervals from 10,000 bootstrap samples of whole replicates (NumPy
# generator with seed cfg$seed + 4500, repeated by numpy_rng() in common.R). These are simulation results, not
# estimates from the fishery data.
# 核对每次重复中三种情形（A2、L1、L2）各年的真实区域丰度相同，且 A2 的拟合与第 21 和 24 步相同。对每个被接受的
# 拟合计算对数指数 MAE 和估计变化误差。图 S14 的 (a) 使用九个拟合均被接受的重复；(b) 使用记录数相同的两种设计的
# 六个拟合均被接受的重复。中位数的 95% 区间来自对整次重复进行的 10,000 次自助抽样（NumPy 生成器，种子
# cfg$seed + 4500，由 common.R 中的 numpy_rng() 重复）。这些是模拟结果，不是渔业数据的估计。
#
# Reads / 读取:  analysis_walkthrough/output/decline_location_fits.csv, decline_location_annual.csv,
#                decline_location_truth.csv (step 26), coverage_simulation.csv (step 21),
#                equal_records_fits.csv (step 24), config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/decline_location_errors.csv, decline_location_summary.csv
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/27_summarize_decline_location.R
# =====================================================================================================

source("analysis_walkthrough/common.R")
cat("Step 27: location of the decline (Section S11) / 第 27 步：下降位置（补充材料 S11）\n")
cfg <- jsonlite::fromJSON("config/directed_fishery.json")
scenarios <- c("periphery_decline", "core_decline", "uniform_decline")
designs <- c("Spain", "Core", "Mixed")
years <- seq(2009, 2023, 2)
keys <- c("rep_id", "scenario", "design")
fits <- read.csv(file.path(out, "decline_location_fits.csv"))
annual <- read.csv(file.path(out, "decline_location_annual.csv"))
areas <- read.csv(file.path(out, "decline_location_truth.csv"))
n_rep <- length(unique(fits$rep_id))
expected <- expand.grid(rep_id=seq_len(n_rep), scenario=scenarios, design=designs, stringsAsFactors=FALSE)
stopifnot(nrow(fits) == nrow(expected), setequal(row_key(fits, keys), row_key(expected, keys)),
          all(fits$records == c(Spain=384, Core=1510, Mixed=1510)[fits$design]))
accepted <- fits[fits$status == "accepted", ]
stopifnot(all(accepted$gradient < cfg$gradient_limit),
          setequal(row_key(annual, keys), row_key(accepted, keys)))

# ---- True abundance / 真实丰度 ----
# Regional abundance is the same in the three cases; the core area is constant in A2, the area outside it is
# constant in L1, and both areas decline by the same proportion in L2.
# 三种情形的区域丰度相同；A2 中核心作业区不变，L1 中核心作业区以外不变，L2 中两个区域按相同比例下降。
target <- sort_rows(unique(areas[c("rep_id", "year")]), c("rep_id", "year"))
for (scenario in scenarios) {
  one <- areas[areas$scenario == scenario, ]
  target[[scenario]] <- one$regional[match(row_key(target, c("rep_id", "year")), row_key(one, c("rep_id", "year")))]
}
stopifnot(max(abs(as.matrix(target[scenarios])/target$periphery_decline - 1), na.rm=TRUE) < 1e-12,
          all_close(areas$core + areas$outside, areas$regional, rtol=1e-12))
for (scenario in scenarios) {
  path <- areas[areas$scenario == scenario, ]
  first <- path[order(path$year), ]
  first <- first[!duplicated(first$rep_id), ]
  at <- match(path$rep_id, first$rep_id)
  core_scale <- path$core/first$core[at]; outside_scale <- path$outside/first$outside[at]
  stopifnot(switch(scenario, uniform_decline=all_close(core_scale, outside_scale, rtol=1e-12),
                   periphery_decline=all_close(core_scale, rep(1, length(core_scale)), rtol=1e-12),
                   core_decline=all_close(outside_scale, rep(1, length(outside_scale)), rtol=1e-12)))
}

# The A2 fits are those of steps 21 (Spanish records alone) and 24 (the two designs with the same number of records).
# A2 的拟合与第 21 步（仅西班牙记录）和第 24 步（记录数相同的两种设计）相同。
original <- rbind(subset(read.csv(file.path(out, "coverage_simulation.csv")), design == "Spain")[c(keys, "status", "estimate")],
                  read.csv(file.path(out, "equal_records_fits.csv"))[c(keys, "status", "estimate")])
original <- sort_rows(subset(original, scenario == "periphery_decline"), keys)
repeat_fits <- sort_rows(subset(fits, scenario == "periphery_decline"), keys)
stopifnot(identical(row_key(original, keys), row_key(repeat_fits, keys)), identical(original$status, repeat_fits$status),
          all_close(original$estimate, repeat_fits$estimate, rtol=1e-7, equal_nan=TRUE))

# ---- Errors of each fit / 各拟合的误差 ----
truth <- data.frame(target[c("rep_id", "year")], regional=target$periphery_decline)
early <- truth$year %in% years[1:2]
reference <- tapply(truth$regional[early], truth$rep_id[early], mean)
truth$scaled <- truth$regional/reference[as.character(truth$rep_id)]
joined <- annual
joined$scaled <- truth$scaled[match(row_key(annual, c("rep_id", "year")), row_key(truth, c("rep_id", "year")))]
joined <- joined[!is.na(joined$scaled), ]
stopifnot(all_close(joined$true_index, joined$scaled, rtol=1e-12))
joined$year_error <- abs(log(joined$estimated_index/joined$scaled))
metrics <- accepted[c(keys, "estimate", "truth")]
rownames(metrics) <- NULL
metrics$year_error <- tapply(joined$year_error, row_key(joined, keys), mean)[row_key(metrics, keys)]
metrics$period_error <- abs(log(metrics$estimate/metrics$truth))
write.csv(metrics, file.path(out, "decline_location_errors.csv"), row.names=FALSE)

# ---- Medians and bootstrap intervals (Fig. S14) / 中位数与自助法区间（图 S14） ----
rng <- numpy_rng(cfg$seed + 4500)
rows <- list()
for (cohort in c("three_designs", "equal_records")) {
  required <- if (cohort == "three_designs") designs else c("Core", "Mixed")
  counts <- table(accepted$rep_id[accepted$design %in% required])
  ids <- sort(as.numeric(names(counts)[counts == length(required)*length(scenarios)]))
  cat(sprintf("  %s: %d replicates with all fits accepted / %d 次重复的拟合全部被接受\n", cohort, length(ids), length(ids)))
  draws <- bootstrap_draws(rng, length(ids))
  record <- function(scenario, series, measure, values) {
    interval <- bootstrap_interval(values, draws)
    rows[[length(rows) + 1L]] <<- data.frame(cohort=cohort, scenario=scenario, series=series, measure=measure,
      n=length(ids), median=median(values), lower=interval[1], upper=interval[2],
      improved=if (series == "Core_minus_Mixed") sum(values > 0) else NA_integer_)
  }
  for (scenario in scenarios) {
    frame <- metrics[metrics$scenario == scenario, ]
    column <- function(design, name) { x <- frame[frame$design == design, ]; x[[name]][match(ids, x$rep_id)] }
    if (cohort == "three_designs") {
      record(scenario, "Truth", "change_percent", 100*(column("Spain", "truth") - 1))
      for (design in designs) record(scenario, design, "change_percent", 100*(column(design, "estimate") - 1))
    }
    for (measure in c("year_error", "period_error")) {
      for (design in required) record(scenario, design, measure, column(design, measure))
      if (cohort == "equal_records") {
        record(scenario, "Core_minus_Mixed", measure, column("Core", measure) - column("Mixed", measure))
      }
    }
  }
}
summary_table <- do.call(rbind, rows)
write.csv(summary_table, file.path(out, "decline_location_summary.csv"), row.names=FALSE)
cat(sprintf("  Accepted %d of %d fits / %d 个拟合中 %d 个被接受\n", nrow(accepted), nrow(fits), nrow(fits), nrow(accepted)))

# ---- Compare with the archived results / 与存档结果比较 ----
if (n_rep < cfg$coverage_replicates) {
  cat(sprintf("  Only %d replicates were run, so the summary differs from the archived run. / 只运行了 %d 次重复，因此汇总与存档结果不同。\n",
              n_rep, n_rep))
}
compare_table("decline_location_errors.csv", keys=keys)
compare_table("decline_location_summary.csv")
