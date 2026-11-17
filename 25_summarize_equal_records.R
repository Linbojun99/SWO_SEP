# =====================================================================================================
# Step 25. Error with the same number of records from the core area or from both areas (Table S8)
# 第 25 步  记录数相同、记录来自核心作业区或两个区域时的误差（表 S8）
#
# Checks the record selections of step 24 (same Spanish records, same fleet and year for every replacement,
# numbers of hooks within a factor of 1.05) and counts the records, hooks and cells of the two designs in each
# replicate. For the replicates in which both fits were accepted it calculates, within each replicate, the MAE
# of the log index and the error of the estimated change of both designs and the reduction in error (core area
# only minus both areas). The median reduction has a 95% interval from 10,000 paired bootstrap samples (NumPy
# generator with seed cfg$seed + 4200, repeated by numpy_rng() in common.R). The fits of the core-area design
# must equal the second design of step 21. These are simulation results, not estimates from the fishery data.
# 核对第 24 步的记录选择（西班牙记录相同、每条替换记录的船队和年份相同、钩数相差不超过 1.05 倍），并统计每次重复
# 两种设计的记录数、钩数和网格数。对两个拟合均被接受的重复，在每次重复内计算两种设计的对数指数 MAE 和估计变化误差，
# 以及误差的减少量（仅核心作业区减去两个区域）。减少量中位数的 95% 区间来自 10,000 次配对自助抽样（NumPy 生成器，
# 种子 cfg$seed + 4200，由 common.R 中的 numpy_rng() 重复）。核心作业区设计的拟合必须与第 21 步的第二种设计相同。
# 这些是模拟结果，不是渔业数据的估计。
#
# Reads / 读取:  analysis_walkthrough/output/equal_records_samples.csv, equal_records_fits.csv,
#                equal_records_annual.csv (step 24), coverage_simulation.csv, coverage_annual_intervals.csv
#                (step 21), config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/equal_records_balance.csv, equal_records_pairs.csv,
#                equal_records_summary.csv (Table S8 / 表 S8)
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/25_summarize_equal_records.R
# =====================================================================================================

source("analysis_walkthrough/common.R")
cat("Step 25: error with the same number of records (Table S8) / 第 25 步：记录数相同时的误差（表 S8）\n")
cfg <- jsonlite::fromJSON("config/directed_fishery.json")
scenarios <- c("stable", "periphery_decline", "periphery_fluctuation")
designs <- c("Core", "Mixed")
years <- seq(2009, 2023, 2)
fits <- read.csv(file.path(out, "equal_records_fits.csv"))
annual <- read.csv(file.path(out, "equal_records_annual.csv"))
n_rep <- length(unique(fits$rep_id))

# ---- Check the record selections / 核对记录选择 ----
# In every replicate both designs keep all Spanish records, and each replacement has the same fleet and year
# as the core-area record it replaces and a number of hooks within a factor of 1.05.
# 每次重复中两种设计都保留全部西班牙记录，每条替换记录与被替换的核心作业区记录船队和年份相同，钩数相差不超过 1.05 倍。
records <- read.csv(file.path(out, "equal_records_samples.csv"))
stopifnot(!anyDuplicated(records[c("Flag", "cell", "Year")]), setequal(records$Year, years))
spanish <- records$Flag == "ESP"
core <- spanish | records$core_cell
stopifnot(identical(core, records$core_design), sum(startsWith(names(records), "core_id_rep_")) == n_rep)
rows <- list()
for (rep in seq_len(n_rep)) {
  original <- records[[sprintf("core_id_rep_%02d", rep)]]
  mixed <- !is.na(original)
  ids <- original[mixed]
  stopifnot(!anyDuplicated(ids), setequal(ids, records$row_id[core]))
  baseline <- records[match(ids, records$row_id), ]
  chosen <- records[mixed, ]
  stopifnot(identical(chosen$Flag, baseline$Flag), all(chosen$Year == baseline$Year))
  same <- chosen$row_id == ids
  stopifnot(all(same[chosen$Flag == "ESP"]), all(!chosen$core_cell[!same]),
            all(pmax(chosen$Hooks, baseline$Hooks) <= 1.05*pmin(chosen$Hooks, baseline$Hooks)),
            all(mixed[spanish]), sum(mixed) == sum(core))
  for (k in seq_along(designs)) {
    mask <- if (k == 1) core else mixed
    included <- records[mask, ]
    selected <- mask & !spanish
    rows[[length(rows) + 1L]] <- data.frame(rep_id=rep, design=designs[k], records=sum(mask),
      spanish_records=sum(mask & spanish), bycatch_records=sum(selected),
      bycatch_outside=sum(selected & !records$core_cell), bycatch_hooks=sum(as.numeric(records$Hooks[selected])),
      core_cells=length(unique(included$cell[included$core_cell])),
      outside_cells=length(unique(included$cell[!included$core_cell])), prediction_cells=141)
  }
}
balance <- do.call(rbind, rows)

# ---- Check the fits / 核对拟合 ----
# The core-area design is the second design of Table 2, so its fits must equal those of step 21.
# 核心作业区设计即表 2 的第二种设计，因此其拟合必须与第 21 步相同。
keys <- c("rep_id", "scenario", "design")
stopifnot(nrow(fits) == n_rep*length(scenarios)*length(designs), !anyDuplicated(row_key(fits, keys)))
at <- match(row_key(fits, c("rep_id", "design")), row_key(balance, c("rep_id", "design")))
stopifnot(!anyNA(at), all(fits$records == balance$records[at]))
main <- subset(read.csv(file.path(out, "coverage_simulation.csv")), design == "Shared")
core_fits <- subset(fits, design == "Core")
at <- match(row_key(core_fits, c("rep_id", "scenario")), row_key(main, c("rep_id", "scenario")))
stopifnot(nrow(core_fits) == n_rep*length(scenarios), !anyNA(at), all(core_fits$status == main$status[at]),
          all_close(core_fits$estimate, main$estimate[at], rtol=1e-5, equal_nan=TRUE))
main_paths <- subset(read.csv(file.path(out, "coverage_annual_intervals.csv")), design == "Shared")
new_paths <- subset(annual, design == "Core")
path_keys <- c("rep_id", "scenario", "year")
at <- match(row_key(new_paths, path_keys), row_key(main_paths, path_keys))
stopifnot(!anyDuplicated(row_key(new_paths, path_keys)), !anyDuplicated(row_key(main_paths, path_keys)),
          nrow(new_paths) == nrow(main_paths), !anyNA(at),
          all_close(new_paths$estimated_index, main_paths$estimated_index[at], rtol=1e-5))
valid <- fits[fits$status == "accepted", ]
stopifnot(all(valid$gradient < cfg$gradient_limit), nrow(annual) == 8*nrow(valid))
write.csv(balance, file.path(out, "equal_records_balance.csv"), row.names=FALSE)

# ---- Errors within each replicate and paired bootstrap / 每次重复内的误差与配对自助抽样 ----
rng <- numpy_rng(cfg$seed + 4200)
contrasts <- list(); summaries <- list()
for (scenario in scenarios) {
  metrics <- list()
  for (rep in seq_len(n_rep)) {
    pair <- fits[fits$rep_id == rep & fits$scenario == scenario, ]
    stopifnot(nrow(pair) == 2, length(unique(pair$truth)) == 1)
    if (!all(pair$status == "accepted")) next
    paths <- annual[annual$rep_id == rep & annual$scenario == scenario, ]
    path <- function(design, column) { x <- paths[paths$design == design, ]; x[[column]][match(years, x$year)] }
    truth <- sapply(designs, path, column="true_index"); est <- sapply(designs, path, column="estimated_index")
    stopifnot(identical(as.numeric(sorted(paths$year)), as.numeric(years)), identical(truth[, "Core"], truth[, "Mixed"]))
    end <- colMeans(est[7:8, ])/colMeans(est[1:2, ])
    true_end <- mean(truth[7:8, "Core"])/mean(truth[1:2, "Core"])
    stopifnot(all_close(true_end, pair$truth[1], rtol=1e-5, atol=1e-8))
    # MAE of the log index (annual) and error of the estimated change (period).
    # 对数指数 MAE（annual）和估计变化误差（period）。
    errors <- list(annual=colMeans(abs(log(est/truth))), period=abs(log(end/true_end)))
    for (name in names(errors)) {
      metrics[[length(metrics) + 1L]] <- data.frame(scenario=scenario, rep_id=rep, measure=name,
        true_period_ratio=true_end, core_period_ratio=end[["Core"]], mixed_period_ratio=end[["Mixed"]],
        core_error=errors[[name]][["Core"]], mixed_error=errors[[name]][["Mixed"]],
        gain=errors[[name]][["Core"]] - errors[[name]][["Mixed"]])
    }
  }
  frame <- do.call(rbind, metrics)
  contrasts <- c(contrasts, metrics)
  n_pairs <- length(unique(frame$rep_id))
  draws <- bootstrap_draws(rng, n_pairs)
  for (name in unique(frame$measure)) {
    group <- frame[frame$measure == name, ]
    group <- group[order(group$rep_id), ]
    interval <- bootstrap_interval(group$gain, draws)
    summaries[[length(summaries) + 1L]] <- data.frame(scenario=scenario, measure=name, pairs=n_pairs,
      core_error=median(group$core_error), mixed_error=median(group$mixed_error), paired_gain=median(group$gain),
      lower=interval[1], upper=interval[2], improved=sum(group$gain > 0),
      true_period_ratio=median(group$true_period_ratio), core_period_ratio=median(group$core_period_ratio),
      mixed_period_ratio=median(group$mixed_period_ratio))
  }
}
write.csv(do.call(rbind, contrasts), file.path(out, "equal_records_pairs.csv"), row.names=FALSE)
summary_table <- do.call(rbind, summaries)
write.csv(summary_table, file.path(out, "equal_records_summary.csv"), row.names=FALSE)
print(summary_table[c("scenario", "measure", "pairs", "core_error", "mixed_error", "paired_gain", "lower", "upper",
                      "improved")], row.names=FALSE)

# ---- Compare with the archived results / 与存档结果比较 ----
if (n_rep < cfg$coverage_replicates) {
  cat(sprintf("  Only %d replicates were run, so the summary differs from the archived run. / 只运行了 %d 次重复，因此汇总与存档结果不同。\n",
              n_rep, n_rep))
}
compare_table("equal_records_balance.csv")
compare_table("equal_records_pairs.csv")
compare_table("equal_records_summary.csv")
