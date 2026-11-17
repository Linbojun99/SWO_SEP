# =====================================================================================================
# Step 31. Medians over all replicates and over the replicates with all fits accepted (Tables S10 and S12)
# 第 31 步  全部重复与拟合全部被接受的重复上的中位数（表 S10 和表 S12）
#
# For the designs with records from more than one fleet, compares the median estimated change and the median
# MAE of the log index over all accepted fits of a design with those over the replicates in which all fits
# being compared were accepted: all three designs of step 21 (Table S10), both designs of step 24, and all seven
# fits of step 30 (last column of Table S12). These are simulation results, not estimates from the fishery data.
# 对使用多支船队记录的设计，比较某设计全部被接受拟合上的估计变化中位数和对数指数 MAE 中位数，与所比较拟合全部被
# 接受的重复上的相应中位数：第 21 步的三种设计（表 S10）、第 24 步的两种设计，以及第 30 步的七个拟合（表 S12 最后
# 一列）。这些是模拟结果，不是渔业数据的估计。
#
# Reads / 读取:  analysis_walkthrough/output/coverage_simulation.csv, coverage_annual_intervals.csv (step 21),
#                equal_records_fits.csv, equal_records_annual.csv (step 24), catchability_trend_fits.csv,
#                catchability_trend_annual.csv (step 29)
# Writes / 写出: analysis_walkthrough/output/coverage_all_vs_complete.csv
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/31_compare_all_and_complete.R
# =====================================================================================================

source("analysis_walkthrough/common.R")
cat("Step 31: medians over all and complete replicates (Table S10) / 第 31 步：全部与完整重复上的中位数（表 S10）\n")
scenarios <- c("stable", "periphery_decline", "periphery_fluctuation")
fit_keys <- c("rep_id", "scenario", "design")
table_of <- function(name) read.csv(file.path(out, paste0(name, ".csv")))

# Estimated and true change in percent and MAE of the log index of each accepted fit.
# 每个被接受拟合以百分比表示的估计与真实变化，以及对数指数 MAE。
fit_metrics <- function(fits, annual, keys) {
  m <- fits[fits$status == "accepted", c(keys, "estimate", "truth")]
  mae <- tapply(abs(log(annual$estimated_index/annual$true_index)), row_key(annual, keys), mean)
  m$annual_abs_log_error <- mae[row_key(m, keys)]
  stopifnot(!anyNA(m$annual_abs_log_error), length(mae) == nrow(m))
  m$change_pct <- 100*(m$estimate - 1)
  m$true_change_pct <- 100*(m$truth - 1)
  m
}

# Replicates in which all the named designs were accepted, by scenario.
# 按情景列出指定设计全部被接受的重复。
complete_ids <- function(fits, need) {
  x <- fits[fits$design %in% need, ]
  lapply(split(x, x$scenario), function(s) {
    accepted <- tapply(s$status == "accepted", s$rep_id, sum); size <- table(s$rep_id)[names(accepted)]
    as.numeric(names(accepted)[accepted == length(need) & size == length(need)])
  })
}

rows <- list()
add <- function(run, scenario, design, subset, definition, frame, rate=NA_real_) {
  rows[[length(rows) + 1L]] <<- data.frame(run=run, scenario=scenario, annual_q_increase=rate, design=design,
    subset=subset, definition=definition, replicates=length(unique(frame$rep_id)),
    median_true_change_pct=median(frame$true_change_pct), median_change_pct=median(frame$change_pct),
    median_annual_abs_log_error=median(frame$annual_abs_log_error))
}

# Three nested designs (step 21; Table S10).
# 三种嵌套设计（第 21 步；表 S10）。
main <- table_of("coverage_simulation")
spain <- main[main$design == "Spain" & main$status == "accepted", ]
spain_ok <- split(spain$rep_id, spain$scenario)
mm <- fit_metrics(main, table_of("coverage_annual_intervals"), fit_keys)
main_complete <- complete_ids(main, c("Spain", "Shared", "Extended"))
for (s in scenarios) {
  for (d in c("Shared", "Extended")) {
    one <- mm[mm$scenario == s & mm$design == d, ]
    add("main", s, d, "all", "all accepted fits of the design", one)
    add("main", s, d, "complete", "replicates with Spain, Shared and Extended accepted",
        one[one$rep_id %in% main_complete[[s]], ])
  }
}

# Two designs with the same number of records (step 24).
# 记录数相同的两种设计（第 24 步）。
er <- table_of("equal_records_fits")
em <- fit_metrics(er, table_of("equal_records_annual"), fit_keys)
er_complete <- complete_ids(er, c("Core", "Mixed"))
for (s in scenarios) {
  for (d in c("Core", "Mixed")) {
    one <- em[em$scenario == s & em$design == d, ]
    add("equal_records", s, d, "all", "all accepted fits of the design", one)
    add("equal_records", s, d, "complete", "replicates with Core and Mixed accepted (Table S8 cohort)",
        one[one$rep_id %in% er_complete[[s]], ])
    add("equal_records", s, d, "spain_accepted",
        "replicates with Core and Mixed accepted and the main-run Spain-only fit accepted",
        one[one$rep_id %in% intersect(er_complete[[s]], spain_ok[[s]]), ])
  }
}

# Increasing catchability (step 29): complete replicates have all seven fits of step 30 accepted.
# 可捕性增加（第 29 步）：完整重复指第 30 步的七个拟合均被接受。
ct <- table_of("catchability_trend_fits")
ckeys <- c(fit_keys, "annual_q_increase")
cm <- fit_metrics(ct, table_of("catchability_trend_annual"), ckeys)
main$annual_q_increase <- 0
status <- rbind(main[c(ckeys, "status")], ct[c(ckeys, "status")])
status <- status[status$status == "accepted", ]
required <- c(paste("Spain", 0), paste(rep(c("Shared", "Extended"), each=3), c(0, 0.01, 0.02)))
ct_complete <- lapply(split(status, status$scenario), function(one) {
  have <- tapply(paste(one$design, one$annual_q_increase), one$rep_id, function(x) all(required %in% x))
  as.numeric(names(have)[have])
})
for (s in scenarios) {
  for (r in c(0.01, 0.02)) {
    for (d in c("Shared", "Extended")) {
      one <- cm[cm$scenario == s & cm$design == d & cm$annual_q_increase == r, ]
      add("catchability_trend", s, d, "all", "all accepted fits of the design at this rate", one, r)
      add("catchability_trend", s, d, "complete",
          "replicates with all seven fits accepted (Spain; Shared and Extended at 0, 1 and 2%)",
          one[one$rep_id %in% ct_complete[[s]], ], r)
    }
  }
}
avc <- do.call(rbind, rows)
write.csv(avc, file.path(out, "coverage_all_vs_complete.csv"), row.names=FALSE)
print(avc[avc$run == "main", c("scenario", "design", "subset", "replicates", "median_change_pct",
                               "median_annual_abs_log_error")], row.names=FALSE)

# ---- Compare with the archived results / 与存档结果比较 ----
compare_table("coverage_all_vs_complete.csv")
