# =====================================================================================================
# Step 09. Rank correlations with Spain's index and precision of the annual indices (Tables S1 and S3)
# 第 09 步  与西班牙指数的秩相关及年度指数的精度（表 S1 和表 S3）
#
# (1) pair_comparisons.csv (Table S1): Spearman's rank correlation between each tuna-targeting longline fleet's
#     index and Spain's index in their area of overlap (two-, three- and four-year thresholds), over the years
#     with estimates for both fleets, and over the years in which both fleets had records in the area of overlap.
# (2) index_precision.csv (Table S3): the standard error of the log index of each baseline index of Fig. 3,
#     recovered from its pointwise 95% interval, which approximates the coefficient of variation of a log-normal
#     index, with the number of records in the cells indexed.
# (1) pair_comparisons.csv（表 S1）：各金枪鱼延绳钓船队指数与西班牙指数在两者重叠区域（两年、三年和四年阈值）内的
#     Spearman 秩相关，分别基于两支船队都有估计的年份，以及两支船队在重叠区域内都有记录的年份。
# (2) index_precision.csv（表 S3）：由逐点 95% 区间还原的图 3 各基线指数的对数指数标准误，近似于对数正态指数的变异
#     系数，并列出所计算网格内的记录数。
#
# Reads / 读取:  analysis_walkthrough/output/indices_baseline.csv (step 07)
# Writes / 写出: analysis_walkthrough/output/pair_comparisons.csv, index_precision.csv
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/09_compare_annual_indices.R
# =====================================================================================================

source("analysis_walkthrough/common.R")
cat("Step 09: rank correlations and index precision (Tables S1 and S3) / 第 09 步：秩相关与指数精度（表 S1 和表 S3）\n")
indices <- read.csv(file.path(out, "indices_baseline.csv"))
baseline <- indices[!indices$annual_only & indices$region == "all" & indices$model == "all", ]

# Spearman's rank correlation: the correlation of the ranks, with tied values given their mean rank.
# Spearman 秩相关：秩的相关系数，相同数值取其平均秩。
spearman <- function(x, y) if (length(x) < 2) NA_real_ else suppressWarnings(cor(rank(x), rank(y)))

# ---- (1) Rank correlations with Spain (Table S1) / 与西班牙的秩相关（表 S1） ----
# Place Spain's and the fleet's indices side by side in each area of overlap over the years both have an estimate.
# 在每个重叠区域内，把西班牙与该船队的指数在两者都有估计的年份并排列出。
others <- baseline[baseline$domain != "ESP", ]
groups <- sort_rows(unique(others[c("domain", "threshold")]), c("domain", "threshold"))
pairs <- list()
for (k in seq_len(nrow(groups))) {
  flag <- groups$domain[k]; threshold <- groups$threshold[k]
  group <- others[others$domain == flag & others$threshold == threshold, ]
  stopifnot(setequal(group$flag, c("ESP", flag)))
  spain <- group[group$flag == "ESP", ]; fleet <- group[group$flag == flag, ]
  years <- sort(intersect(spain$year[!is.na(spain$index)], fleet$year[!is.na(fleet$index)]))
  comparison <- data.frame(year=years, ESP=spain$index[match(years, spain$year)],
                           fleet=fleet$index[match(years, fleet$year)])
  observed <- !is.na(spain$records[match(years, spain$year)]) & spain$records[match(years, spain$year)] > 0 &
    !is.na(fleet$records[match(years, fleet$year)]) & fleet$records[match(years, fleet$year)] > 0
  # All years with both estimates, and the years in which both fleets had records in the area of overlap.
  # 两者都有估计的全部年份，以及两支船队在重叠区域内都有记录的年份。
  for (mode in c("paired_years", "local_records")) {
    series <- if (mode == "paired_years") comparison else comparison[observed, ]
    pairs[[length(pairs) + 1L]] <- data.frame(flag=flag, threshold=threshold, years=nrow(series), support=mode,
      cells=group$domain_cells[1], rho=spearman(series$ESP, series$fleet))
  }
}
pairs <- do.call(rbind, pairs)
write.csv(pairs, file.path(out, "pair_comparisons.csv"), row.names=FALSE)
print(reshape(transform(pairs[c("flag", "support", "threshold")], rho=round(pairs$rho, 2)), idvar=c("flag", "support"),
              timevar="threshold", direction="wide"), row.names=FALSE)

# ---- (2) Precision of the annual indices (Table S3) / 年度指数的精度（表 S3） ----
# The interval is exp(log(est) -/+ 1.96 se), so (log(upper) - log(lower)) / (2 x 1.96) recovers the standard error se of the log index.
# 区间为 exp(log(est) -/+ 1.96 se)，因此 (log(upper) - log(lower)) / (2 × 1.96) 还原对数指数的标准误 se。
z <- 1.959964
rows <- baseline[baseline$threshold == 3, ]
stopifnot(!anyDuplicated(rows[c("domain", "flag", "year")]))
rows$sd_log <- (log(rows$upper) - log(rows$lower))/(2*z)
precision <- sort_rows(rows[c("domain", "flag", "year", "records", "index", "lower", "upper", "sd_log")],
                       c("domain", "flag", "year"))
write.csv(precision, file.path(out, "index_precision.csv"), row.names=FALSE)
summary_rows <- lapply(split(precision, row_key(precision, c("domain", "flag"))), function(g)
  data.frame(domain=g$domain[1], flag=g$flag[1], years=nrow(g), median_sd_log=median(g$sd_log),
             min_sd_log=min(g$sd_log), max_sd_log=max(g$sd_log), median_records=median(g$records)))
summary_table <- sort_rows(do.call(rbind, summary_rows), c("domain", "flag"))
summary_table[4:7] <- round(summary_table[4:7], 2)
print(summary_table, row.names=FALSE)
own <- precision[precision$domain == precision$flag & precision$domain != "ESP", ]
own_rho <- sapply(sorted(own$flag), function(f) with(own[own$flag == f, ], round(cor(sd_log, records, method="spearman"), 2)))
cat("Spearman correlation of records and standard error, own series / 记录数与标准误的 Spearman 相关（各船队自身序列）:",
    paste0(names(own_rho), " ", own_rho, collapse=", "), "\n")

# ---- Compare with the archived results / 与存档结果比较 ----
compare_table("pair_comparisons.csv")
compare_table("index_precision.csv")
