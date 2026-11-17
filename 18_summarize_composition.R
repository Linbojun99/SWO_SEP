# =====================================================================================================
# Step 18. Change in the annual index after adjustment for catch composition (Section 2.3, Table 3)
# 第 18 步  按渔获组成校正后年度指数的变化（第 2.3 节，表 3）
#
# For each tuna-targeting longline fleet in its area of overlap with Spain, divides each index from a model with
# composition covariates (the adjusted index) by the baseline index fitted to the same records in the same year.
# The annual percentage change is 100 (I_a / I_b - 1); Table 3 gives, for the cluster and PC models with
# swordfish excluded from the composition, the median absolute percentage change over all years with estimates.
# 对每支金枪鱼延绳钓船队在其与西班牙的重叠区域内，把含组成协变量模型的指数（校正指数）除以同年以同样记录拟合的基线
# 指数。年度百分比变化为 100 (I_a / I_b - 1)；表 3 给出组成中排除剑鱼的聚类和 PC 模型在所有有估计年份上的绝对百分比
# 变化中位数。
#
# Reads / 读取:  analysis_walkthrough/output/indices.csv, fits.csv (step 17)
# Writes / 写出: analysis_walkthrough/output/composition_adjustment_summary.csv (Table 3 / 表 3),
#                composition_changes.csv
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/18_summarize_composition.R
# =====================================================================================================

source("analysis_walkthrough/common.R")
cat("Step 18: change after adjustment for catch composition (Table 3) / 第 18 步：按渔获组成校正后的变化（表 3）\n")
indices <- read.csv(file.path(out, "indices.csv"))
fits <- read.csv(file.path(out, "fits.csv"))
rows <- indices[!indices$annual_only & indices$region == "all", ]

# The baseline, cluster and PC models of each fleet (swordfish excluded) used the same records.
# 各船队的基线、聚类和 PC 模型（排除剑鱼）使用相同的记录。
adjust <- rows[rows$model %in% c("common", "pc_exclude", "cluster_exclude"), ]
keys <- c("flag", "threshold", "year")
wide <- sort_rows(unique(adjust[keys]), keys)
for (m in c("cluster_exclude", "common", "pc_exclude")) {
  one <- adjust[adjust$model == m, ]
  stopifnot(!anyDuplicated(row_key(one, keys)))
  wide[[m]] <- one$index[match(row_key(wide, keys), row_key(one, keys))]
}
main_fits <- fits[fits$id %in% adjust$id, ]
stopifnot(all(main_fits$status == "accepted"), all(table(main_fits$flag) == 3),
          all(tapply(main_fits$records, main_fits$flag, function(x) length(unique(x))) == 1))

# Table 3: median absolute percentage change at the three-year threshold.
# 表 3：三年阈值下绝对百分比变化的中位数。
main <- wide[wide$threshold == 3, ]
values <- as.matrix(main[c("cluster_exclude", "common", "pc_exclude")])
stopifnot(all(is.finite(values)), all(values > 0), all(adjust$flag == adjust$domain))
cluster_change <- abs((main$cluster_exclude/main$common - 1)*100)
pc_change <- abs((main$pc_exclude/main$common - 1)*100)
flags <- sorted(main$flag)
typical <- data.frame(flag=flags, years=as.vector(table(main$flag)[flags]),
                      cluster_median_abs_pct=as.vector(tapply(cluster_change, main$flag, median)[flags]),
                      pc_median_abs_pct=as.vector(tapply(pc_change, main$flag, median)[flags]))
write.csv(typical, file.path(out, "composition_adjustment_summary.csv"), row.names=FALSE)
print(transform(typical, cluster_median_abs_pct=round(cluster_median_abs_pct, 1), pc_median_abs_pct=round(pc_median_abs_pct, 1)),
      row.names=FALSE)

# Annual ratio, percentage change and absolute log ratio of each adjusted index, with the records in the area of
# overlap; a year without both indices has no ratio.
# 各校正指数的年度比值、百分比变化和绝对对数比，并附重叠区域内的记录数；没有两个指数的年份没有比值。
two <- rep(seq_len(nrow(wide)), each=2)
ratios <- data.frame(wide[two, keys], model=rep(c("cluster_exclude", "pc_exclude"), nrow(wide)),
                     ratio=as.vector(rbind(wide$cluster_exclude, wide$pc_exclude))/wide$common[two], row.names=NULL)
ratios <- ratios[!is.na(ratios$ratio), ]
ratios$change_pct <- (ratios$ratio - 1)*100
ratios$abs_log_ratio <- abs(log(ratios$ratio))
counts <- adjust[adjust$model == "common", c(keys, "records")]
stopifnot(!anyDuplicated(row_key(counts, keys)))
ratios$records <- counts$records[match(row_key(ratios, keys), row_key(counts, keys))]
ratios <- ratios[!is.na(ratios$records), ]
write.csv(ratios, file.path(out, "composition_changes.csv"), row.names=FALSE)

# ---- Compare with the archived results / 与存档结果比较 ----
compare_table("composition_adjustment_summary.csv")
compare_table("composition_changes.csv")
