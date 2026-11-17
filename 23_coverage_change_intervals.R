# =====================================================================================================
# Step 23. Median estimated change with bootstrap intervals (Section 2.4, Fig. 6b at constant catchability)
# 第 23 步  估计变化的中位数及其自助法区间（第 2.4 节，图 6b 中可捕性恒定的部分）
#
# Over the replicates in which all three fits were accepted, calculates the median true change and the median
# estimated change of each design, as percentage change from the early mean (2009 and 2011) to the mean of
# 2021 and 2023. The 95% intervals come from 10,000 bootstrap samples of whole replicates within each
# scenario, so that the designs remain paired (R seed cfg$seed + 3100, scenarios in the order A2, A1, A3).
# These are simulation results, not estimates from the fishery data.
# 在三个拟合均被接受的重复上，计算真实变化的中位数和各设计估计变化的中位数，即由早期均值（2009 和 2011 年）到
# 2021 和 2023 年均值的百分比变化。95% 区间来自在每个情景内对整次重复进行的 10,000 次自助抽样，使各设计保持配对
# （R 种子 cfg$seed + 3100，情景顺序为 A2、A1、A3）。这些是模拟结果，不是渔业数据的估计。
#
# Reads / 读取:  analysis_walkthrough/output/coverage_simulation.csv (step 21), coverage_paired_errors.csv
#                (step 22), config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/coverage_change_intervals.csv
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/23_coverage_change_intervals.R
# =====================================================================================================

out <- "analysis_walkthrough/output"
cfg <- jsonlite::fromJSON("config/directed_fishery.json")
cat("Step 23: median estimated change with bootstrap intervals / 第 23 步：估计变化的中位数及其自助法区间\n")

# Keep the fits of the replicates in which all three fits were accepted (those listed in coverage_paired_errors.csv).
# 保留三个拟合均被接受的重复的拟合（即 coverage_paired_errors.csv 列出的重复）。
trials <- read.csv(file.path(out,"coverage_simulation.csv"))
paired <- read.csv(file.path(out,"coverage_paired_errors.csv"))
selected <- merge(trials,paired[c("scenario","rep_id")],by=c("scenario","rep_id"))
stopifnot(all(selected$status=="accepted"),!anyDuplicated(selected[c("scenario","rep_id","design")]))

# One random stream for the three scenarios, in the order A2, A1, A3; the draws are shared by the designs.
# 三个情景共用一个随机数流，顺序为 A2、A1、A3；各设计共用同一组抽样。
set.seed(cfg$seed+3100L)
designs <- c("Spain","Shared","Extended")
boot_order <- c("periphery_decline","stable","periphery_fluctuation")
stopifnot(setequal(unique(selected$scenario),boot_order))
coverage_stats <- do.call(rbind,lapply(split(selected,selected$scenario)[boot_order],function(one) {
  ids <- sort(unique(one$rep_id))
  stopifnot(all(table(one$rep_id)==length(designs)),length(ids)>=2L)
  draws <- matrix(sample.int(length(ids),length(ids)*10000L,replace=TRUE),nrow=length(ids))
  truth <- unique(one[c("rep_id","truth")])
  stopifnot(nrow(truth)==length(ids))
  truth_row <- data.frame(scenario=one$scenario[1],design="Truth",pairs=length(ids),
                          change=100*(median(truth$truth)-1),lower=NA_real_,upper=NA_real_)
  estimates <- lapply(designs,function(design) {
    subset_design <- one[one$design==design,]
    values <- 100*(subset_design$estimate[match(ids,subset_design$rep_id)]-1)
    boot <- apply(matrix(values[draws],nrow=length(ids)),2,median)
    ci <- unname(quantile(boot,c(.025,.975)))
    data.frame(scenario=one$scenario[1],design,pairs=length(ids),change=median(values),
               lower=ci[1],upper=ci[2])
  })
  do.call(rbind,c(list(truth_row),estimates))
}))
print(coverage_stats,row.names=FALSE)
write.csv(coverage_stats,file.path(out,"coverage_change_intervals.csv"),row.names=FALSE)

# ---- Compare with the archived results / 与存档结果比较 ----
source("analysis_walkthrough/common.R")
compare_table("coverage_change_intervals.csv")
