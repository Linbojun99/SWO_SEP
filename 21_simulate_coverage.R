# =====================================================================================================
# Step 21. Simulation of the three abundance scenarios under the three sampling designs (Section 2.4)
# 第 21 步  三种丰度情景在三种取样设计下的模拟（第 2.4 节）
#
# For each of 100 independent replicates of scenarios A1 to A3 (Table 1), generates catches once at the
# reported cells and numbers of hooks and fits the estimation model to the three nested sampling designs of
# Table 2 (900 fits). Each accepted fit gives the bias-corrected regional index over the 141 cells and its
# pointwise 95% interval, scaled to the early mean, beside the true regional abundance. These are simulation
# results, not estimates from the fishery data.
# 对情景 A1 至 A3（表 1）的 100 次独立重复，每次在报告的网格和钩数上生成一次渔获，并以估计模型拟合表 2 的三种嵌套
# 取样设计（900 个拟合）。每个被接受的拟合给出 141 个网格上经偏差校正的区域指数及其逐点 95% 区间（除以早期均值），
# 并与真实区域丰度并列。这些是模拟结果，不是渔业数据的估计。
#
# Reads / 读取:  analysis_walkthrough/output/input.rds (step 04), config/directed_fishery.json,
#                data/raw/longline/PublicLLTunaBillfishNum.csv
# Writes / 写出: analysis_walkthrough/output/coverage_simulation.csv (one row per fit / 每个拟合一行),
#                coverage_annual_intervals.csv (estimated and true annual index / 估计与真实年度指数),
#                coverage_uncorrected_check.csv (index without bias correction, kept as a check / 未经偏差校正的指数，作为核对保留),
#                replicates/main/rep_<replicate>.rds
#
# Run from the project folder; the first line makes R use one core per process and read the Chinese folder
# names. Add --replicates 10 for a quick run.
# 在项目文件夹中运行；第一行让 R 每个进程只用一个核心，并能读取中文文件夹名。快速运行时加上
# --replicates 10。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/21_simulate_coverage.R --workers 8
# =====================================================================================================

cat("Step 21: simulation of scenarios A1 to A3 under three sampling designs\n")
cat("第 21 步：情景 A1 至 A3 在三种取样设计下的模拟\n")
source("analysis_walkthrough/simulation_functions.R")

# Run the replicates; each one simulates the three scenarios and fits the three designs to each.
# 运行各次重复；每次模拟三个情景，并对每个情景拟合三种设计。
parts <- cached_replicates("main",seq_len(n_rep),function(i) {
  ans <- experiment(i)
  message(sprintf("Replicate %03d: %d/%d accepted",i,sum(ans$fits$status=="accepted"),nrow(ans$fits)))
  ans
})

# Combine the replicates: 9 fits per replicate and 8 annual rows per accepted fit.
# 合并各次重复：每次重复 9 个拟合，每个被接受的拟合 8 行年度结果。
total <- do.call(rbind,lapply(parts,"[[","fits"))
annual <- do.call(rbind,lapply(parts,"[[","annual"))
stopifnot(nrow(total)==n_rep*3*length(scenarios),!any(total$status=="error"),
          nrow(annual)==nt*sum(total$status=="accepted"))
write.csv(total[c("rep_id","scenario","design","status","uncorrected","estimate")],
          file.path(out,"coverage_uncorrected_check.csv"),row.names=FALSE)
total$uncorrected <- NULL
write.csv(total,file.path(out,"coverage_simulation.csv"),row.names=FALSE)
write.csv(annual,file.path(out,"coverage_annual_intervals.csv"),row.names=FALSE)
print(table(total$scenario,total$design,total$status))

# ---- Compare with the archived results / 与存档结果比较 ----
# A run with fewer replicates is compared with the same replicates of the archived run.
# 重复次数较少的运行与存档结果中相同的重复比较。
source("analysis_walkthrough/common.R")
same_reps <- function(path) subset(read.csv(path), rep_id <= n_rep)
compare_table("coverage_simulation.csv", old=same_reps(file.path(archive,"coverage_simulation.csv")))
compare_table("coverage_annual_intervals.csv", old=same_reps(file.path(archive,"coverage_annual_intervals.csv")))
compare_table("coverage_uncorrected_check.csv", old=same_reps("tmp/directed_fishery/coverage_uncorrected_check.csv"))
