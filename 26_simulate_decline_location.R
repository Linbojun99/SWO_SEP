# =====================================================================================================
# Step 26. Simulation of the location of the decline (Section 2.4, Section S11)
# 第 26 步  下降位置的模拟（第 2.4 节，补充材料 S11）
#
# Uses the regional abundance trajectory of scenario A2 in each replicate and places the decline outside the
# core area only (A2), inside it only (L1) or in both areas at the same proportional rate (L2), so that the
# true regional change is the same in the three cases. Each case is fitted to Spanish records alone and to
# the two designs with the same number of records of step 24, with the same record selection (900 fits).
# These are simulation results, not estimates from the fishery data.
# 采用每次重复中情景 A2 的区域丰度轨迹，并把下降分别放在仅核心作业区以外（A2）、仅核心作业区以内（L1）或两个区域按
# 相同比例下降（L2），使三种情形的真实区域变化相同。每种情形分别拟合仅西班牙记录，以及第 24 步中记录数相同的两种
# 设计，记录选择相同（900 个拟合）。这些是模拟结果，不是渔业数据的估计。
#
# Reads / 读取:  analysis_walkthrough/output/input.rds (step 04), equal_records_samples.csv (step 24),
#                config/directed_fishery.json, data/raw/longline/PublicLLTunaBillfishNum.csv
# Writes / 写出: analysis_walkthrough/output/decline_location_fits.csv, decline_location_annual.csv,
#                decline_location_truth.csv (regional, core-area and outside abundance / 区域、核心作业区及以外的丰度),
#                replicates/decline_location/rep_<replicate>.rds
#
# Run from the project folder; the first line makes R use one core per process and read the Chinese folder
# names. Add --replicates 10 for a quick run.
# 在项目文件夹中运行；第一行让 R 每个进程只用一个核心，并能读取中文文件夹名。快速运行时加上
# --replicates 10。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/26_simulate_decline_location.R --workers 8
# =====================================================================================================

cat("Step 26: simulation of the location of the decline (A2, L1, L2)\n")
cat("第 26 步：下降位置的模拟（A2、L1、L2）\n")
source("analysis_walkthrough/simulation_functions.R")
locations <- c("periphery_decline","core_decline","uniform_decline")

# Draw the same record selections as step 24 and check them against the records written there.
# 抽取与第 24 步相同的记录选择，并与该步写出的记录核对。
samples <- equal_record_samples(d, cfg$seed + 4100L)
selected <- read.csv(file.path(out,"equal_records_samples.csv"))
stopifnot(identical(as.character(d$Flag),selected$Flag), identical(d$Year,selected$Year),
          identical(d$cell,selected$cell), all(d$Hooks == selected$Hooks))

# Fit Spanish records alone and the two designs with the same number of records in each of the three cases.
# 在三种情形下分别拟合仅西班牙记录以及记录数相同的两种设计。
parts <- cached_replicates("decline_location", seq_along(samples), function(i) {
  mixed <- which(!is.na(selected[[sprintf("core_id_rep_%02d",i)]]))
  rows <- c(list(Spain=design_rows$Spain),lapply(samples[[i]],sort))
  stopifnot(identical(rows$Core,which(selected$core_design)), identical(rows$Mixed,mixed))
  ans <- experiment(i, run=locations, sampling_rows=rows)
  message(sprintf("Replicate %03d: %d/%d accepted",i,sum(ans$fits$status=="accepted"),nrow(ans$fits)))
  ans
})
collect <- function(name) do.call(rbind,lapply(parts,"[[",name))
fits <- collect("fits")
annual <- collect("annual")
truth <- collect("area_truth")
stopifnot(nrow(fits)==n_rep*3L*3L, !anyDuplicated(fits[c("rep_id","scenario","design")]),
          !any(fits$status=="error"), nrow(annual)==nt*sum(fits$status=="accepted"),
          nrow(truth)==n_rep*3L*nt, !anyDuplicated(truth[c("rep_id","scenario","year")]))
write.csv(fits,file.path(out,"decline_location_fits.csv"),row.names=FALSE)
write.csv(annual,file.path(out,"decline_location_annual.csv"),row.names=FALSE)
write.csv(truth,file.path(out,"decline_location_truth.csv"),row.names=FALSE)
print(table(fits$scenario, fits$design, fits$status))

# ---- Compare with the archived results / 与存档结果比较 ----
source("analysis_walkthrough/common.R")
same_reps <- function(name) subset(read.csv(file.path(archive, name)), rep_id <= n_rep)
for (name in c("decline_location_fits.csv","decline_location_annual.csv","decline_location_truth.csv"))
  compare_table(name, old=same_reps(name))
