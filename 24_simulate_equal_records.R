# =====================================================================================================
# Step 24. Simulation of two designs with the same number of records (Section 2.4, Section S8)
# 第 24 步  记录数相同的两种设计的模拟（第 2.4 节，补充材料 S8）
#
# Compares two designs with the same 384 Spanish records and the same number of records from the tuna-
# targeting longline fleets. One adds all 1,126 core-area records of these fleets (the second design of
# Table 2); in the other, each of those records is replaced, where possible, with an unused record from
# outside the core area for the same fleet and year with a similar number of hooks. The selection is drawn
# once per replicate before the counts are generated and used for all three abundance scenarios (600 fits).
# These are simulation results, not estimates from the fishery data.
# 比较两种设计：两者都含相同的 384 条西班牙记录和相同数目的金枪鱼延绳钓船队记录。一种加入这些船队在核心作业区的
# 全部 1,126 条记录（表 2 的第二种设计）；另一种尽可能把其中每条记录替换为同一船队同一年份、核心作业区以外尚未
# 使用且钩数相近的记录。每次重复在生成尾数之前抽取一次记录选择，并用于全部三种丰度情景（600 个拟合）。这些是模拟
# 结果，不是渔业数据的估计。
#
# Reads / 读取:  analysis_walkthrough/output/input.rds (step 04), config/directed_fishery.json,
#                data/raw/longline/PublicLLTunaBillfishNum.csv
# Writes / 写出: analysis_walkthrough/output/equal_records_samples.csv (the records of both designs in each
#                replicate / 每次重复两种设计的记录), equal_records_fits.csv, equal_records_annual.csv,
#                replicates/equal_records/rep_<replicate>.rds
#
# Run from the project folder; the first line makes R use one core per process and read the Chinese folder
# names. Add --replicates 10 for a quick run.
# 在项目文件夹中运行；第一行让 R 每个进程只用一个核心，并能读取中文文件夹名。快速运行时加上
# --replicates 10。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/24_simulate_equal_records.R --workers 8
# =====================================================================================================

cat("Step 24: simulation of two designs with the same number of records\n")
cat("第 24 步：记录数相同的两种设计的模拟\n")
source("analysis_walkthrough/simulation_functions.R")

# Draw the record selections of all replicates from one random stream (seed cfg$seed + 4100), before fitting.
# 在拟合之前，从同一个随机数流（种子 cfg$seed + 4100）抽取全部重复的记录选择。
samples <- equal_record_samples(d, cfg$seed + 4100L)

# Fit the two designs (Core: core-area records; Mixed: records from both areas) in each replicate.
# 在每次重复中拟合两种设计（Core：核心作业区记录；Mixed：两个区域的记录）。
parts <- cached_replicates("equal_records", seq_along(samples), function(i) {
  ans <- experiment(i, sampling_rows=lapply(samples[[i]], sort))
  message(sprintf("Replicate %03d: %d/6 accepted", i, sum(ans$fits$status == "accepted")))
  ans
})
fits <- do.call(rbind, lapply(parts, "[[", "fits"))
annual <- do.call(rbind, lapply(parts, "[[", "annual"))
stopifnot(nrow(fits) == n_rep * length(scenarios) * 2L,
          nrow(annual) == nt * sum(fits$status == "accepted"))

# Write the selected records: for each replicate, the core-area record that each record of the Mixed design replaced.
# 写出所选记录：对每次重复，给出 Mixed 设计中每条记录所替换的核心作业区记录。
selections <- d[c("Flag","Year","cell","Hooks","spanish")]
names(selections)[5] <- "core_cell"
selections$row_id <- seq_len(nrow(d))
selections$core_design <- seq_len(nrow(d)) %in% samples[[1]]$Core
for (i in seq_along(samples)) {
  original_id <- rep(NA_integer_, nrow(d))
  original_id[samples[[i]]$Mixed] <- samples[[i]]$Core
  selections[[sprintf("core_id_rep_%02d", i)]] <- original_id
}
write.csv(selections, file.path(out, "equal_records_samples.csv"), row.names=FALSE)
write.csv(fits, file.path(out, "equal_records_fits.csv"), row.names=FALSE)
write.csv(annual, file.path(out, "equal_records_annual.csv"), row.names=FALSE)
print(table(fits$scenario, fits$design, fits$status))

# ---- Compare with the archived results / 与存档结果比较 ----
source("analysis_walkthrough/common.R")
same_reps <- function(name) subset(read.csv(file.path(archive, name)), rep_id <= n_rep)
old_samples <- read.csv(file.path(archive, "equal_records_samples.csv"))
compare_table("equal_records_samples.csv", old=old_samples[names(selections)])
compare_table("equal_records_fits.csv", old=same_reps("equal_records_fits.csv"))
compare_table("equal_records_annual.csv", old=same_reps("equal_records_annual.csv"))
