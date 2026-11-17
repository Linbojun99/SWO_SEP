# =====================================================================================================
# Step 29. Simulation of an increase in the catchability of the tuna-targeting longline fleets (Section 2.4, Section S12)
# 第 29 步  金枪鱼延绳钓船队可捕性增加的模拟（第 2.4 节，补充材料 S12）
#
# Repeats the three abundance scenarios with the catchability of the five tuna-targeting longline fleets
# multiplied by (1 + r)^(year - 2009) for r = 0.01 and 0.02 (scenarios Q1 to Q6), the annual increases in
# longline catchability of the 2024 benchmark assessment of bigeye tuna (Xu et al., 2024b). Spain's
# catchability, the replicate seeds and the uniform draws are those of step 21, so the Spanish catches do not
# change and r = 0 is the run of step 21. The second and third designs of Table 2 are refitted with the
# estimation model of step 21, whose fleet effects are constant (1,200 fits). These are simulation results,
# not estimates from the fishery data.
# 在三种丰度情景下，把五支金枪鱼延绳钓船队的可捕性乘以 (1 + r)^(年份 - 2009)，r = 0.01 和 0.02（情景 Q1 至 Q6），
# 即 2024 年东太平洋大眼金枪鱼基准评估所用的延绳钓可捕性年增幅（Xu et al., 2024b）。西班牙的可捕性、各重复的随机
# 种子和均匀随机数与第 21 步相同，因此西班牙渔获不变，r = 0 即第 21 步的模拟。表 2 的第二和第三种设计以第 21 步的
# 估计模型重新拟合，其船队效应恒定（1,200 个拟合）。这些是模拟结果，不是渔业数据的估计。
#
# Reads / 读取:  analysis_walkthrough/output/input.rds (step 04), config/directed_fishery.json,
#                data/raw/longline/PublicLLTunaBillfishNum.csv
# Writes / 写出: analysis_walkthrough/output/catchability_trend_fits.csv, catchability_trend_annual.csv,
#                catchability_trend_count_checks.csv (Spanish counts unchanged / 西班牙尾数不变),
#                replicates/catchability_trend/rep_<replicate>.rds
#
# Run from the project folder; the first line makes R use one core per process and read the Chinese folder
# names. Add --replicates 10 for a quick run.
# 在项目文件夹中运行；第一行让 R 每个进程只用一个核心，并能读取中文文件夹名。快速运行时加上
# --replicates 10。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/29_simulate_catchability_trend.R --workers 8
# =====================================================================================================

cat("Step 29: simulation of an annual increase in catchability (Q1 to Q6)\n")
cat("第 29 步：可捕性逐年增加的模拟（Q1 至 Q6）\n")
source("analysis_walkthrough/simulation_functions.R")
rates <- c(0.01, 0.02)
rows <- design_rows[c("Shared", "Extended")]
esp <- d$Flag == "ESP"
first <- d$Year == min(years)

# For each replicate, generate the counts at r = 0 without fitting, then at r = 0.01 and 0.02 with fitting.
# Spanish counts and all 2009 counts must equal those at r = 0; later counts of the other fleets cannot fall.
# 对每次重复，先在 r = 0 下生成尾数而不拟合，再在 r = 0.01 和 0.02 下生成并拟合。西班牙尾数和 2009 年的全部尾数
# 必须与 r = 0 相同；其他船队之后年份的尾数不会减少。
fit_one <- function(i) {
  base <- experiment(i, run=scenarios, sampling_rows=list())
  stopifnot(is.null(base$fits), identical(names(base$counts), scenarios))
  parts <- lapply(rates, function(r) {
    ans <- experiment(i, run=scenarios, sampling_rows=rows, q_trend=r)
    check <- do.call(rbind, lapply(scenarios, function(s) {
      new <- ans$counts[[s]]; old <- base$counts[[s]]
      stopifnot(identical(new[esp], old[esp]), identical(new[first], old[first]), all(new >= old))
      data.frame(rep_id=i, scenario=s, annual_q_increase=r, esp_rows=sum(esp),
                 esp_counts_identical=TRUE, nonesp_rows=sum(!esp),
                 nonesp_rows_changed=sum(new[!esp] != old[!esp]),
                 nonesp_count_ratio=sum(new[!esp]) / sum(old[!esp]))
    }))
    ans$fits$annual_q_increase <- r
    if (!is.null(ans$annual)) ans$annual$annual_q_increase <- r
    ans$counts <- NULL
    c(ans, list(count_check=check))
  })
  message(sprintf("Replicate %03d: %d/%d accepted", i,
                  sum(vapply(parts, function(x) sum(x$fits$status == "accepted"), 0)),
                  sum(vapply(parts, function(x) nrow(x$fits), 0))))
  parts
}
parts <- unlist(cached_replicates("catchability_trend", seq_len(n_rep), fit_one), recursive=FALSE)
collect <- function(name) do.call(rbind, lapply(parts, "[[", name))
fits <- collect("fits")
annual <- collect("annual")
check <- collect("count_check")
# 12 fits and 6 count checks per replicate.
# 每次重复 12 个拟合和 6 项尾数核对。
stopifnot(nrow(fits) == n_rep * 12L, !anyDuplicated(fits[c("rep_id", "scenario", "design", "annual_q_increase")]),
          !any(fits$status == "error"), nrow(annual) == nt * sum(fits$status == "accepted"),
          nrow(check) == n_rep * 6L, all(check$esp_counts_identical))
fits$uncorrected <- NULL
write.csv(fits, file.path(out, "catchability_trend_fits.csv"), row.names=FALSE)
write.csv(annual, file.path(out, "catchability_trend_annual.csv"), row.names=FALSE)
write.csv(check, file.path(out, "catchability_trend_count_checks.csv"), row.names=FALSE)
print(table(paste(fits$scenario, fits$annual_q_increase), fits$design, fits$status))

# ---- Compare with the archived results / 与存档结果比较 ----
source("analysis_walkthrough/common.R")
same_reps <- function(path) subset(read.csv(path), rep_id <= n_rep)
for (name in c("catchability_trend_fits.csv","catchability_trend_annual.csv"))
  compare_table(name, old=same_reps(file.path(archive, name)))
compare_table("catchability_trend_count_checks.csv", old=same_reps("tmp/directed_fishery/catchability_trend_count_checks.csv"))
