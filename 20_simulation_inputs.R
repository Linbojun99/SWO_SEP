# =====================================================================================================
# Step 20. Simulation inputs: sampling designs and catch levels (Section 2.4, Tables 2 and S7)
# 第 20 步  模拟输入：取样设计与渔获水平（第 2.4 节，表 2 和表 S7）
#
# Counts, for each of the three nested sampling designs, the records fitted and the cells sampled in and
# outside Spain's core area (Table 2). Compares the relative catchability of each fleet in the simulation with
# its nominal CPUE (total swordfish divided by total hooks) in the public domain records of the simulated
# cells and years, and gives the median expected count at mean density (Table S7).
# 对三种嵌套取样设计分别统计拟合的记录数和在西班牙核心作业区内外取样的网格数（表 2）。将模拟中各船队的相对可捕
# 性与其在模拟网格和年份的公开记录中的名义 CPUE（剑鱼总尾数除以总钩数）比较，并给出平均密度下期望尾数的中位数
# （表 S7）。
#
# Reads / 读取:  analysis_walkthrough/output/input.rds (step 04), config/directed_fishery.json,
#                data/raw/longline/PublicLLTunaBillfishNum.csv
# Writes / 写出: analysis_walkthrough/output/coverage_sampling_designs.csv, simulated_catch_levels.csv
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/20_simulation_inputs.R
# =====================================================================================================

cat("Step 20: simulation inputs (Tables 2 and S7) / 第 20 步：模拟输入（表 2 和表 S7）\n")
source("analysis_walkthrough/simulation_functions.R")

# ---- Table 2: records and cells of each sampling design / 表 2：各取样设计的记录数与网格数 ----
designs <- do.call(rbind,lapply(names(design_rows),function(design) {
  train <- d[design_rows[[design]],]
  data.frame(design,records=nrow(train),sampled_core=length(unique(train$cell[train$spanish])),
             sampled_outside=length(unique(train$cell[!train$spanish])),
             spanish_outside_records=sum(train$Flag=="ESP" & !train$spanish),
             prediction_cells=nrow(g))
}))
print(designs)
write.csv(designs,file.path(out,"coverage_sampling_designs.csv"),row.names=FALSE)

# ---- Table S7: catchability and nominal CPUE / 表 S7：可捕性与名义 CPUE ----
# Sum hooks and retained swordfish of the eligible records by cell, year and fleet in the simulated cells and years.
# 在模拟网格和年份内，按网格、年份和船队汇总合格记录的钩数和剑鱼留存尾数。
catch <- aggregate(cbind(Hooks,SWOn)~cell+Year+Flag,subset(z$d,Year %in% years & cell %in% g$cell),sum)
phi <- 5
# At mean density (1) the expected count of a record is catchability x hooks/1000, and the NB2 coefficient of
# variation of a count is sqrt(1/mu + 1/phi).
# 在平均密度（1）下，一条记录的期望尾数为可捕性 × 钩数/1000；NB2 尾数的变异系数为 sqrt(1/mu + 1/phi)。
catch_levels <- do.call(rbind,lapply(cfg$flags,function(f) {
  x <- catch[catch$Flag==f,]
  mu <- q[[f]]*x$Hooks/1000
  m <- median(mu)
  data.frame(fleet=f,records=nrow(x),hooks=sum(x$Hooks),swordfish=sum(x$SWOn),
             nominal_cpue=sum(x$SWOn)/sum(x$Hooks)*1000,simulated_catchability=q[[f]],
             median_expected_count=m,share_expected_below_5=mean(mu<5),
             cv_at_median=sqrt(1/m+1/phi),cv_lower_limit=sqrt(1/phi))
}))
print(catch_levels,row.names=FALSE)
write.csv(catch_levels,file.path(out,"simulated_catch_levels.csv"),row.names=FALSE)

# ---- Compare with the archived results / 与存档结果比较 ----
source("analysis_walkthrough/common.R")
for (name in c("coverage_sampling_designs.csv","simulated_catch_levels.csv")) compare_table(name)
