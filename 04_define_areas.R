# =====================================================================================================
# Step 04. Define the core area, the areas of overlap and the cells by assessment area (Section 2.1, Table S1)
# 第 04 步  确定核心作业区、重叠区域及各评估分区的网格（第 2.1 节，表 S1）
#
# Counts, for each cell, the years to 2023 in which Spain and each fleet both reported effort. Spain's own
# count defines Spain's core area (effort in at least three years) and each pair's count defines the area of
# overlap between the two fleets (effort by both in at least three of the same years); two-year and four-year
# thresholds are kept for Table S1. Keeps the eligible records, marks the records used for the composition
# models, finds the 15 years with estimates for all six fleets, lists the 31 fits, lists the cells with
# positive swordfish counts (Fig. 2) and counts the cells in each assessment area of H2 (Table S1).
# 对每个网格，统计 2023 年以前西班牙与各船队同年都报告努力量的年数。西班牙自身的年数确定西班牙核心作业区（至少
# 三年有努力量），每对船队的年数确定两支船队的重叠区域（两者至少在三个相同年份都有努力量）；两年和四年阈值保留
# 供表 S1 使用。保留合格记录，标记组成模型所用的记录，找出六支船队都有估计的 15 个年份，列出 31 个拟合，列出有
# 剑鱼正尾数的网格（图 2），并统计 H2 各评估分区的网格数（表 S1）。
#
# Reads / 读取:  analysis_walkthrough/output/records.rds (step 03), config/directed_fishery.json,
#                config/targeting_public_h2.json
# Writes / 写出: analysis_walkthrough/output/support.csv, footprint.csv, tasks.csv, area_structure.csv,
#                input.rds (the model input of the later steps / 后续步骤的模型输入)
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/04_define_areas.R
# =====================================================================================================

out <- "analysis_walkthrough/output"
cfg <- jsonlite::fromJSON("config/directed_fishery.json")
cat("Step 04: defining the core area, the areas of overlap and the assessment areas\n")
cat("第 04 步：确定核心作业区、重叠区域与评估分区\n")
records <- readRDS(file.path(out, "records.rds"))
d <- records$d; grid <- records$grid

# ---- (1) Years with effort in each cell (support.csv) / 每个网格有努力量的年数 ----
# List the fleet, cell and year combinations to 2023 from all H2 records, so that the areas follow reported
# effort whether or not swordfish were reported.
# 用全部 H2 记录列出 2023 年以前的船队、网格和年份组合，使区域按报告的努力量确定，而不论是否报告了剑鱼。
available <- unique(d[d$Year<=cfg$support_end, c("Flag","cell","Year")])
sp <- subset(available, Flag=="ESP", select=c(cell,Year))

# For each cell, count the years in which Spain and the fleet both reported; for Spain these are the years Spain reported.
# 对每个网格，统计西班牙与该船队同年都有报告的年数；对西班牙而言即其报告的年数。
membership <- list()
for (f in cfg$flags) {
  common <- merge(sp, subset(available, Flag==f, select=c(cell,Year)), by=c("cell","Year"))
  counts <- table(factor(common$cell, levels=grid$cell))
  membership[[f]] <- data.frame(domain=f, cell=grid$cell, shared_years=as.integer(counts))
}
coverage <- do.call(rbind, membership)
write.csv(merge(coverage, grid[c("cell","LonC5","LatC5")], by="cell"), file.path(out, "support.csv"), row.names=FALSE)

# Cells of the core area (Spain) and of each area of overlap at the two-, three- and four-year thresholds (Table S1).
# 两年、三年和四年阈值下核心作业区（西班牙）和各重叠区域的网格数（表 S1）。
print(sapply(2:4, function(k) tapply(coverage$shared_years>=k, factor(coverage$domain, levels=cfg$flags), sum)))

# ---- (2) Eligible records and the records of the composition models / 合格记录与组成模型的记录 ----
# Keep the eligible records, those of fleet-years with a positive swordfish count.
# 保留合格记录，即有剑鱼正尾数的船队年份的记录。
d <- subset(d, eligible)

# Mark the records with at least one retained albacore, bigeye tuna, yellowfin tuna, blue marlin or striped
# marlin; the composition models and the baseline model fitted to the same records use these records.
# 标记至少留存一尾长鳍金枪鱼、大眼金枪鱼、黄鳍金枪鱼、蓝枪鱼或条纹四鳍旗鱼的记录；组成模型以及以同样记录拟合的
# 基线模型使用这些记录。
d$common <- rowSums(d[cfg$companion])>0

# ---- (3) The 15 years with estimates for all six fleets / 六支船队都有估计的 15 个年份 ----
# Each index is divided by its mean over the years present in the model data of all six fleets (Spain's
# eligible records and the records of the composition models for the other fleets): the 15-year mean.
# 每个指数都除以六支船队模型数据（西班牙为合格记录，其他船队为组成模型的记录）都有的年份上的均值，即 15 年均值。
refs <- Reduce(intersect, lapply(cfg$flags, function(f) unique(d$Year[d$Flag==f & (f=="ESP"|d$common)])))
stopifnot(length(refs)>=3, all(refs<=cfg$support_end))
cat("Years of the 15-year mean / 15 年均值的年份:", sort(refs), "\n")

# ---- (4) The 31 fits / 31 个拟合 ----
# One baseline model per fleet fitted to all eligible records, then for each tuna-targeting longline fleet the
# baseline model fitted to the records of the composition models and four composition models (cluster or PC,
# with swordfish excluded from or included in the composition).
# 每支船队一个用全部合格记录拟合的基线模型；每支金枪鱼延绳钓船队另有以组成模型记录拟合的基线模型和四个组成模型
# （聚类或 PC，组成中排除或包含剑鱼）。
flags <- cfg$flags
tasks <- rbind(data.frame(flag=flags, model="all"), expand.grid(flag=flags[-1],
  model=c("common","cluster_exclude","pc_exclude","cluster_include","pc_include"), stringsAsFactors=FALSE))
tasks$id <- paste(tasks$flag, tasks$model, sep="_")
write.csv(tasks, file.path(out, "tasks.csv"), row.names=FALSE)

# ---- (5) Cells with positive swordfish counts (Fig. 2) / 有剑鱼正尾数的网格（图 2） ----
foot <- aggregate(SWOn~Flag+cell+LonC5+LatC5, d[d$SWOn>0,], length)
names(foot)[5] <- "positive_records"
write.csv(foot, file.path(out, "footprint.csv"), row.names=FALSE)

# Save the model input used by the later steps.
# 保存后续步骤使用的模型输入。
saveRDS(list(d=d, grid=grid, coverage=coverage, refs=sort(refs), tasks=tasks), file.path(out, "input.rds"))

# ---- (6) Cells in each assessment area of H2 (Table S1) / H2 各评估分区的网格数（表 S1） ----
# Areas 2 and 3 lie west of 100W and areas 4 and 5 east of it, each split at 20S; a cell belongs to the area of its centre.
# 分区 2、3 位于 100°W 以西，分区 4、5 位于其以东，各以 20°S 分南北；网格按其中心归入分区。
h2 <- jsonlite::fromJSON("config/targeting_public_h2.json")$assessment_areas
s <- read.csv(file.path(out, "support.csv")); f <- read.csv(file.path(out, "footprint.csv"))
assessment_area <- function(lon,lat) {
  stopifnot(!any(lon==h2$split_lon),!any(lat==h2$split_lat))
  ifelse(lat>=10,1L,ifelse(lon<h2$split_lon,ifelse(lat>h2$split_lat,2L,3L),ifelse(lat>h2$split_lat,4L,5L)))
}
tab <- function(set,flag,definition,cells) {
  cells <- unique(cells[c("cell","LonC5","LatC5")])
  n <- table(factor(assessment_area(cells$LonC5,cells$LatC5),levels=1:5))
  stopifnot(n[["1"]]==0)
  data.frame(set,flag,definition,area2=n[["2"]],area3=n[["3"]],area4=n[["4"]],area5=n[["5"]],total=nrow(cells))
}
# Count the core area, the cells with Spanish catches, the cells where only the tuna-targeting longline fleets
# reported swordfish, and each area of overlap ("shared_grounds" in the table).
# 统计核心作业区、有西班牙渔获的网格、只有金枪鱼延绳钓船队报告剑鱼的网格，以及各重叠区域（表中为 "shared_grounds"）。
threshold <- 3L
spf <- subset(f,Flag=="ESP"); oth <- subset(f,Flag!="ESP")
rows <- list(
  tab("core_area","ESP","cells with Spanish records in at least 3 years to 2023 (support.csv, domain ESP, shared_years >= 3)",
      subset(s,domain=="ESP" & shared_years>=threshold)),
  tab("spain_footprint","ESP","cells with at least one positive Spanish swordfish record (footprint.csv; Fig. 2a)",spf),
  tab("added_cells","JPN+KOR+CHN+TWN+VUT","cells with positive records from the five tuna-targeting fleets but none from Spain (Fig. 2b)",
      subset(oth,!cell %in% spf$cell)))
for (fl in c("JPN","KOR","CHN","TWN","VUT"))
  rows[[length(rows)+1L]] <- tab("shared_grounds",fl,
    "cells where Spain and the fleet both report in at least 3 common years to 2023 (support.csv, shared_years >= 3)",
    subset(s,domain==fl & shared_years>=threshold))
res <- do.call(rbind,rows)
stopifnot(all(res$area2+res$area3+res$area4+res$area5==res$total))
print(res[-3])
write.csv(res, file.path(out, "area_structure.csv"), row.names=FALSE)

# ---- Compare with the archived results / 与存档结果比较 ----
source("analysis_walkthrough/common.R")
for (name in c("support.csv","footprint.csv","tasks.csv","area_structure.csv")) compare_table(name)
same <- identical(readRDS(file.path(out,"input.rds")), readRDS("tmp/directed_fishery/input.rds"))
cat("  input.rds:", if (same) "matches the archived model input / 与存档模型输入一致\n" else "DIFFERS from the archived model input / 与存档模型输入不一致\n")
