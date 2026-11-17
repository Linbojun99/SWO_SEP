# =====================================================================================================
# Step 03. Prepare the public domain records (Section 2.1)
# 第 03 步  整理公开记录（第 2.1 节）
#
# Keeps the six fleets from 2005 to 2023, removes every record with positive catch weight for any species and
# every record without hooks, marks the fleet-years with at least one positive swordfish count anywhere in the
# eastern Pacific, keeps the H2 study area (south of 10N, east of 150W) and adds the projected position and
# the spherical area of each 5-degree cell. Each record is one fleet in one 5-degree cell in one calendar month.
# 保留 2005-2023 年的六支船队，删除任一鱼种有正渔获重量的记录和没有钩数的记录，标记在东太平洋任一位置至少有一个
# 剑鱼正尾数的船队年份，保留 H2 研究海域（10°N 以南、150°W 以东），并为每个 5° 网格加上投影坐标和球面面积。每条
# 记录代表一支船队在某个日历月份、某个 5° 网格的作业。
#
# Reads / 读取:  config/directed_fishery.json, data/raw/longline/PublicLLTunaBillfishNum.csv
# Writes / 写出: analysis_walkthrough/output/records.rds (d = H2 records before the selection of fleet-years,
#                grid = the 5-degree cells / d = 船队年份筛选前的 H2 记录，grid = 5° 网格)
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/03_prepare_records.R
# =====================================================================================================

out <- "analysis_walkthrough/output"
dir.create(out, recursive=TRUE, showWarnings=FALSE)

# Read the analysis settings: fleets, years, area limits, map projection and random seeds.
# 读取分析设置：船队、年份、区域边界、地图投影和随机种子。
cfg <- jsonlite::fromJSON("config/directed_fishery.json")
cat("Step 03: preparing the public domain records from", cfg$num, "\n")
cat("第 03 步：由", cfg$num, "整理公开记录\n")

# Read the raw file, which is never changed, and number its rows.
# 读取原始文件（从不修改），并为各行编号。
raw <- read.csv(cfg$num)
raw$source_row <- seq_len(nrow(raw))
cat("Raw rows / 原始行数:", nrow(raw), "\n")

# Keep Spain, Japan, Korea, China, Chinese Taipei and Vanuatu from 2005 to 2023.
# 保留 2005-2023 年的西班牙、日本、韩国、中国、中华台北和瓦努阿图。
raw <- subset(raw, Flag %in% cfg$flags & Year >= min(cfg$years) & Year <= max(cfg$years))

# Each fleet, year, month and cell appears once, with no missing hooks or counts.
# 每个船队、年、月和网格只出现一次，钩数和尾数都没有缺失。
stopifnot(!anyDuplicated(raw[c("Flag","Year","Month","LonC5","LatC5")]),
          !anyNA(raw[c("Hooks","SWOn",cfg$companion)]))

# Remove records with positive catch weight for any species (columns ending in "mt") and records without
# hooks, so that numerical and weight catches are never combined and every record gives CPUE in numbers.
# 删除任一鱼种有正渔获重量（列名以 "mt" 结尾）的记录和没有钩数的记录，使尾数渔获与重量渔获从不合并，且每条记录
# 都能给出以尾数计的 CPUE。
w <- raw[grepl("mt$", names(raw))]
stopifnot(!anyNA(w))
raw <- raw[rowSums(w>0)==0 & raw$Hooks>0,]

# The public domain totals do not distinguish all non-reporting from zero retained catch, so a fleet's records
# in a year are used only if at least one of them has a positive swordfish count.
# 公开汇总资料不能区分所有未报告与留存渔获为零的情形，因此只有当某船队某年至少一条记录的剑鱼尾数为正时，才使用
# 该船队该年的记录。
raw$eligible <- ave(raw$SWOn, raw$Flag, raw$Year, FUN=max)>0

# Keep the H2 study area: cell centres south of 10N and east of 150W.
# 保留 H2 研究海域：网格中心位于 10°N 以南、150°W 以东。
d <- subset(raw, LatC5 < cfg$north_boundary & LonC5 > cfg$west_boundary)
d$cell <- paste(d$LonC5, d$LatC5, sep="_")

# List the 5-degree cells, sorted by longitude and then latitude.
# 列出 5° 网格，先按经度、再按纬度排序。
grid <- unique(d[c("cell","LonC5","LatC5")])
grid <- grid[order(grid$LonC5, grid$LatC5),]

# Project the cell centres to an azimuthal equidistant projection in kilometres, centred at 15S, 110W.
# 把网格中心投影到以 15°S、110°W 为中心的方位等距投影，单位为千米。
xy <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(grid, coords=c("LonC5","LatC5"), crs=4326), cfg$projection))
grid$x <- xy[,1]; grid$y <- xy[,2]

# Spherical area of each full 5-degree cell in square kilometres, later used as the index weight.
# 每个完整 5° 网格的球面面积（平方千米），之后用作指数权重。
grid$area <- (sin((grid$LatC5+2.5)*pi/180)-sin((grid$LatC5-2.5)*pi/180))*5*pi/180*6371.0088^2

# Add the projected position to every record.
# 把投影坐标加到每条记录上。
d$x <- grid$x[match(d$cell, grid$cell)]; d$y <- grid$y[match(d$cell, grid$cell)]
cat("H2 records / H2 海域记录数:", nrow(d), "| cells / 网格数:", nrow(grid), "\n")

# Save the records and cells for step 04.
# 保存记录与网格，供第 04 步使用。
saveRDS(list(d=d, grid=grid), file.path(out, "records.rds"))

# ---- Compare with the archived results / 与存档结果比较 ----
# The archived model input holds the same cells and, after the selection of fleet-years in step 04, the same records.
# 存档的模型输入包含相同的网格，以及经第 04 步船队年份筛选后相同的记录。
archived <- readRDS("tmp/directed_fishery/input.rds")
same_grid <- identical(grid, archived$grid)
same_records <- identical(subset(d, eligible), archived$d[setdiff(names(archived$d), "common")])
cat("  cells / 网格:", if (same_grid) "match the archived model input / 与存档模型输入一致" else "DIFFER from the archived model input / 与存档模型输入不一致", "\n")
cat("  eligible records / 合格记录:", if (same_records) "match the archived model input / 与存档模型输入一致" else "DIFFER from the archived model input / 与存档模型输入不一致", "\n")
