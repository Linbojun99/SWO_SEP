# =====================================================================================================
# Step 06. Cells with positive swordfish counts by period (Section 2.1, Section 3.1, Table S1 notes)
# 第 06 步  各时段有剑鱼正尾数的网格（第 2.1 节、第 3.1 节、表 S1 注释）
#
# Counts, from the public file, the cells with positive retained swordfish counts reported by Spain, by all six
# fleets and only by the five tuna-targeting longline fleets, for 2005 to 2023 and 2019 to 2023, and lists the
# years with positive counts in each cell where only the tuna-targeting longline fleets reported swordfish. A
# record with a positive swordfish count belongs to an eligible fleet-year by definition.
# 由公开文件统计 2005-2023 年和 2019-2023 年西班牙、全部六支船队以及只有五支金枪鱼延绳钓船队报告剑鱼留存正尾数的
# 网格数，并列出只有金枪鱼延绳钓船队报告剑鱼的每个网格中有正尾数的年份。剑鱼尾数为正的记录按定义属于合格船队年份。
#
# Reads / 读取:  config/directed_fishery.json, data/raw/longline/PublicLLTunaBillfishNum.csv
# Writes / 写出: analysis_walkthrough/output/period_coverage.csv, added_cell_positive_years.csv
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/06_count_catch_cells.R
# =====================================================================================================

source("analysis_walkthrough/common.R")
cat("Step 06: cells with swordfish catches by period / 第 06 步：各时段有剑鱼渔获的网格\n")
cfg <- jsonlite::fromJSON("config/directed_fishery.json")

# Records in H2 with positive hooks, a positive swordfish count and no catch in weight; a missing value never
# counts as positive.
# H2 内钩数为正、剑鱼尾数为正且无重量渔获的记录；缺失值从不算作正值。
raw <- read.csv(cfg$num)
weight <- raw[grepl("mt$", names(raw))]
fish <- raw[which(raw$Flag %in% cfg$flags & raw$Hooks > 0 & raw$SWOn > 0 & raw$LatC5 < cfg$north_boundary &
                  raw$LonC5 > cfg$west_boundary & rowSums(weight > 0, na.rm=TRUE) == 0), ]

# Count the cells for 2005-2023 and for the last five years, 2019-2023; "additional_west" lies west of 100W.
# 统计 2005-2023 年和最后五年 2019-2023 年的网格数；"additional_west" 指 100°W 以西的网格。
cells <- function(x) unique(x[c("LonC5", "LatC5")])
footprints <- list()
for (period_years in list(cfg$years, c(cfg$years[2] - 4, cfg$years[2]))) {
  period <- fish[fish$Year >= period_years[1] & fish$Year <= period_years[2], ]
  spanish <- cells(period[period$Flag == "ESP", ])
  other <- cells(period[period$Flag != "ESP", ])
  added <- other[!row_key(other, names(other)) %in% row_key(spanish, names(spanish)), ]
  combined <- unique(rbind(spanish, other))
  footprints[[length(footprints) + 1L]] <- data.frame(start=period_years[1], end=period_years[2],
    spain=nrow(spanish), combined=nrow(combined), additional=nrow(added), additional_west=sum(added$LonC5 < -100))
  if (identical(as.numeric(period_years), as.numeric(cfg$years))) {
    # A year counts once even if several fleets or months reported a positive catch.
    # 即使有多支船队或多个月份报告正渔获，一个年份也只计一次。
    cell_years <- unique(period[period$Flag != "ESP", c("LonC5", "LatC5", "Year")])
    added_years <- cell_years[row_key(cell_years, c("LonC5", "LatC5")) %in% row_key(added, c("LonC5", "LatC5")), ]
    write.csv(sort_rows(added_years, c("LonC5", "LatC5", "Year")), file.path(out, "added_cell_positive_years.csv"),
              row.names=FALSE)
  }
}
table_cells <- do.call(rbind, footprints)
write.csv(table_cells, file.path(out, "period_coverage.csv"), row.names=FALSE)
print(table_cells, row.names=FALSE)
years <- table(row_key(read.csv(file.path(out, "added_cell_positive_years.csv")), c("LonC5", "LatC5")))
cat(sprintf("Added cells with positive counts in at least three years / 至少三年有正尾数的新增网格: %d; in all 19 years / 全部 19 年: %d\n",
            sum(years >= 3), sum(years == 19)))

# ---- Compare with the archived results / 与存档结果比较 ----
compare_table("period_coverage.csv")
compare_table("added_cell_positive_years.csv")
