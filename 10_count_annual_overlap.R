# =====================================================================================================
# Step 10. Annual overlap of the eligible records of each fleet and Spain (Section 2.2, Section S1)
# 第 10 步  各船队与西班牙合格记录的逐年重叠（第 2.2 节，补充材料 S1）
#
# For each tuna-targeting longline fleet and year, counts the cells of the area of overlap with eligible
# records from Spain, from the fleet and from both, and the combinations of cell and month with records from
# both. For China and Spain, recalculates the rank correlation of the two baseline indices of Fig. 3 over all
# years with estimates for both fleets, over the years with at least one cell with records from both fleets,
# and over the years with at least one cell and month with records from both. No model is fitted.
# 对每支金枪鱼延绳钓船队和每个年份，统计重叠区域内有西班牙合格记录、有该船队合格记录以及两者都有记录的网格数，
# 以及两者都有记录的网格与月份组合数。对中国与西班牙，分别在两支船队都有估计的全部年份、至少一个网格两者都有记录
# 的年份、至少一个网格与月份组合两者都有记录的年份上，重新计算图 3 两条基线指数的秩相关。不拟合模型。
#
# Reads / 读取:  analysis_walkthrough/output/input.rds (step 04), indices_baseline.csv (step 07),
#                config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/pair_overlap.csv, china_spain_overlap_correlations.csv
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/10_count_annual_overlap.R
# =====================================================================================================

out <- "analysis_walkthrough/output"
cat("Step 10: annual overlap of the eligible records / 第 10 步：合格记录的逐年重叠\n")
z <- readRDS(file.path(out, "input.rds"))
cfg <- jsonlite::fromJSON("config/directed_fishery.json")
ind <- read.csv(file.path(out, "indices_baseline.csv"))
years <- seq.int(cfg$years[1], cfg$years[2])
threshold <- 3L

# For each fleet and year, count the cells and cell-months of the area of overlap with records from Spain, the fleet or both.
# 对每支船队和每个年份，统计重叠区域内有西班牙、该船队或两者记录的网格数及网格与月份组合数。
others <- setdiff(cfg$flags, "ESP")
rows <- list(); paired <- list()
for (f in others) {
  cells <- z$coverage$cell[z$coverage$domain == f & z$coverage$shared_years >= threshold]
  es <- z$d[z$d$Flag == "ESP" & z$d$cell %in% cells, c("Year","Month","cell")]
  fl <- z$d[z$d$Flag == f & z$d$cell %in% cells, c("Year","Month","cell")]
  for (y in years) {
    a <- unique(es$cell[es$Year == y]); b <- unique(fl$cell[fl$Year == y])
    am <- unique(paste(es$cell, es$Month)[es$Year == y]); bm <- unique(paste(fl$cell, fl$Month)[fl$Year == y])
    rows[[length(rows) + 1L]] <- data.frame(fleet=f, year=y, spain_cells=length(a), fleet_cells=length(b),
      both_cells=length(intersect(a, b)), both_cell_months=length(intersect(am, bm)),
      spain_records=sum(es$Year == y), fleet_records=sum(fl$Year == y))
  }
  # Years with estimates for both fleets: years with eligible records of each fleet anywhere in H2.
  # 两支船队都有估计的年份：两支船队在 H2 任一位置都有合格记录的年份。
  paired[[f]] <- intersect(sort(unique(z$d$Year[z$d$Flag == "ESP"])), sort(unique(z$d$Year[z$d$Flag == f])))
  # The record counts equal the records of the indices of step 07 in the same area of overlap.
  # 记录数与第 07 步同一重叠区域内指数的记录数相同。
  r <- do.call(rbind, rows); r <- r[r$fleet == f, ]
  for (g in c("ESP", f)) {
    s <- ind[ind$flag == g & ind$domain == f & ind$threshold == threshold & ind$region == "all", ]
    col <- if (g == "ESP") "spain_records" else "fleet_records"
    ccol <- if (g == "ESP") "spain_cells" else "fleet_cells"
    m <- match(s$year, r$year)
    stopifnot(nrow(s) > 0, !anyNA(m), all(r[[col]][m] == s$records), all(r[[ccol]][m] == s$observed_cells),
              all(r[[col]][!(r$year %in% s$year)] == 0))
  }
}
res <- do.call(rbind, rows)
stopifnot(all(res$both_cells <= pmin(res$spain_cells, res$fleet_cells)))
write.csv(res, file.path(out, "pair_overlap.csv"), row.names=FALSE)
print(subset(res, fleet == "CHN" & year %in% c(2008, 2015, 2020)), row.names=FALSE)

# Recalculate the China-Spain rank correlation over the years that meet each criterion, keeping the annual
# estimates and the 49 cells of the area of overlap unchanged.
# 在满足各项条件的年份上重新计算中国与西班牙的秩相关，年度估计值和重叠区域的 49 个网格保持不变。
base <- subset(ind, domain == "CHN" & region == "all" & model == "all" & threshold == 3L)
joint <- merge(subset(base, flag == "ESP", c(year, index)), subset(base, flag == "CHN", c(year, index)),
               by="year", suffixes=c("_spain", "_china"))
joint <- merge(joint, subset(res, fleet == "CHN"), by="year")
stopifnot(identical(joint$year, paired[["CHN"]]), all(base$domain_cells == 49))
filters <- list(all_paired=rep(TRUE, nrow(joint)), same_cell=joint$both_cells > 0,
                same_cell_month=joint$both_cell_months > 0)
correlations <- do.call(rbind, lapply(names(filters), function(rule) {
  d <- joint[filters[[rule]], ]
  data.frame(year_filter=rule, n_years=nrow(d), years=paste(d$year, collapse=" "),
             spearman_rho=cor(d$index_spain, d$index_china, method="spearman"))
}))
print(correlations[c("year_filter","n_years","spearman_rho")], row.names=FALSE)
write.csv(correlations, file.path(out, "china_spain_overlap_correlations.csv"), row.names=FALSE)

# ---- Compare with the archived results / 与存档结果比较 ----
source("analysis_walkthrough/common.R")
compare_table("pair_overlap.csv")
compare_table("china_spain_overlap_correlations.csv")
