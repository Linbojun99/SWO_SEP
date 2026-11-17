# =====================================================================================================
# Step 14. Change from 2019 to 2020 in the same cells and quarter for China and Spain (Section 2.2, Section S9)
# 第 14 步  中国与西班牙在相同网格和季度内由 2019 到 2020 年的变化（第 2.2 节，补充材料 S9）
#
# In the area of overlap between China and Spain, finds the 28 cells with eligible records from both fleets
# in 2019 and 2020 and the 13 of them with records from both fleets in the same calendar quarter of both
# years. For each fleet, refits the baseline model with a separate intercept for every observed year and
# quarter and constant monthly differences within each quarter, keeping the monthly records, offset for the
# number of hooks, mesh, spatial field and independent annual spatiotemporal fields of the baseline model, and
# refits the baseline model itself. From 4,000 joint draws of each fit it calculates the 2020/2019 change in
# log CPUE in each cell and quarter for the same calendar month, the area-weighted China minus Spain difference
# in each quarter, the mean of the four quarters with equal weights, and simultaneous intervals for the six
# comparisons between quarters.
# 在中国与西班牙的重叠区域内，找出 2019 和 2020 年两支船队都有合格记录的 28 个网格，以及其中两年在同一日历季度都有
# 两支船队记录的 13 个网格。对每支船队，以每个有观测的年份与季度组合各一个截距、各季度内月份差异恒定的形式重新拟合
# 基线模型，月度记录、钩数偏置项、三角网格、空间随机场和逐年独立的时空随机场与基线模型相同；并重新拟合基线模型本身。
# 由每个拟合的 4,000 次联合抽样，计算同一日历月份下各网格和季度 2020/2019 年对数 CPUE 的变化、各季度按面积加权的
# 中国减西班牙之差、四个季度等权重的均值，以及六项季度间比较的同时区间。
#
# Reads / 读取:  analysis_walkthrough/output/input.rds (step 04), fits_baseline.csv (step 07),
#                config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/quarterly_matched_cells.csv, quarterly_fit_checks.csv,
#                quarterly_shared_change.csv (Fig. 4b / 图 4b), quarterly_between_seasons.csv
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/14_quarterly_change.R
# =====================================================================================================

suppressPackageStartupMessages(library(sdmTMB))
out <- "analysis_walkthrough/output"
cat("Step 14: change from 2019 to 2020 in the same cells and quarter / 第 14 步：相同网格和季度内 2019 到 2020 年的变化\n")
cfg <- jsonlite::fromJSON("config/directed_fishery.json")
z <- readRDS(file.path(out, "input.rds"))
baseline <- read.csv(file.path(out, "fits_baseline.csv"))
flags <- c("CHN", "ESP")
years <- 2019:2020
n_draw <- 4000L
d <- subset(z$d, Flag %in% flags)
d$quarter <- (d$Month - 1L) %/% 3L + 1L
shared <- subset(z$coverage, domain == "CHN" & shared_years >= 3)$cell
target <- subset(d, Year %in% years & cell %in% shared)

# A cell (or cell and quarter) is compared only if each fleet has eligible records there in both years.
# 只有两支船队在两年都有合格记录的网格（或网格与季度）才参与比较。
common_units <- function(keys) {
  observed <- unique(target[c(keys, "Flag", "Year")])
  hits <- aggregate(rep(1L, nrow(observed)), observed[keys], sum)
  hits[hits$x == 4L, keys, drop = FALSE]
}
annual_cells <- common_units("cell")
matched <- common_units(c("cell", "quarter"))
stopifnot(length(shared) == 49L, nrow(annual_cells) == 28L, length(unique(matched$cell)) == 13L)
cat("Cells in both years / 两年都有记录的网格:", nrow(annual_cells), "| in the same quarter / 同一季度:",
    length(unique(matched$cell)), "| cells per quarter / 各季度网格数:", tabulate(matched$quarter, nbins = 4), "\n")
matched <- merge(matched, z$grid[c("cell", "LonC5", "LatC5", "area")], by = "cell")
annual_cells <- merge(annual_cells, z$grid[c("cell", "area")], by = "cell")

# Weights: within each quarter, cell areas summing to one; the 4-quarter mean gives each quarter one quarter of
# the weight; the check on the choice of cells uses all 28 cells with the baseline model.
# 权重：各季度内为网格面积，和为 1；四季度均值给每个季度四分之一的权重；对网格选择的核对用基线模型和全部 28 个网格。
units <- unique(rbind(matched[c("cell", "quarter")], data.frame(cell = annual_cells$cell, quarter = 1L)))
units <- merge(units, z$grid, by = "cell")
unit_key <- paste(units$cell, units$quarter)
weight_rows <- function(cells, quarter) {
  pos <- match(paste(cells$cell, quarter), unit_key)
  stopifnot(!anyNA(pos))
  w <- numeric(nrow(units))
  w[pos] <- cells$area / sum(cells$area)
  w
}
weights <- t(sapply(1:4, function(q) weight_rows(subset(matched, quarter == q), q)))
rownames(weights) <- paste0("Q", 1:4)
weights <- rbind(weights, equal_quarters = colMeans(weights), annual_28 = weight_rows(annual_cells, 1))
stopifnot(max(abs(rowSums(weights) - 1)) < 1e-12)

# Fit one fleet with year and month effects (season = FALSE) or year-quarter intercepts and within-quarter month
# indicators (season = TRUE), and return the weighted 2020/2019 changes in log CPUE for every joint draw.
# 以年份和月份效应（season = FALSE）或年份与季度截距加季内月份指示变量（season = TRUE）拟合一支船队，并返回每次联合
# 抽样中按权重计算的 2020/2019 年对数 CPUE 变化。
fit_changes <- function(flag, season) {
  fleet <- subset(d, Flag == flag)
  old <- baseline[baseline$id == paste0(flag, "_all"), ]
  stopifnot(nrow(old) == 1, old$status == "accepted", old$records == nrow(fleet))
  series_years <- sort(unique(fleet$Year))
  levels_yq <- sort(unique(4L * fleet$Year + fleet$quarter))
  # January, April, July and October are the reference months of the four quarters.
  # 1、4、7 和 10 月为四个季度的参照月份。
  month_terms <- paste0("m", c(2, 3, 5, 6, 8, 9, 11, 12))
  design <- function(x) {
    x$t <- match(x$Year, series_years)
    x$year_f <- factor(x$Year, levels = series_years)
    x$month_f <- factor(x$Month, levels = sort(unique(fleet$Month)))
    x$year_quarter <- factor(4L * x$Year + x$quarter, levels = levels_yq)
    for (term in month_terms) x[[term]] <- as.integer(x$Month == as.integer(sub("m", "", term)))
    x
  }
  fleet <- design(fleet)
  form <- reformulate(if (season) c("year_quarter", month_terms) else c("year_f", "month_f"), "SWOn")
  X <- model.matrix(form, fleet)
  stopifnot(qr(X)$rank == ncol(X))
  loc <- subset(z$grid, cell %in% fleet$cell)
  geo <- fmesher::fm_mesh_2d(loc = as.matrix(loc[c("x","y")]), cutoff = 100,
                              max.edge = c(850, 1800), offset = c(650, 1400))
  mesh <- make_mesh(fleet, c("x", "y"), mesh = geo)
  cat(" ", flag, deparse(form), "\n")
  fit <- sdmTMB(form, data = fleet, mesh = mesh, time = "t", family = nbinom2(),
                offset = log(fleet$Hooks / 1000), spatial = "on", spatiotemporal = "iid",
                share_range = TRUE, control = sdmTMBcontrol(backend = "tmb", newton_loops = 2,
                                                           get_joint_precision = TRUE))
  gradient <- max(abs(fit$tmb_obj$gr(fit$model$par)))
  stopifnot(fit$model$convergence == 0L, isTRUE(fit$sd_report$pdHess), gradient < cfg$gradient_limit)
  # Predict each cell and quarter in 2019 and 2020 for the first month of the quarter at 1,000 hooks, and draw
  # 4,000 joint realizations of the fixed effects and random fields (seed cfg$seed + 100 season + fleet number).
  # 以 1,000 钩、季度第一个月预测 2019 和 2020 年的各网格与季度，并抽取固定效应和随机场的 4,000 次联合实现（种子
  # cfg$seed + 100 × season + 船队编号）。
  p <- units[rep(seq_len(nrow(units)), times = 2), ]
  p$Year <- rep(years, each = nrow(units))
  p$Month <- 3L * p$quarter - 2L
  p <- design(p)
  estimate <- predict(fit, newdata = p, type = "link", offset = rep(0, nrow(p)))$est
  seed <- cfg$seed + 100L * as.integer(season) + match(flag, flags)
  set.seed(seed)
  draws <- predict(fit, newdata = p, type = "link", offset = rep(0, nrow(p)), nsim = n_draw)
  values <- cbind(estimate, draws)
  stopifnot(nrow(values) == nrow(p), ncol(values) == n_draw + 1L, all(is.finite(values)))
  delta <- values[p$Year == 2020, , drop = FALSE] - values[p$Year == 2019, , drop = FALSE]
  w <- if (season) weights[rownames(weights) != "annual_28", , drop = FALSE] else weights
  metrics <- data.frame(flag, model = if (season) "year_quarter" else "annual",
                        records = nrow(fleet), years = length(series_years), fixed_effect_columns = ncol(X),
                        optimizer = fit$model$convergence,
                        pd_hessian = fit$sd_report$pdHess, max_abs_gradient = gradient,
                        AIC = AIC(fit), seed)
  list(changes = w %*% delta, metrics = metrics)
}

support_sizes <- c(setNames(tabulate(matched$quarter, nbins = 4), paste0("Q", 1:4)),
                   equal_quarters = length(unique(matched$cell)), annual_28 = nrow(annual_cells))
# Estimate and pointwise 95% interval on the ratio scale from the 2.5th and 97.5th percentiles of the draws.
# 由抽样的第 2.5 和 97.5 百分位数得到比值尺度上的估计和逐点 95% 区间。
contrast_ci <- function(x) {
  stopifnot(all(is.finite(x)))
  q <- quantile(x[-1], c(.025, .975), names = FALSE)
  setNames(exp(c(x[1], q)), c("estimate", "lower", "upper"))
}
fits <- list()
tables <- list()
for (m in c("annual", "year_quarter")) {
  f <- setNames(lapply(flags, fit_changes, season = m == "year_quarter"), flags)
  fits[[m]] <- do.call(rbind, lapply(f, "[[", "metrics"))
  # China minus Spain, within each joint draw; the draws of the two fleets are independent.
  # 在每次联合抽样内计算中国减西班牙；两支船队的抽样相互独立。
  contrast <- f$CHN$changes - f$ESP$changes
  periods <- rownames(contrast)
  for (name in c(flags, "CHN_over_ESP")) {
    vals <- if (name == "CHN_over_ESP") contrast else f[[name]]$changes
    ci <- t(apply(vals, 1, contrast_ci))
    tables[[paste(m, name)]] <- data.frame(model = m, series = name, period = periods,
                                           n_cells = support_sizes[periods], ci, row.names = NULL)
  }
}

# Simultaneous intervals for the six comparisons between quarters: the 95th percentile of the maximum absolute
# standardized deviation across the six differences.
# 六项季度间比较的同时区间：六个差值的最大绝对标准化偏差的第 95 百分位数。
pairs <- combn(1:4, 2)
diff <- contrast[pairs[1, ], ] - contrast[pairs[2, ], ]
sds <- apply(diff[, -1, drop = FALSE], 1, sd)
critical <- unname(quantile(apply(abs(sweep(diff[, -1] - diff[, 1], 1, sds, "/")), 2, max), .95))
between <- data.frame(numerator = paste0("Q", pairs[1, ]), denominator = paste0("Q", pairs[2, ]),
                      ratio = exp(diff[, 1]), lower_simultaneous = exp(diff[, 1] - critical * sds),
                      upper_simultaneous = exp(diff[, 1] + critical * sds))
outputs <- list(quarterly_matched_cells = matched, quarterly_fit_checks = do.call(rbind, fits),
                quarterly_shared_change = do.call(rbind, tables), quarterly_between_seasons = between)
for (name in names(outputs)) write.csv(outputs[[name]], file.path(out, paste0(name, ".csv")), row.names = FALSE)
print(subset(outputs$quarterly_shared_change, series == "CHN_over_ESP"), digits = 3, row.names = FALSE)

# ---- Compare with the archived results / 与存档结果比较 ----
source("analysis_walkthrough/common.R")
for (name in names(outputs)) compare_table(paste0(name, ".csv"), ignore = "source_rows_sha256")
