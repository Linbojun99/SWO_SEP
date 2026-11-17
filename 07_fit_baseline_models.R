# =====================================================================================================
# Step 07. Fit the baseline model of each fleet to all eligible records (Section 2.2)
# 第 07 步  以全部合格记录拟合各船队的基线模型（第 2.2 节）
#
# Fits the first six rows of tasks.csv: for each fleet, the spatiotemporal NB2 model of the retained swordfish
# count with year and month effects, the offset for the number of hooks, a spatial field and independent
# annual spatiotemporal fields (Eqs. 1 and 2, f = 0). Each accepted fit gives the bias-corrected annual
# indices over Spain's core area or the area of overlap (thresholds of two, three and four years), scaled to
# the 15-year mean. The same fits give the parameter estimates, randomized quantile residuals and proportions
# of zeros (Section S3), the residuals by cell and year (Section S4) and the relative CPUE in each cell of the
# areas of overlap (Section S9).
# 拟合 tasks.csv 的前六行：对每支船队，拟合剑鱼留存尾数的时空 NB2 模型，含年份效应、月份效应、钩数偏置项、空间
# 随机场和逐年独立的时空随机场（式 1 和式 2，f = 0）。每个被接受的拟合给出西班牙核心作业区或重叠区域上（两年、
# 三年和四年阈值）经偏差校正的年度指数，并除以 15 年均值。同一批拟合还给出参数估计值、随机分位数残差和零值比例
# （补充材料 S3）、按网格和年份汇总的残差（补充材料 S4），以及重叠区域各网格的相对 CPUE（补充材料 S9）。
#
# Reads / 读取:  analysis_walkthrough/output/input.rds (step 04), config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/fits_baseline.csv, indices_baseline.csv and, for each fleet,
#                fits/<id>_fit.csv, _series.csv, _diagnostics.csv, _parameters.csv, _qq.csv,
#                fits/<flag>_residual_cells.csv, fits/<flag>_cpue_map.csv
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/07_fit_baseline_models.R
# =====================================================================================================

source("analysis_walkthrough/model_functions.R")
out <- "analysis_walkthrough/output"
dir <- file.path(out, "fits")
dir.create(dir, recursive=TRUE, showWarnings=FALSE)
cat("Step 07: fitting the six baseline models to all eligible records\n")
cat("第 07 步：以全部合格记录拟合六个基线模型\n")
z <- readRDS(file.path(out, "input.rds"))
cfg <- jsonlite::fromJSON("config/directed_fishery.json")

# The six baseline models fitted to all eligible records are rows 1-6 of tasks.csv.
# 用全部合格记录拟合的六个基线模型是 tasks.csv 的第 1-6 行。
rows <- which(z$tasks$model=="all")
stopifnot(identical(rows, 1:6))
fits <- do.call(rbind, lapply(rows, function(i) {
  cat(sprintf("  fitting / 正在拟合 %s ...\n", z$tasks$id[i]))
  fit_task(i, z, cfg, dir, diagnostics=TRUE, spatial=TRUE)
}))
print(fits[c("id","records","cells","optimizer","pdHess","gradient","status")])

# Combine the fit summaries and the annual indices of the six fits.
# 合并六个拟合的汇总和年度指数。
read_parts <- function(suffix) do.call(rbind, lapply(z$tasks$id[rows], function(x) read.csv(file.path(dir, paste0(x, suffix)))))
series <- read_parts("_series.csv")
stopifnot(!anyDuplicated(series[c("id","domain","threshold","region","year")]))
write.csv(fits, file.path(out, "fits_baseline.csv"), row.names=FALSE)
write.csv(series, file.path(out, "indices_baseline.csv"), row.names=FALSE)

# ---- Compare with the archived results / 与存档结果比较 ----
# The archived tables hold all 31 fits; these six are compared with their rows.
# 存档表包含全部 31 个拟合；这六个拟合与其中相应的行比较。
source("analysis_walkthrough/common.R")
ids <- z$tasks$id[rows]
rows_of <- function(name) subset(read.csv(file.path(archive, name)), id %in% ids)
compare_table("fits_baseline.csv", old=rows_of("fits.csv"), ignore="rows_sha256")
compare_table("indices_baseline.csv", old=rows_of("indices.csv"))
compare_table("fits/*_diagnostics.csv", new=read_parts("_diagnostics.csv"), old=rows_of("model_diagnostics.csv"), ignore=c("rows_sha256","delta_AIC_vs_common"))
compare_table("fits/*_parameters.csv", new=read_parts("_parameters.csv"), old=rows_of("model_parameters.csv"))
compare_table("fits/*_qq.csv", new=read_parts("_qq.csv"), old=rows_of("residual_quantiles.csv"))
