# =====================================================================================================
# Step 17. Fit the composition models and the baseline models fitted to the same records (Section 2.3, Tables S2 and S6)
# 第 17 步  拟合组成模型及以同样记录拟合的基线模型（第 2.3 节，表 S2 和表 S6）
#
# Fits rows 7-31 of tasks.csv. For each tuna-targeting longline fleet, three models are fitted to the same
# records (those with at least one retained albacore, bigeye tuna, yellowfin tuna, blue marlin or striped
# marlin): the baseline model (f = 0), the PC model with the retained component scores of step 16, and the
# cluster model with the clusters of step 16 (reference cluster with effect zero); the PC and cluster models are
# fitted with swordfish excluded from and included in the composition. Indices are predicted at PC scores of
# zero or for the reference cluster, with the mesh, numerical checks, bias correction and 15-year mean of
# step 07. Then combines the tables of all 31 fits, with the change in AIC from the baseline model fitted to
# the same records (Table S6).
# 拟合 tasks.csv 的第 7-31 行。对每支金枪鱼延绳钓船队，以相同记录（至少留存一尾长鳍金枪鱼、大眼金枪鱼、黄鳍金枪鱼、
# 蓝枪鱼或条纹四鳍旗鱼的记录）拟合三个模型：基线模型（f = 0）、含第 16 步保留主成分得分的 PC 模型，以及含第 16 步
# 聚类的聚类模型（参照聚类的效应为零）；PC 和聚类模型分别在组成中排除和包含剑鱼的情况下拟合。指数在 PC 得分为零或
# 参照聚类处预测，三角网格、数值检查、偏差校正和 15 年均值与第 07 步相同。随后合并全部 31 个拟合的结果表，并给出
# 相对于以同样记录拟合的基线模型的 AIC 变化（表 S6）。
#
# Reads / 读取:  analysis_walkthrough/output/input.rds (step 04), fits_baseline.csv and fits/ (step 07),
#                composition_selection.csv (step 16), config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/fits/<id>_*.csv of the 25 fits, and for all 31 fits fits.csv,
#                indices.csv, model_diagnostics.csv (Tables S2 and S6 / 表 S2 和 S6), model_parameters.csv
#                (Table S2 / 表 S2), residual_quantiles.csv (Fig. S2 / 图 S2)
#
# Run from the project folder; the first line makes R use one core per process and read the Chinese folder
# names. The 25 fits run on --workers N processes (default 4); the fits use their own random seeds, so the
# results do not depend on the number of workers.
# 在项目文件夹中运行；第一行让 R 每个进程只用一个核心，并能读取中文文件夹名。25 个拟合在 --workers N 个进程上运行
# （默认 4 个）；各拟合使用自己的随机种子，因此结果与进程数无关。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/17_fit_composition_models.R --workers 4
# =====================================================================================================

source("analysis_walkthrough/model_functions.R")
out <- "analysis_walkthrough/output"
dir <- file.path(out, "fits")
cat("Step 17: fitting the composition models / 第 17 步：拟合组成模型\n")
z <- readRDS(file.path(out, "input.rds"))
cfg <- jsonlite::fromJSON("config/directed_fishery.json")
argv <- commandArgs(trailingOnly=TRUE)
workers <- if ("--workers" %in% argv) as.integer(argv[match("--workers", argv) + 1L]) else 4L

# Rows 7-31 of tasks.csv: the baseline model fitted to the same records and four composition models per fleet.
# tasks.csv 第 7-31 行：每支船队以同样记录拟合的基线模型和四个组成模型。
rows <- which(z$tasks$model!="all")
stopifnot(identical(rows, 7:31))
fits <- parallel::mclapply(rows, function(i) fit_task(i, z, cfg, dir, diagnostics=TRUE),
                           mc.cores=workers, mc.preschedule=FALSE)
failed <- !vapply(fits, is.data.frame, logical(1))
if (any(failed)) stop("Fits not accepted / 拟合未被接受: ", paste(z$tasks$id[rows][failed], collapse=", "), "\n",
                      paste(unlist(fits[failed]), collapse="\n"))
fits <- do.call(rbind, fits)

# The numbers of retained components and clusters equal those of step 16.
# 保留的主成分数和聚类数与第 16 步相同。
selection <- subset(read.csv(file.path(out, "composition_selection.csv")), selected)
composition <- subset(fits, grepl("cluster|pc", model))
key <- match(paste(sub(".*_", "", composition$model), composition$flag), paste(selection$variant, selection$flag))
stopifnot(!anyNA(key), all(composition$k==selection$k[key]), all(composition$pcs==selection$retained_pcs[key]))

# ---- Tables of all 31 fits / 全部 31 个拟合的结果表 ----
read_parts <- function(suffix) do.call(rbind, lapply(z$tasks$id, function(x) read.csv(file.path(dir, paste0(x, suffix)))))
all_fits <- rbind(read.csv(file.path(out, "fits_baseline.csv")), fits)
all_fits <- all_fits[match(z$tasks$id, all_fits$id),]
stopifnot(all(all_fits$status=="accepted"))
write.csv(all_fits, file.path(out, "fits.csv"), row.names=FALSE)
write.csv(read_parts("_series.csv"), file.path(out, "indices.csv"), row.names=FALSE)
# Change in AIC from the baseline model fitted to the same records, within each fleet (Table S6).
# 各船队内相对于以同样记录拟合的基线模型的 AIC 变化（表 S6）。
diagnostics <- read_parts("_diagnostics.csv")
common_aic <- with(subset(diagnostics, model=="common"), setNames(AIC, flag))
diagnostics$delta_AIC_vs_common <- ifelse(diagnostics$model=="all", NA, diagnostics$AIC - common_aic[diagnostics$flag])
write.csv(diagnostics, file.path(out, "model_diagnostics.csv"), row.names=FALSE)
write.csv(read_parts("_parameters.csv"), file.path(out, "model_parameters.csv"), row.names=FALSE)
write.csv(read_parts("_qq.csv"), file.path(out, "residual_quantiles.csv"), row.names=FALSE)
table_s6 <- reshape(subset(diagnostics, model!="all", c(flag, model, delta_AIC_vs_common)), idvar="flag",
                    timevar="model", direction="wide")
table_s6$baseline_AIC <- common_aic[table_s6$flag]
print(table_s6, digits=6, row.names=FALSE)

# ---- Compare with the archived results / 与存档结果比较 ----
source("analysis_walkthrough/common.R")
compare_table("fits.csv", ignore="rows_sha256")
compare_table("indices.csv")
compare_table("model_diagnostics.csv", ignore="rows_sha256")
compare_table("model_parameters.csv")
compare_table("residual_quantiles.csv")
