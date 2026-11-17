# =====================================================================================================
# Step 19. Annual indices with and without composition covariates (Section 2.3, Fig. 5, Fig. S9)
# 第 19 步  含与不含组成协变量的年度指数（第 2.3 节，图 5、图 S9）
#
# Fig. 5 shows, for each tuna-targeting longline fleet in its area of overlap with Spain, the baseline, cluster
# and PC indices fitted to the same records, with swordfish excluded from the composition, and their pointwise
# 95% confidence intervals. Fig. S9 shows the same with swordfish included in the composition.
# 图 5 显示各金枪鱼延绳钓船队在其与西班牙重叠区域内以同样记录拟合的基线、聚类和 PC 指数（组成中排除剑鱼）及其逐点
# 95% 置信区间。图 S9 显示组成中包含剑鱼时的相应指数。
# Same drawing as the Fig. 5 and Fig. S9 part of R/plot_directed_fishery.R.
# 与 R/plot_directed_fishery.R 中图 5 和图 S9 部分的绘图相同。
#
# Reads / 读取:  analysis_walkthrough/output/indices.csv (step 17), config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/Fig5_composition.pdf/.png, FigS9_composition_inclusive.pdf/.png
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/19_plot_composition_indices.R
# =====================================================================================================

source("analysis_walkthrough/plot_functions.R")
cat("Step 19: drawing Fig. 5 and Fig. S9 / 第 19 步：绘制图 5 和图 S9\n")
idx <- read.csv(file.path(out,"indices.csv"))
idx$annual_only <- as.logical(idx$annual_only)
primary <- subset(idx,threshold==3 & !annual_only)
panelnames <- index_panels(subset(primary,model=="all" & region=="all"))
clabel <- c(common="Baseline",cluster_exclude="Cluster",pc_exclude="PC",cluster_include="Cluster",pc_include="PC")
# Five fleets fill five of six panels; the legend sits in the empty sixth panel.
# 五支船队占六个面板中的五个；图例放在空出的第六个面板。
legend_in_slot <- theme(legend.position="inside",legend.position.inside=c(.75,.165),
                        legend.justification.inside=c(.5,.5),legend.direction="vertical")

# Fig. 5: swordfish excluded from the composition. Fig. S9: swordfish included.
# 图 5：组成中排除剑鱼。图 S9：组成中包含剑鱼。
for (variant in c("exclude","include")) {
  models <- c("common",paste0(c("cluster_","pc_"),variant))
  composition <- subset(primary,region=="all" & model %in% models)
  composition$series <- unname(clabel[composition$model])
  composition$panel <- factor(panelnames[match(composition$flag,cfg$flags)],levels=panelnames[-1])
  save_fig(plot_lines(composition,composition_style) + legend_in_slot,
           if (variant=="exclude") "Fig5_composition" else "FigS9_composition_inclusive",152)
}

# ---- Compare with the archived figures / 与存档图比较 ----
source("analysis_walkthrough/common.R")
for (name in c("Fig5_composition.png","FigS9_composition_inclusive.png")) compare_figure(name)
