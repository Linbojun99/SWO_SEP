# =====================================================================================================
# Step 02. Study design (Fig. 1)
# 第 02 步  研究设计（图 1）
#
# Draws the study design: the reported catch and effort of Spain and the five tuna-targeting longline fleets
# lead to three empirical comparisons (Sections 2.2 and 2.3), and the reported cells and numbers of hooks are
# used to simulate swordfish catches and to compare regional indices under alternative sampling designs
# with the simulated abundance (Section 2.4). The diagram uses no data.
# 绘制研究设计：西班牙和五支金枪鱼延绳钓船队报告的渔获与努力量用于三项实证比较（第 2.2 和 2.3 节）；报告的
# 网格和钩数用于模拟剑鱼渔获，并将不同取样设计下的区域指数与模拟丰度比较（第 2.4 节）。示意图不使用数据。
# Same drawing as R/plot_study_design.R. / 与 R/plot_study_design.R 的绘图相同。
#
# Reads / 读取:  nothing / 无
# Writes / 写出: analysis_walkthrough/output/Fig1_study_design.pdf, Fig1_study_design.png
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/02_draw_study_design.R
# =====================================================================================================

library(grid)
out <- "analysis_walkthrough/output"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
cat("Step 02: drawing the study design (Fig. 1) / 第 02 步：绘制研究设计（图 1）\n")
width_mm <- 190
height_mm <- 96

# The diagram separates the data (left) from the four comparisons (right), after Winker et al. (2018, Fig. 1)
# and Silva et al. (2023, Fig. 1).
# 示意图把数据（左）与四项比较（右）分开，参照 Winker et al.（2018，图 1）和 Silva et al.（2023，图 1）。
draw_design <- function() {
  grid.newpage()
  pushViewport(viewport(xscale=c(0, width_mm), yscale=c(0, height_mm)))
  label <- function(x, y, value, size, bold=FALSE, align="left") {
    grid.text(value, x=x, y=y, default.units="native", just=c(align,"centre"),
              gp=gpar(fontfamily="Arial", fontsize=size, fontface=if (bold) "bold" else "plain",
                      col="#202020", lineheight=1.1))
  }
  link <- function(x, y, tip=FALSE) {
    grid.lines(x=x, y=y, default.units="native",
               arrow=if (tip) arrow(length=unit(1.7,"mm"), type="closed") else NULL,
               gp=gpar(col="#404040", fill="#404040", lwd=1.25))
  }
  comparison <- function(y, title, method) {
    grid.rect(x=66, y=y, width=122, height=20.5, just=c("left","centre"), default.units="native",
              gp=gpar(col="#363636", fill="white", lwd=1.25))
    label(69.5,y+3.7,title,10.6,TRUE)
    label(69.5,y-3.0,method,9.4)
  }

  # Reported catch and effort lead to the three empirical comparisons.
  # 报告的渔获与努力量通向三项实证比较。
  label(28.5,66.5,"Reported catch and effort",10.5,TRUE,"centre")
  label(28.5,58.5,"Spain and five tuna-targeting fleets\nMonthly records in 5° cells",9.3,align="centre")
  link(c(56,60),c(62,62))
  link(c(60,60),c(39,85))
  for (y in c(85,62,39)) link(c(60,65.4),c(y,y),TRUE)
  comparison(85,"Annual indices in the area of overlap",
             "Compare Spain’s annual index with each tuna-targeting fleet’s index.")
  comparison(62,"CPUE in cells fished by both fleets",
             "Compare relative CPUE in cells with records from both fleets.")
  comparison(39,"Sensitivity to catch composition",
             "Compare indices with and without catch composition covariates.")

  # The reported cells and numbers of hooks lead to the simulated catches and the comparison of sampling designs.
  # 报告的网格和钩数通向模拟渔获和取样设计的比较。
  link(c(28.5,28.5),c(51,22),TRUE)
  label(28.5,16.5,"Simulated swordfish catches",10.5,TRUE,"centre")
  label(28.5,8.5,"Reported cells and numbers of hooks\nSpecified abundance and catchability",9.0,align="centre")
  link(c(56,65.4),c(12,12),TRUE)
  comparison(12,"Regional index under alternative sampling designs",
             "Compare each regional index with the simulated abundance trend.")
  popViewport()
}

# Write the PDF and a 600 dpi PNG.
# 写出 PDF 和 600 dpi 的 PNG。
for (ext in c("pdf","png")) {
  file <- file.path(out,paste0("Fig1_study_design.",ext))
  if (ext == "pdf") {
    grDevices::cairo_pdf(file,width=width_mm/25.4,height=height_mm/25.4)
  } else {
    grDevices::png(file,width=width_mm,height=height_mm,units="mm",res=600,type="cairo")
  }
  draw_design()
  invisible(dev.off())
}

# ---- Compare with the archived figure / 与存档图比较 ----
source("analysis_walkthrough/common.R")
compare_figure("Fig1_study_design.png")
