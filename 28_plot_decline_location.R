# =====================================================================================================
# Step 28. Location of the decline (Fig. S14)
# 第 28 步  下降位置（图 S14）
#
# Panel (a) shows the median true and estimated percentage change from the early mean to the mean of 2021
# and 2023 for the three locations of the decline (A2, L1, L2) and the three designs, over the replicates in
# which all nine fits were accepted. Panel (b) shows the median reduction in the MAE of the log index when
# part of the added records came from outside the core area, with the same number of records. Bars are 95%
# bootstrap intervals. These are simulation results, not estimates from the fishery data.
# (a) 显示三种下降位置（A2、L1、L2）和三种设计下，由早期均值到 2021 和 2023 年均值的真实与估计百分比变化的中位数，
# 使用九个拟合均被接受的重复。(b) 显示记录数相同时，部分新增记录来自核心作业区以外使对数指数 MAE 减少的中位数。
# 误差线为 95% 自助法区间。这些是模拟结果，不是渔业数据的估计。
# Same drawing as R/plot_decline_location.R. / 与 R/plot_decline_location.R 的绘图相同。
#
# Reads / 读取:  analysis_walkthrough/output/decline_location_summary.csv (step 27)
# Writes / 写出: analysis_walkthrough/output/FigS14_decline_location.pdf, FigS14_decline_location.png
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/28_plot_decline_location.R
# =====================================================================================================

suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
out <- "analysis_walkthrough/output"
cat("Step 28: drawing Fig. S14 / 第 28 步：绘制图 S14\n")
tab <- read.csv(file.path(out,"decline_location_summary.csv"))

# Scenario codes of Table 1: A2 outside the core area, L1 core area, L2 both areas.
# 表 1 的情景代码：A2 核心作业区以外，L1 核心作业区，L2 两个区域。
placement <- c(periphery_decline="A2 Outside the core area only",core_decline="L1 Core area only",uniform_decline="L2 Both areas")
tab$location <- factor(tab$scenario,levels=rev(names(placement)),labels=rev(placement))
trend <- subset(tab,cohort=="three_designs" & measure=="change_percent")
gain <- subset(tab,cohort=="equal_records" & measure=="year_error" & series=="Core_minus_Mixed")
stopifnot(nrow(trend)==12L,nrow(gain)==3L,!anyNA(tab$location))
styles <- c(Truth="grey42",Spain="black",Core="#3C5488",Mixed="#E64B35")
symbols <- c(Truth=18,Spain=16,Core=17,Mixed=15)
series_labels <- c(Truth="True change",Spain="Spain only",Core="Spain + core-area records",
                   Mixed="Spain + records from both areas")
offset <- c(Truth=.26,Spain=.085,Core=-.085,Mixed=-.26)
trend$y <- as.integer(trend$location) + offset[trend$series]
trend$series <- factor(trend$series,levels=names(styles))
fr <- theme_bw(base_size=9,base_family="Arial") +
  theme(panel.grid=element_blank(),panel.border=element_rect(colour="grey25",linewidth=.5),
        legend.position="bottom", legend.title=element_blank(), plot.title=element_text(face="bold",size=9),
        plot.margin=margin(4,6,3,4))
y_scale <- scale_y_continuous(breaks=1:3,labels=rev(unname(placement)),limits=c(.49,3.51),
                              expand=expansion(mult=0))

# (a) Median true and estimated change with bootstrap intervals.
# (a) 真实与估计变化的中位数及自助法区间。
a <- ggplot(trend,aes(median,y,colour=series)) +
  geom_vline(xintercept=0,colour="grey75",linewidth=.45) +
  geom_errorbar(data=subset(trend,series!="Truth"),aes(xmin=lower,xmax=upper),
                orientation="y",width=.07,linewidth=.7,show.legend=FALSE) +
  geom_point(aes(shape=series),size=2.25) + y_scale + scale_x_continuous(labels=scales::label_number(style_negative="minus")) +
  scale_colour_manual(name="Series",values=styles,labels=series_labels) +
  scale_shape_manual(name="Series",values=symbols,breaks=names(styles),labels=unname(series_labels)) +
  guides(colour="none",shape=guide_legend(nrow=1,override.aes=list(colour=unname(styles)))) +
  labs(title="(a) Regional abundance change",x="Change from early to late years (%)",y="Location of decline") + fr

# (b) Median reduction in the MAE of the log index, core area only minus both areas.
# (b) 对数指数 MAE 减少量的中位数（仅核心作业区减去两个区域）。
b <- ggplot(gain,aes(median,as.integer(location))) +
  geom_vline(xintercept=0,colour="grey60",linewidth=.5) +
  geom_errorbar(aes(xmin=lower,xmax=upper),width=.10,orientation="y",colour="#3C5488",linewidth=.75) +
  geom_point(colour="#3C5488",size=2.35) + y_scale +
  labs(title="(b) Both areas instead of core area only",x="Reduction in MAE of the log index",y=NULL) + fr +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
plot <- (a | b) + plot_layout(widths=c(1.13,1),guides="collect") &
  theme(legend.position="bottom",legend.box.spacing=grid::unit(1.5,"mm"),legend.margin=margin(0,0,0,0))
for (ext in c("pdf","png")) {
  filename <- file.path(out,paste0("FigS14_decline_location.",ext))
  ggsave(filename,plot,width=190,height=86,units="mm",dpi=600,bg="white",
         device=if (ext=="pdf") grDevices::cairo_pdf else "png")
}

# ---- Compare with the archived figure / 与存档图比较 ----
source("analysis_walkthrough/common.R")
compare_figure("FigS14_decline_location.png")
