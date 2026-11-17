# =====================================================================================================
# Step 32. Estimated and true regional abundance in the simulation (Fig. 6, Fig. S13)
# 第 32 步  模拟中估计与真实的区域丰度（图 6、图 S13）
#
# Fig. 6a shows, for the replicates in which all three fits were accepted, the median annual estimate, the
# interquartile range across replicates and the median true abundance, each relative to the early mean.
# Fig. 6b shows the median percentage change from the early mean to the mean of 2021 and 2023 with paired
# bootstrap 95% intervals at constant catchability (step 23) and with an annual increase of 1% or 2% in the
# catchability of the tuna-targeting longline fleets (step 30). Fig. S13 shows replicate 1 with pointwise
# 95% intervals. These are simulation results, not estimates from the fishery data.
# 图 6a 显示三个拟合均被接受的重复上，年度估计的中位数、各重复间的四分位距和真实丰度的中位数，均相对于早期均值。
# 图 6b 显示由早期均值到 2021 和 2023 年均值的百分比变化中位数及配对自助法 95% 区间，分别为可捕性恒定（第 23 步）
# 和金枪鱼延绳钓船队可捕性每年增加 1% 或 2%（第 30 步）。图 S13 显示第 1 次重复及其逐点 95% 区间。这些是模拟结果，
# 不是渔业数据的估计。
# Same drawing as the Fig. 6 and Fig. S13 part of R/plot_directed_fishery.R.
# 与 R/plot_directed_fishery.R 中图 6 和图 S13 部分的绘图相同。
#
# Reads / 读取:  analysis_walkthrough/output/coverage_simulation.csv, coverage_annual_intervals.csv (step 21),
#                coverage_annual_paired_summary.csv (step 22), coverage_change_intervals.csv (step 23),
#                catchability_trend_summary.csv (step 30), config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/Fig6_coverage_experiment.pdf/.png, FigS13_single_replicate.pdf/.png
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/32_plot_coverage_experiment.R
# =====================================================================================================

source("analysis_walkthrough/plot_functions.R")
cat("Step 32: drawing Fig. 6 and Fig. S13 / 第 32 步：绘制图 6 和图 S13\n")
designs <- c("Spain","Shared","Extended")
design_labels <- c(Spain="Spain only",Shared="+ core area",Extended="+ outside core area")
# Scenario codes A1 to A3 of Table 1.
# 表 1 的情景代码 A1 至 A3。
scenario_labels <- c(stable="A1 Stable",periphery_decline="A2 Decline outside core area",
                     periphery_fluctuation="A3 Decline and recovery\noutside core area")
scenario_factor <- function(x) factor(x,levels=names(scenario_labels),labels=scenario_labels)
trials <- read.csv(file.path(out,"coverage_simulation.csv"))
coverage_stats <- read.csv(file.path(out,"coverage_change_intervals.csv"))

# The summary of step 30 must hold, at constant catchability, the same medians, intervals and replicates as step 23.
# 第 30 步的汇总在可捕性恒定时，必须与第 23 步有相同的中位数、区间和重复。
trend <- read.csv(file.path(out,"catchability_trend_summary.csv"))
stopifnot(nrow(trend)==3L*3L*3L,!anyDuplicated(trend[c("scenario","annual_q_increase","design")]))
base <- subset(trend,annual_q_increase==0)
ref <- coverage_stats[match(paste(base$scenario,base$design),paste(coverage_stats$scenario,coverage_stats$design)),]
ref_truth <- subset(coverage_stats,design=="Truth")
stopifnot(!anyNA(ref$change), all(abs(base$change_pct-ref$change)<=1e-8), all(abs(base$change_lower-ref$lower)<=1e-8),
          all(abs(base$change_upper-ref$upper)<=1e-8), all(base$replicates==ref$pairs),
          all(abs(trend$true_change_pct-ref_truth$change[match(trend$scenario,ref_truth$scenario)])<=1e-8))

# ---- Annual panels / 年度面板 ----
series <- read.csv(file.path(out,"coverage_annual_intervals.csv"))
example <- subset(series,rep_id==1)
stopifnot(nrow(example)==8L*sum(trials$rep_id==1 & trials$status=="accepted"),
          all(table(example$scenario,example$design)==8L))
truth_colour <- "grey42"
annual_panel <- function(data,title,band_label,truth_linetype) {
  data$panel <- scenario_factor(data$scenario)
  data$design <- factor(data$design,levels=names(design_labels),labels=design_labels)
  ggplot(data,aes(year,estimate)) +
    facet_grid(design~panel,switch="y") +
    geom_ribbon(aes(ymin=lower,ymax=upper,fill=band_label),alpha=.20,colour=NA) +
    geom_line(aes(colour="Estimate",linetype="Estimate"),linewidth=.82) +
    geom_line(aes(y=truth,colour="Truth",linetype="Truth"),linewidth=.72) +
    geom_point(colour=npg["navy"],size=1.2) +
    scale_colour_manual(name="Series",labels=c(Estimate="Estimated index",Truth="True abundance"),
                        values=c(Estimate=unname(npg["navy"]),Truth=truth_colour)) +
    scale_linetype_manual(name="Series",labels=c(Estimate="Estimated index",Truth="True abundance"),
                          values=c(Estimate="solid",Truth=truth_linetype)) +
    scale_fill_manual(name=NULL,values=setNames(unname(npg["navy"]),band_label)) +
    guides(colour=guide_legend(order=1),linetype=guide_legend(order=1),fill=guide_legend(order=2)) +
    scale_x_continuous(breaks=c(2009,2013,2017,2023),expand=expansion(add=1.25)) +
    scale_y_continuous(breaks=scales::breaks_pretty(n=4)) +
    labs(title=title,x="Year",y="Relative abundance") + fr +
    theme(strip.text.y.left=element_text(angle=0,hjust=1,face="plain"),
          strip.placement="outside",legend.key.width=grid::unit(1.8,"lines"),
          legend.key.height=grid::unit(.7,"lines"))
}
# Fig. 6a: medians and interquartile ranges over the complete replicates. Fig. S13: replicate 1.
# 图 6a：完整重复上的中位数和四分位距。图 S13：第 1 次重复。
annual_stats <- read.csv(file.path(out,"coverage_annual_paired_summary.csv"))
expected_pairs <- unique(coverage_stats[c("scenario","pairs")])
stopifnot(nrow(annual_stats)==8L*length(designs)*nrow(expected_pairs),
          all(annual_stats$pairs==expected_pairs$pairs[match(annual_stats$scenario,expected_pairs$scenario)]))
annual <- annual_panel(transform(annual_stats,estimate=median_estimate,truth=median_truth,lower=q25,upper=q75),
                       "(a)  Median annual indices","Interquartile range",truth_linetype="dashed")
single <- annual_panel(transform(example,estimate=estimated_index,truth=true_index),
                       NULL,"Pointwise 95% CI",truth_linetype="dashed")
save_fig(single + legend_tight,"FigS13_single_replicate",129)

# ---- Fig. 6b: change at annual catchability increases of 0%, 1% and 2% / 图 6b：可捕性年增幅 0%、1% 和 2% 时的变化 ----
# The Spain-only estimate does not depend on the rate, so that row has its 0% point only.
# 仅西班牙的估计与增幅无关，因此该行只有 0% 的点。
rates <- c("0"="0%","0.01"="1%","0.02"="2%")
trend$panel <- scenario_factor(trend$scenario)
trend$rate <- factor(rates[as.character(trend$annual_q_increase)],levels=rates)
trend_rows <- rev(c("Truth",names(design_labels)))
est <- subset(trend,design!="Spain" | annual_q_increase==0)
stopifnot(nrow(est)==3L*(1L+2L*3L))
# Rates share a row unless two intervals, widened by 1.2 percentage points for the symbols, overlap; such a row is dodged vertically.
# 各增幅共用一行，除非两个区间（为符号各加宽 1.2 个百分点）重叠；此时该行在垂直方向错开。
overlap <- function(lo,hi) { o <- order(lo); any(head(hi[o],-1)+1.2 > tail(lo[o],-1)-1.2) }
crowded <- vapply(split(est,list(est$scenario,est$design),drop=TRUE),
                  function(x) nrow(x)>1L && overlap(x$change_lower,x$change_upper),logical(1))
est$dodge <- ifelse(crowded[paste(est$scenario,est$design,sep=".")],c(.2,0,-.2)[as.integer(est$rate)],0)
est$y <- match(est$design,trend_rows) + est$dodge
trend_truth <- unique(trend[c("scenario","panel","true_change_pct")])
trend_truth$y <- match("Truth",trend_rows)
# Value labels for the true change and the 0% estimates, at least five percentage points from the zero line.
# 真实变化和 0% 估计的数值标签，距零线至少五个百分点。
fmt <- scales::label_number(accuracy=.1,suffix="%",style_positive="plus",style_negative="minus")
labels <- rbind(data.frame(panel=trend_truth$panel,change=trend_truth$true_change_pct,y=trend_truth$y),
                with(subset(est,annual_q_increase==0),data.frame(panel,change=change_pct,y)))
labels$label <- fmt(labels$change)
labels$label_x <- with(labels,ifelse(change < 0,pmin(change,-5),pmax(change,5)))
rate_title <- "Annual catchability increase, tuna-targeting fleets"
p4 <- ggplot(est,aes(change_pct,y)) + facet_wrap(~panel,nrow=1) +
  geom_vline(xintercept=0,colour="grey75",linewidth=.42) +
  geom_errorbar(aes(xmin=change_lower,xmax=change_upper,colour=rate,linewidth=rate),orientation="y",width=.14) +
  geom_point(aes(colour=rate,shape=rate,size=rate)) +
  geom_point(data=trend_truth,aes(true_change_pct,y),colour=truth_colour,shape=18,size=3.1) +
  geom_text(data=labels,aes(x=label_x,y=y+.40,label=label),size=2.8,colour="grey20") +
  scale_colour_manual(values=setNames(unname(npg[c("navy","green","red")]),rates),name=rate_title) +
  scale_shape_manual(values=setNames(c(16,17,15),rates),name=rate_title) +
  scale_size_manual(values=setNames(c(2.3,1.9,1.8),rates),name=rate_title) +
  scale_linewidth_manual(values=setNames(c(.78,.6,.6),rates),name=rate_title) +
  scale_x_continuous(breaks=scales::breaks_width(10),labels=scales::label_number(style_negative="minus"),
                     expand=expansion(mult=.12)) +
  scale_y_continuous(breaks=seq_along(trend_rows),labels=rev(c("True change",unname(design_labels))),
                     limits=c(.35,length(trend_rows)+.65),expand=expansion(mult=0)) +
  labs(title="(b)  Change in regional abundance with constant and increasing catchability",
       x="Regional abundance change (%)",y=NULL) +
  fr + theme(legend.position="bottom",legend.direction="horizontal",legend.title=element_text(size=9),
             legend.title.position="left",legend.key.width=grid::unit(1.8,"lines")) + legend_tight
# Panel (b) does not reserve aligned space on its left, so its axis stays next to the panel.
# (b) 左侧不保留对齐空间，使其坐标轴紧靠面板。
figure6 <- (annual + legend_tight) / free(p4,type="space",side="l") + plot_layout(heights=c(2.02,1.15)) &
  theme(plot.title=element_text(size=9,face="bold",margin=margin(b=3)))
save_fig(figure6,"Fig6_coverage_experiment",190)

# ---- Compare with the archived figures / 与存档图比较 ----
source("analysis_walkthrough/common.R")
for (name in c("Fig6_coverage_experiment.png","FigS13_single_replicate.png")) compare_figure(name)
