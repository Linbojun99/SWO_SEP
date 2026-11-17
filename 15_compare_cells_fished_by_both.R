# =====================================================================================================
# Step 15. CPUE in cells fished by both fleets (Section 2.2, Section S9, Figs. 4, S11 and S12, Table S11)
# 第 15 步  两支船队均有作业的网格内的 CPUE（第 2.2 节，补充材料 S9，图 4、S11、S12，表 S11）
#
# Takes the relative CPUE of step 07: in each cell of each area of overlap with Spain, each fleet's predicted
# CPUE divided by its own mean over the 15 years in that cell, with the median and pointwise 95% interval of
# 1,000 joint draws. For each year it counts the cells fished by both fleets (cells with eligible records from
# both), the cells in which the two medians lie on opposite sides of 1 and the cells in which the two
# intervals lie wholly on opposite sides of 1 (Table S11). Fig. 4a and Fig. S12 map the ratio of the named
# fleet's relative CPUE to Spain's in 2019-2020 and 2021-2023, Fig. 4b shows the quarterly China/Spain change
# ratios of step 14, Fig. S11 maps each fleet's relative CPUE in 2019 and 2020, and the supplementary maps
# show Fig. S11 for every year, two years per page.
# 使用第 07 步的相对 CPUE：在与西班牙每个重叠区域的每个网格内，各船队的预测 CPUE 除以其在该网格 15 年中的均值，并
# 给出 1,000 次联合抽样的中位数和逐点 95% 区间。逐年统计两支船队均有作业的网格（两者都有合格记录的网格）、两个中位数
# 位于 1 两侧的网格，以及两个区间完全位于 1 两侧的网格（表 S11）。图 4a 和图 S12 绘制 2019-2020 年和 2021-2023 年所列
# 船队相对 CPUE 与西班牙相对 CPUE 之比，图 4b 显示第 14 步各季度中国与西班牙的变化比，图 S11 绘制 2019 和 2020 年各
# 船队的相对 CPUE，补充地图按每页两年给出每一年的图 S11。
#
# Reads / 读取:  analysis_walkthrough/output/fits/<flag>_cpue_map.csv (step 07), support.csv (step 04),
#                pair_overlap.csv (step 10), quarterly_shared_change.csv (step 14), config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/shared_cell_cpue_by_year.csv, shared_cell_direction_by_year.csv
#                (Table S11 / 表 S11), shared_cell_cpc_spain_ratios.csv, Fig4_local_contrasts.pdf/.png,
#                FigS11_shared_cell_cpue.pdf/.png, FigS12_local_contrasts_2021_2023.pdf/.png,
#                shared_cell_cpue_atlas.pdf (Supplementary maps / 补充地图)
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/15_compare_cells_fished_by_both.R
# =====================================================================================================

source("analysis_walkthrough/plot_functions.R")
cat("Step 15: CPUE in cells fished by both fleets (Figs. 4, S11, S12; Table S11)\n")
cat("第 15 步：两支船队均有作业的网格内的 CPUE（图 4、S11、S12；表 S11）\n")
calendar <- seq.int(cfg$years[1],cfg$years[2])

# ---- Relative CPUE of both fleets in each area of overlap / 各重叠区域内两支船队的相对 CPUE ----
# For each area of overlap, keep Spain's and the named fleet's relative CPUE in its cells; the record counts
# must equal the annual overlap of step 10.
# 对每个重叠区域，保留其网格内西班牙与所列船队的相对 CPUE；记录数必须与第 10 步的逐年重叠一致。
mapped <- do.call(rbind,lapply(cfg$flags,function(f) read.csv(file.path(out,"fits",paste0(f,"_cpue_map.csv")))))
stopifnot(!anyDuplicated(mapped[c("flag","cell","year")]),all(mapped$lower>0 & mapped$lower<mapped$ratio & mapped$ratio<mapped$upper))
overlaps <- read.csv(file.path(out,"pair_overlap.csv"))
support <- subset(read.csv(file.path(out,"support.csv")),domain!="ESP" & shared_years>=3)
paired_maps <- do.call(rbind,lapply(sort(unique(support$domain)),function(flag) {
  cells <- support$cell[support$domain==flag]
  pair <- do.call(rbind,lapply(c("ESP",flag),function(side) {
    sub <- mapped[mapped$flag==side & mapped$cell %in% cells,]
    stopifnot(all(table(sub$year)==length(cells)))
    sub$domain <- flag
    sub
  }))
  counts <- reshape(pair[c("year","cell","flag","records")],idvar=c("year","cell"),timevar="flag",direction="wide")
  counts[is.na(counts)] <- 0
  check <- overlaps[overlaps$fleet==flag,]
  both <- tapply(counts[[paste0("records.",flag)]]>0 & counts$records.ESP>0,counts$year,sum)
  stopifnot(all(both==check$both_cells[match(as.integer(names(both)),check$year)]))
  pair
}))
write.csv(paired_maps,file.path(out,"shared_cell_cpue_by_year.csv"),row.names=FALSE)

# ---- Directions by year (Table S11) / 逐年方向（表 S11） ----
# Count only where both fleets reported effort that year; years without records from one fleet remain missing.
# 只统计两支船队当年都报告努力量的网格；某船队没有记录的年份记为缺失。
skeleton <- merge(support[c("domain","cell","LonC5","LatC5")],expand.grid(year=calendar,side=c("Spain","Paired fleet")))
skeleton$flag <- ifelse(skeleton$side=="Spain","ESP",skeleton$domain)
pairs <- merge(skeleton,paired_maps,by=c("domain","cell","LonC5","LatC5","flag","year"),all.x=TRUE)
pair_cell <- interaction(pairs$domain,pairs$cell,pairs$year,drop=TRUE)
pairs$estimable <- ave(is.finite(pairs$ratio),pair_cell,FUN=all)
pairs$both <- ave(!is.na(pairs$records) & pairs$records>0,pair_cell,FUN=all)
pairs$fleet <- factor(pairs$domain,levels=names(fleets)[-1],labels=unname(fleets[-1]))
pairs$clear <- pairs$upper<1 | pairs$lower>1
stopifnot(all(pairs$estimable[pairs$both]))
by_pair <- do.call(rbind,lapply(split(pairs,interaction(pairs$domain,pairs$year,drop=TRUE)),function(d) {
  a <- subset(d,side=="Spain" & both)
  b <- subset(d,side=="Paired fleet" & both)
  b <- b[match(a$cell,b$cell),]
  stopifnot(identical(a$cell,b$cell))
  measured <- nrow(a)>0
  data.frame(domain=d$domain[1],year=d$year[1],paired_estimates=all(d$estimable),both_cells=nrow(a),
             opposite_medians=if (measured) sum((a$ratio-1)*(b$ratio-1)<0) else NA_integer_,
             spain_below_fleet_above=if (measured) sum(a$ratio<1 & b$ratio>1) else NA_integer_,
             opposite_pointwise_intervals=if (measured) sum((a$upper<1 & b$lower>1)|(a$lower>1 & b$upper<1)) else NA_integer_)
}))
write.csv(by_pair,file.path(out,"shared_cell_direction_by_year.csv"),row.names=FALSE)
print(subset(by_pair,year %in% 2019:2020)[c("domain","year","both_cells","opposite_medians","opposite_pointwise_intervals")],row.names=FALSE)

# ---- Fleet/Spain ratio of relative CPUE / 所列船队与西班牙相对 CPUE 之比 ----
# The ratio is calculated in cells fished by both fleets; a dot marks cells where one fleet's interval lies
# wholly above 1 and the other's wholly below 1.
# 比值在两支船队均有作业的网格内计算；一支船队的区间完全高于 1、另一支完全低于 1 的网格以点标出。
main_years <- 2019:2020
followup_years <- 2021:2023
map_years <- subset(pairs,year %in% c(main_years,followup_years))
cell_keys <- c("domain","cell","year","LonC5","LatC5")
spain <- subset(map_years,side=="Spain",c(cell_keys,"ratio","lower","upper"))
fleet_rows <- subset(map_years,side=="Paired fleet")
contrast <- merge(fleet_rows,spain,by=cell_keys,suffixes=c("","_spain"))
stopifnot(nrow(contrast)==nrow(fleet_rows), !anyDuplicated(contrast[cell_keys]))
contrast$supported <- with(contrast,both & ((lower>1 & upper_spain<1) | (upper<1 & lower_spain>1)))
observed <- subset(contrast,both)
observed$cpc_over_spain <- observed$ratio / observed$ratio_spain
stopifnot(all(is.finite(observed$cpc_over_spain)))
map_values <- with(observed,data.frame(domain,year,cell,LonC5,LatC5,cpc_relative_cpue=ratio,
                                       spain_relative_cpue=ratio_spain,cpc_over_spain,
                                       opposite_within_fleet_intervals=supported))
write.csv(map_values,file.path(out,"shared_cell_cpc_spain_ratios.csv"),row.names=FALSE)
ratio_map <- function(years) {
  d <- subset(contrast,year %in% years)
  local_rates <- subset(observed,year %in% years)
  missing_pairs <- unique(d[!d$estimable,c("domain","fleet","year")])
  missing_pairs$label <- with(missing_pairs,ifelse(domain=="VUT" & year==2020,
                              "Catch reported\nby weight only","No paired estimates"))
  ggplot(d,aes(LonC5,LatC5)) +
    geom_tile(data=subset(d,estimable & !both),width=5,height=5,fill="white",colour="grey60",linewidth=.16) +
    geom_tile(data=local_rates,aes(fill=log2(cpc_over_spain)),width=5,height=5,colour="grey60",linewidth=.16) +
    geom_vline(xintercept=-100,colour="grey52",linewidth=.32) +
    geom_hline(yintercept=-20,colour="grey52",linewidth=.32) +
    geom_point(data=subset(local_rates,supported),aes(shape="95% intervals on opposite sides of 1"),
               size=.8,colour="black",show.legend=TRUE) +
    map_land_layer +
    geom_label(data=missing_pairs,aes(x=-112,y=-16,label=label),inherit.aes=FALSE,
               linewidth=0,fill="white",size=2.46,lineheight=.95) +
    facet_grid(year~fleet,switch="y") +
    scale_fill_gradient2(name="Fleet/Spain",low="#326998",mid="#F4E9C9",high="#B74A28",
                         midpoint=0,limits=c(-3,3),oob=scales::squish,breaks=-3:3,
                         labels=c("≤0.125","0.25","0.5","1","2","4","≥8"),
                         guide=guide_colourbar(order=1,barwidth=grid::unit(78,"mm"),barheight=grid::unit(3.5,"mm"))) +
    scale_shape_manual(name=NULL,limits=c("95% intervals on opposite sides of 1","Missing records for one or both fleets"),
                       values=c("95% intervals on opposite sides of 1"=16,"Missing records for one or both fleets"=22),
                       guide=guide_legend(order=2,ncol=1,override.aes=list(size=c(1.7,2.8),
                                         colour=c("black","grey60"),fill=c("black","white")))) +
    coord_sf(default_crs=st_crs(4326),xlim=c(-151,-68),ylim=c(-42,11),expand=FALSE) +
    scale_x_continuous(breaks=c(-140,-100),labels=lon_labels) +
    scale_y_continuous(breaks=c(-30,0),labels=lat_labels) +
    labs(x=NULL,y=NULL) + map_theme +
    theme(axis.text=element_text(size=7),legend.box="horizontal",legend.spacing.y=grid::unit(0,"pt"),
          legend.margin=margin(0,0,0,0),legend.key.size=grid::unit(3.8,"mm"))
}

# ---- Fig. 4: (a) ratio maps for 2019 and 2020, (b) quarterly change ratios / 图 4：(a) 2019 和 2020 年比值地图，(b) 季度变化比 ----
quarter_results <- read.csv(file.path(out, "quarterly_shared_change.csv"))
quarter_results <- subset(quarter_results, model == "year_quarter" & series == "CHN_over_ESP")
quarter_results$period <- factor(quarter_results$period, levels = c(paste0("Q", 1:4), "equal_quarters"),
                                 labels = c("1st quarter", "2nd quarter", "3rd quarter", "4th quarter", "4-quarter mean"))
stopifnot(nrow(quarter_results) == 5L, !anyNA(quarter_results$period))
quarter_panel <- ggplot(quarter_results, aes(period, estimate)) +
  geom_hline(yintercept = 1, linewidth = .47, colour = "grey65") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .12, linewidth = .72) +
  geom_point(aes(shape = period == "4-quarter mean"), size = 2.1) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
  scale_x_discrete(labels = setNames(paste0(quarter_results$period, "\n(", quarter_results$n_cells, " cells)"),
                                      quarter_results$period)) +
  scale_y_log10(breaks = c(1, 2, 4, 8)) +
  labs(x = NULL, y = "Ratio of CPUE changes
(China/Spain)") + map_theme
panel_label <- theme(plot.title = element_text(face = "bold", size = 9, hjust = 0, margin = margin(b = 2)),
                     plot.title.position = "plot")
main_map <- ratio_map(main_years) + labs(title = "(a)") + panel_label
quarter_panel <- quarter_panel + labs(title = "(b)") + panel_label
figure4 <- patchwork::wrap_plots(main_map, quarter_panel, ncol = 1, heights = c(1.2, 1))
save_fig(figure4, "Fig4_local_contrasts", 120)

# ---- Fig. S12: ratio maps for 2021 to 2023 / 图 S12：2021 至 2023 年比值地图 ----
save_fig(ratio_map(followup_years),"FigS12_local_contrasts_2021_2023",98)

# ---- Fig. S11 and the supplementary maps: relative CPUE of each fleet / 图 S11 与补充地图：各船队的相对 CPUE ----
blue <- "#2166AC"; red <- "#B2182B"; middle <- "#F7F7F2"
cpue_map <- function(years) {
  d <- subset(pairs,year %in% years)
  columns <- as.vector(outer(c("Spain","Paired fleet"),years,paste,sep=", "))
  d$column <- factor(paste(d$side,d$year,sep=", "),levels=columns,
                     labels=sub("Paired fleet","Tuna-targeting fleet",columns,fixed=TRUE))
  missing <- unique(d[!d$estimable,c("domain","fleet","column","year","side")])
  observed <- subset(d,both)
  plot <- ggplot(d,aes(LonC5,LatC5)) +
    geom_tile(data=subset(d,estimable & !both),width=5,height=5,fill="grey82",colour="grey65",linewidth=.14) +
    geom_tile(data=observed,aes(fill=log2(ratio)),width=5,height=5,colour="grey65",linewidth=.15) +
    geom_vline(xintercept=-100,colour="grey52",linewidth=.32) +
    geom_hline(yintercept=-20,colour="grey52",linewidth=.32) +
    geom_point(data=subset(observed,clear),aes(shape="95% interval excludes 1"),size=.65,colour="black") +
    map_land_layer + facet_grid(fleet~column,switch="y") +
    scale_fill_gradient2(name="Standardized CPUE / 15-year mean in the cell",low=blue,mid=middle,high=red,
                         midpoint=0,limits=c(-2,2),oob=scales::squish,breaks=-2:2,
                         labels=c("≤0.25","0.5","1","2","≥4"),guide=guide_colourbar(barwidth=grid::unit(58,"mm"),order=1)) +
    scale_shape_manual(name=NULL,values=c("95% interval excludes 1"=16),
                       guide=guide_legend(order=2,override.aes=list(size=1.6))) +
    coord_sf(default_crs=st_crs(4326),xlim=c(-151,-68),ylim=c(-46,11),expand=FALSE) +
    scale_x_continuous(breaks=c(-140,-100,-80),labels=lon_labels) +
    scale_y_continuous(breaks=c(-40,0),labels=lat_labels) +
    labs(x=NULL,y=NULL) + map_theme + theme(axis.text=element_text(size=7))
  if (nrow(missing)) {
    missing$label <- "No paired estimates"
    weight_only <- missing$domain=="VUT" & missing$year==2020
    missing$label[weight_only & missing$side=="Spain"] <- "No Vanuatu\ncount index"
    missing$label[weight_only & missing$side=="Paired fleet"] <- "Catch reported\nby weight only"
    plot <- plot + geom_label(data=missing,aes(label=label),x=-115,y=-19,inherit.aes=FALSE,
                              linewidth=0,size=2.46,fill="white",label.padding=grid::unit(.6,"mm"))
  }
  plot
}
save_fig(cpue_map(main_years),"FigS11_shared_cell_cpue",175)
grDevices::cairo_pdf(file.path(out,"shared_cell_cpue_atlas.pdf"),width=190/25.4,height=175/25.4,onefile=TRUE)
for (first in seq(min(calendar),max(calendar),by=2)) {
  plot <- cpue_map(intersect(first:(first+1L),calendar)) +
    labs(title="Swordfish CPUE in the area of overlap with Spain",
         caption=paste("Dots mark 95% intervals wholly above (red) or below (blue) the fleet’s 15-year mean.",
                       "Grey cells lack records from one or both fleets that year.",sep="\n")) +
    theme(plot.title=element_text(size=10,face="bold"),plot.caption=element_text(size=7,hjust=0))
  print(plot)
}
invisible(grDevices::dev.off())

# ---- Compare with the archived results / 与存档结果比较 ----
source("analysis_walkthrough/common.R")
compare_table("shared_cell_cpue_by_year.csv", keys=c("domain","flag","cell","year"))
compare_table("shared_cell_direction_by_year.csv")
compare_table("shared_cell_cpc_spain_ratios.csv")
for (name in c("Fig4_local_contrasts.png","FigS11_shared_cell_cpue.png","FigS12_local_contrasts_2021_2023.png")) compare_figure(name)
cat("  shared_cell_cpue_atlas.pdf: not compared (a PDF records its creation time) / 未比较（PDF 记录了生成时间）\n")
