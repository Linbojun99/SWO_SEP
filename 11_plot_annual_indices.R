# =====================================================================================================
# Step 11. Annual indices in the core area and the areas of overlap (Fig. 3, Fig. S1)
# 第 11 步  核心作业区和各重叠区域的年度指数（图 3、图 S1）
#
# Fig. 3 maps the cells of each index (Spain's core area and each area of overlap) beside Spain's index and
# the named fleet's index, each scaled to its 15-year mean, with pointwise 95% confidence intervals. Open
# symbols mark years in which that fleet had no records in the cells shown. Fig. S1 shows Spain's index for the
# core-area cells in each assessment area beside the index of the whole core area.
# 图 3 左侧为各指数所用的网格（西班牙核心作业区和各重叠区域），右侧为西班牙指数和所列船队的指数，均除以各自的
# 15 年均值，并附逐点 95% 置信区间。空心符号标出该船队在所示网格内没有记录的年份。图 S1 显示各评估分区内核心作业区
# 网格的西班牙指数，并与整个核心作业区的指数并列。
# Same drawing as the Fig. 3 and Fig. S1 part of R/plot_directed_fishery.R.
# 与 R/plot_directed_fishery.R 中图 3 和图 S1 部分的绘图相同。
#
# Reads / 读取:  analysis_walkthrough/output/indices_baseline.csv (step 07), support.csv (step 04),
#                spain_area_indices.csv (step 08), config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/Fig3_annual.pdf/.png, FigS1_spain_areas.pdf/.png
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/11_plot_annual_indices.R
# =====================================================================================================

source("analysis_walkthrough/plot_functions.R")
cat("Step 11: drawing Fig. 3 and Fig. S1 / 第 11 步：绘制图 3 和图 S1\n")

# Baseline indices at the three-year threshold over Spain's core area and each area of overlap.
# 三年阈值下西班牙核心作业区和各重叠区域的基线指数。
idx <- read.csv(file.path(out,"indices_baseline.csv"))
idx$annual_only <- as.logical(idx$annual_only)
compare <- subset(idx,threshold==3 & !annual_only & model=="all" & region=="all")
compare$series <- ifelse(compare$flag!="ESP","Tuna-targeting fleet","Spain")
compare$series <- factor(compare$series,levels=c("Tuna-targeting fleet","Spain"))
panelnames <- index_panels(compare)
compare$panel <- factor(panelnames[match(compare$domain,cfg$flags)],levels=panelnames)
p2 <- plot_lines(compare,line_style(c("Spain","Tuna-targeting fleet"),
                                    c("black",npg["navy"]),c(16,17),c(21,24)))
# Name the years without an estimate for the fleet.
# 标出该船队没有估计的年份。
gaps <- data.frame(domain=cfg$flags,stringsAsFactors=FALSE)
gaps$label <- vapply(gaps$domain,function(f) {
  missing <- setdiff(seq.int(cfg$years[1],cfg$years[2]),compare$year[compare$flag==f & compare$domain==f])
  if(length(missing)) paste("Missing:",paste(missing,collapse=", ")) else ""
},character(1))
gaps$panel <- factor(panelnames,levels=panelnames)
p2 <- p2 + geom_text(data=subset(gaps,label!=""),aes(x=2005,y=Inf,label=label),inherit.aes=FALSE,
                      size=2.46,hjust=0,vjust=1.4,colour="grey35")

# Maps of the cells indexed: the core area and the areas of overlap, which follow effort, not positive catch.
# 所计算网格的地图：核心作业区和各重叠区域，按努力量而非正渔获确定。
support <- subset(read.csv(file.path(out,"support.csv")),shared_years>=3)
expected <- unique(compare[c("domain","domain_cells")])
counts <- table(factor(support$domain,levels=expected$domain))
stopifnot(!anyDuplicated(support[c("domain","cell")]),all(counts==expected$domain_cells),
          all(support$cell %in% support$cell[support$domain=="ESP"]))
support$panel <- factor(panelnames[match(support$domain,cfg$flags)],levels=panelnames)
map_labels <- sprintf("%s\n%d %s", fleets[cfg$flags],
  expected$domain_cells[match(cfg$flags,expected$domain)],
  ifelse(cfg$flags=="ESP","core-area cells","cells in the\narea of overlap"))
map3 <- ggplot(support,aes(LonC5,LatC5)) +
  geom_tile(width=5,height=5,fill=unname(npg["navy"]),colour="white",linewidth=.11) +
  geom_sf(data=st_sf(geometry=land_mask),inherit.aes=FALSE,fill="grey83",colour=NA) +
  facet_grid(panel~.,switch="y",labeller=labeller(panel=setNames(map_labels,panelnames))) +
  coord_sf(default_crs=st_crs(4326),xlim=c(-151,-68),ylim=c(-42,12),expand=FALSE) +
  scale_x_continuous(breaks=c(-140,-100,-80),labels=lon_labels) +
  scale_y_continuous(breaks=c(-30,0),labels=lat_labels) +
  labs(x=NULL,y=NULL,title="Fishing grounds") + fr +
  theme(strip.placement="outside",strip.text.y.left=element_text(angle=0,hjust=0),
        axis.text=element_text(size=7),axis.ticks.length=grid::unit(1.1,"pt"),
        panel.spacing.y=grid::unit(1.7,"mm"),plot.title=element_text(size=9,face="bold"))
annual3 <- p2 + facet_grid(panel~.) + labs(title="Relative swordfish index") +
  theme(strip.text.y=element_blank(),axis.title.y=element_blank(),
        panel.spacing.y=grid::unit(1.7,"mm"),plot.title=element_text(size=9,face="bold"))
figure3 <- (map3 | annual3) + plot_layout(widths=c(NA,2.6),guides="collect") &
  theme(legend.position="bottom") & legend_tight
save_fig(figure3,"Fig3_annual",181)

# ---- Fig. S1: Spain's index in each assessment area / 图 S1：各评估分区的西班牙指数 ----
area_idx <- subset(read.csv(file.path(out,"spain_area_indices.csv")),threshold==3)
core <- subset(compare,flag=="ESP" & domain=="ESP")
stopifnot(nrow(area_idx)==4L*nrow(core),all(area_idx$id=="ESP_all"))
area_names <- c(area2="oceanic north",area3="oceanic south",area4="coastal north",area5="coastal south")
area_panels <- sapply(names(area_names),function(a) sprintf("Area %s, %s (%d cells)",sub("area","",a),area_names[[a]],
  unique(area_idx$domain_cells[area_idx$region==a])))
stopifnot(sum(unique(area_idx[c("region","domain_cells")])$domain_cells)==unique(core$domain_cells))
area_idx$series <- "Spain, assessment area"
area_idx$panel <- area_panels[area_idx$region]
core_ref <- do.call(rbind,lapply(area_panels,function(pn) transform(core,panel=pn,series="Spain, core area")))
s5 <- rbind(area_idx[c("year","index","lower","upper","records","series","panel")],
            core_ref[c("year","index","lower","upper","records","series","panel")])
s5$panel <- factor(s5$panel,levels=area_panels)
s5$series <- factor(s5$series,levels=c("Spain, core area","Spain, assessment area"))
pS5 <- plot_lines(s5,line_style(c("Spain, core area","Spain, assessment area"),c("black",npg["green"]),
                                c(16,15),c(21,22)))
save_fig(pS5 + legend_tight,"FigS1_spain_areas",120)

# ---- Compare with the archived figures / 与存档图比较 ----
source("analysis_walkthrough/common.R")
for (name in c("Fig3_annual.png","FigS1_spain_areas.png")) compare_figure(name)
