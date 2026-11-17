# =====================================================================================================
# Shared drawing settings for the figures of the paper
# 论文各图共用的绘图设定
#
# Sourced by the plotting steps. Defines the fleet names, the colours (Nature Publishing Group palette of
# ggsci), the plot theme, the axis labels in degrees, the Natural Earth 1:50m land, and plot_lines(), which
# draws annual indices with their pointwise 95% intervals (Figs. 3, 5, S1 and S9).
# 由各绘图步骤读入。定义船队名称、颜色（ggsci 的 Nature Publishing Group 配色）、绘图主题、以度为单位的坐标轴
# 标注、Natural Earth 1:50m 陆地，以及 plot_lines()：绘制年度指数及其逐点 95% 区间（图 3、5、S1 和 S9）。
# Same drawing as the first part of R/plot_directed_fishery.R and of R/plot_spatial_fleet_checks.R.
# 与 R/plot_directed_fishery.R 和 R/plot_spatial_fleet_checks.R 开头部分的绘图相同。
#
# Reads / 读取:  config/directed_fishery.json
# Writes / 写出: nothing by itself / 本身不写出文件
# =====================================================================================================

suppressPackageStartupMessages({library(ggplot2); library(sf); library(patchwork)})
out <- "analysis_walkthrough/output"
cfg <- jsonlite::fromJSON("config/directed_fishery.json")
fleets <- c(ESP="Spain",JPN="Japan",KOR="Korea",CHN="China",TWN="Chinese Taipei",VUT="Vanuatu")
npg <- c(red="#E64B35",cyan="#4DBBD5",green="#00A087",navy="#3C5488",salmon="#F39B7F",slate="#8491B4")
stopifnot(identical(unname(npg),sub("FF$","",ggsci::pal_npg("nrc")(6))))

# Each series is a solid line; white-filled versions of its symbol mark years without records in the cells indexed.
# 每条序列为实线；同形状的白色空心符号标出所计算网格内没有记录的年份。
line_style <- function(series,colour,shape,open)
  data.frame(series,colour=unname(colour),shape,open,stringsAsFactors=FALSE)
composition_style <- line_style(c("Baseline","Cluster","PC"),c("black",npg["red"],npg["navy"]),
                                c(16,17,15),c(21,24,22))
fr <- theme_bw(base_size=9,base_family="Arial") +
  theme(panel.grid=element_blank(), panel.border=element_rect(colour="grey25",linewidth=.50),
        strip.background=element_blank(),strip.text=element_text(hjust=0,face="bold"),
        legend.position="bottom",legend.title=element_blank(),legend.key.width=grid::unit(1.1,"lines"),
        plot.margin=margin(4,5,3,5))
legend_tight <- theme(legend.box.spacing=grid::unit(1.5,"mm"),legend.margin=margin(0,0,0,0))
# Save each figure as a 190 mm wide PDF and a 600 dpi PNG.
# 每幅图保存为 190 毫米宽的 PDF 和 600 dpi 的 PNG。
save_fig <- function(p,name,h) {
  ggsave(file.path(out,paste0(name,".pdf")),p,width=190,height=h,units="mm",device=cairo_pdf)
  ggsave(file.path(out,paste0(name,".png")),p,width=190,height=h,units="mm",dpi=600,bg="white")
}
lon_labels <- function(x) ifelse(x<0,paste0(-x,"°W"),ifelse(x>0,paste0(x,"°E"),"0°"))
lat_labels <- function(y) ifelse(y<0,paste0(-y,"°S"),ifelse(y>0,paste0(y,"°N"),"0°"))

# Land from Natural Earth 1:50m country polygons, dissolved, for the maps of Figs. 2 and 3.
# Natural Earth 1:50m 国家多边形合并后的陆地，用于图 2 和图 3 的地图。
land <- rnaturalearth::ne_countries(scale=50,returnclass="sf")
sf_use_s2(FALSE)
coast <- suppressMessages(suppressWarnings(st_crop(st_make_valid(land),xmin=-155,xmax=-65,ymin=-59,ymax=22)))
land_mask <- suppressMessages(st_union(st_geometry(coast)))

# Annual indices with pointwise 95% intervals on a logarithmic scale. A series is broken where a year has no estimate.
# 对数尺度上的年度指数及其逐点 95% 区间。没有估计的年份处序列断开。
plot_lines <- function(data,style) {
  data <- data[order(data$panel,data$series,data$year),]
  data$segment <- ave(data$year,data$panel,data$series,FUN=function(y)cumsum(c(1,diff(y)!=1)))
  p <- ggplot(data,aes(year,index,colour=series,group=interaction(panel,series,segment))) +
    facet_wrap(~panel,ncol=2) + scale_colour_manual(values=setNames(style$colour,style$series)) +
    scale_shape_manual(values=setNames(style$shape,style$series)) +
    scale_y_log10() + scale_x_continuous(breaks=unique(c(seq(cfg$years[1],cfg$years[2],by=5),cfg$years[2])),limits=cfg$years+c(-.5,.6)) +
    labs(x="Year",y="Relative swordfish index") + fr
  stopifnot(all(is.finite(data$lower) & data$lower>0 & data$lower<=data$index & data$index<=data$upper),
            all(is.finite(data$upper)))
  p <- p + geom_ribbon(aes(ymin=lower,ymax=upper,fill=series),colour=NA,alpha=.17,show.legend=FALSE) +
    scale_fill_manual(values=setNames(style$colour,style$series),guide="none") +
    geom_line(linewidth=.6) +
    geom_point(data=subset(data,records>0),aes(shape=series),size=1.3)
  for (i in seq_len(nrow(style))) {
    zero <- subset(data,records==0 & series==style$series[i])
    if (nrow(zero)) p <- p + geom_point(data=zero,shape=style$open[i],fill="white",size=1.56,stroke=.55,show.legend=FALSE)
  }
  p
}

# Panel names "Fleet (n cells)" for the baseline indices of Fig. 3, in the order of the fleets.
# 图 3 基线指数的面板名称"船队（n 个网格）"，按船队顺序排列。
index_panels <- function(compare)
  sapply(cfg$flags,function(f) sprintf("%s (%d cells)",fleets[f],
    unique(compare$domain_cells[compare$domain==f & compare$region=="all"])))

# ---- Map style of Figs. 4, S3, S11 and S12 / 图 4、S3、S11 和 S12 的地图样式 ----
map_land_mask <- st_sf(geometry=suppressMessages(st_union(st_geometry(
  suppressWarnings(suppressMessages(st_crop(st_make_valid(land),xmin=-158,xmax=-65,ymin=-64,ymax=16)))))))
map_land_layer <- geom_sf(data=map_land_mask,inherit.aes=FALSE,colour=NA,fill="grey80")
map_theme <- theme_bw(base_size=9,base_family="Arial") +
  theme(panel.grid=element_blank(),panel.border=element_rect(colour="grey35",linewidth=.48),
        strip.background=element_blank(),strip.text=element_text(face="bold"),
        strip.text.y.left=element_text(angle=0,hjust=1),strip.placement="outside",
        legend.position="bottom",legend.title.position="top",
        panel.spacing=grid::unit(1.5,"mm"),plot.margin=margin(4,4,3,4))
