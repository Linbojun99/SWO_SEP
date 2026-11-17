# =====================================================================================================
# Step 05. Cells with positive retained swordfish counts (Section 2.1, Fig. 2)
# 第 05 步  有剑鱼留存正尾数的网格（第 2.1 节，图 2）
#
# Maps the 5-degree cells with a positive retained swordfish count in any month from 2005 to 2023: the cells of
# the Spanish fleet (a) and of all six fleets (b), with the cells where only the five tuna-targeting longline
# fleets reported swordfish in red. Only the marine portion of a cell is shaded; the cells are clipped at the
# Natural Earth 1:50m coastline and the marine fraction of each cell is written to coastal_grid_check.csv.
# 绘制 2005-2023 年任一月份剑鱼留存尾数为正的 5° 网格：西班牙船队的网格（a）和全部六支船队的网格（b），只有五支
# 金枪鱼延绳钓船队报告剑鱼的网格以红色表示。只为网格的海洋部分着色；网格按 Natural Earth 1:50m 海岸线裁剪，各网格的
# 海洋面积比例写入 coastal_grid_check.csv。
# Same drawing as the Fig. 2 part of R/plot_directed_fishery.R. / 与 R/plot_directed_fishery.R 中图 2 部分的绘图相同。
#
# Reads / 读取:  analysis_walkthrough/output/footprint.csv (step 04), config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/Fig2_coverage.pdf/.png, coastal_grid_check.csv
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/05_map_catch_cells.R
# =====================================================================================================

source("analysis_walkthrough/plot_functions.R")
cat("Step 05: drawing the cells with swordfish catches (Fig. 2) / 第 05 步：绘制有剑鱼渔获的网格（图 2）\n")

# Cells with Spanish catches, and cells where only the tuna-targeting longline fleets reported swordfish.
# 有西班牙渔获的网格，以及只有金枪鱼延绳钓船队报告剑鱼的网格。
foot <- read.csv(file.path(out,"footprint.csv"))
sp <- subset(foot,Flag=="ESP",select=c(cell,LonC5,LatC5))
others <- unique(subset(foot,Flag!="ESP",select=c(cell,LonC5,LatC5)))
add <- subset(others,!cell %in% sp$cell)
cat("Cells with Spanish catches / 有西班牙渔获的网格:", nrow(sp),
    "| cells added by the other fleets / 其他船队新增的网格:", nrow(add), "\n")
sp$class <- "Spain"; add$class <- "Tuna-targeting longline fleets only"
mapdata <- rbind(transform(sp,panel="(a)  Spanish swordfish fleet"),
                 transform(rbind(sp,add),panel="(b)  All six fleets"))

# Clip each 5-degree cell at the coastline; a cell centre is not a fishing location.
# 按海岸线裁剪每个 5° 网格；网格中心并非作业位置。
locations <- unique(mapdata[c("cell","LonC5","LatC5")])
boxes <- lapply(seq_len(nrow(locations)),function(i) {
  x <- locations$LonC5[i]; y <- locations$LatC5[i]
  st_polygon(list(cbind(x+c(-2.5,2.5,2.5,-2.5,-2.5),y+c(-2.5,-2.5,2.5,2.5,-2.5))))
})
rectangles <- st_sf(locations,geometry=st_sfc(boxes,crs=4326))
marine <- suppressWarnings(suppressMessages(st_difference(rectangles,land_mask)))
stopifnot(setequal(marine$cell,locations$cell),!any(st_is_empty(marine)))
centres <- st_as_sf(locations,coords=c("LonC5","LatC5"),crs=4326)
on_land <- lengths(suppressMessages(st_intersects(centres,land_mask)))>0
check <- locations
check$centre_on_land <- on_land
area_m2 <- function(x) as.numeric(st_area(st_transform(x,6933)))
check$marine_fraction <- area_m2(marine)[match(check$cell,marine$cell)]/area_m2(rectangles)
stopifnot(all(check$marine_fraction>0),max(check$marine_fraction)<1.00001)
write.csv(check,file.path(out,"coastal_grid_check.csv"),row.names=FALSE)
marine <- merge(marine,unique(mapdata[c("cell","class","panel")]),by="cell")

# Draw both panels with the northern limit of H2 at 10N and the assessment areas 2 to 5 split at 100W and 20S.
# 绘制两个面板，标出 H2 北界 10°N，以及以 100°W 和 20°S 划分的评估分区 2 至 5。
p1 <- ggplot() + geom_sf(data=land,fill="grey86",colour=NA) +
  geom_sf(data=marine,aes(fill=class),colour="white",linewidth=.13) +
  geom_hline(yintercept=10,linewidth=.55,colour="grey35") +
  geom_segment(data=data.frame(x=c(-100,-151),xend=c(-100,-70.2),y=c(-52,-20),yend=c(10,-20)),
               aes(x=x,xend=xend,y=y,yend=yend),linewidth=.38,colour="grey30") +
  geom_label(data=data.frame(x=c(-125,-125,-88,-88),y=c(-6,-33,-11,-33),label=c("2","3","4","5")),
             aes(x=x,y=y,label=label),size=2.6,fill="white",colour="black",linewidth=0,
             label.r=grid::unit(0,"pt"),label.padding=grid::unit(1.2,"pt")) +
  annotate("text",x=-144,y=14,label="Area 1",size=2.6) +
  annotate("text",x=-85,y=12.8,label="H2 limit",size=2.6) +
  scale_fill_manual(values=c("Spain"="grey60","Tuna-targeting longline fleets only"=unname(npg["red"])),
                    breaks=c("Spain","Tuna-targeting longline fleets only")) +
  facet_wrap(~panel,nrow=1) +
  coord_sf(xlim=c(-151,-69),ylim=c(min(mapdata$LatC5)-3.5,17),expand=FALSE) +
  scale_x_continuous(breaks=c(-140,-120,-100,-80),labels=lon_labels) +
  scale_y_continuous(breaks=c(-40,-20,0,10),labels=lat_labels) +
  labs(x=NULL,y=NULL) + fr + theme(strip.text=element_text(size=9))
save_fig(p1,"Fig2_coverage",95)

# ---- Compare with the archived results / 与存档结果比较 ----
source("analysis_walkthrough/common.R")
compare_table("coastal_grid_check.csv")
compare_figure("Fig2_coverage.png")
