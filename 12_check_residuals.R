# =====================================================================================================
# Step 12. Fitted count distributions and randomized quantile residuals (Section 2.2, Sections S3 and S4)
# 第 12 步  拟合的尾数分布与随机分位数残差（第 2.2 节，补充材料 S3 和 S4）
#
# Prints the NB2 dispersion, the shared Matern range, the standard deviations of the two fields and the
# reported and fitted proportions of zeros of the six baseline models fitted to all eligible records (Table
# S2). Draws the quantiles of their randomized quantile residuals (Fig. S2) and maps the mean residual by
# fleet, 5-degree cell and period (Fig. S3). Calculates the Pearson correlations of the mean residuals of
# cells sharing an edge in the same year and of the same cell in consecutive years (Section S4).
# 输出用全部合格记录拟合的六个基线模型的 NB2 离散参数、共用 Matern 范围参数、两个随机场的标准差，以及报告与拟合的
# 零值比例（表 S2）。绘制其随机分位数残差的分位数（图 S2），并按船队、5° 网格和时段绘制平均残差地图（图 S3）。计算
# 同一年共边网格之间以及同一网格相邻年份之间平均残差的 Pearson 相关（补充材料 S4）。
#
# Reads / 读取:  analysis_walkthrough/output/fits/<id>_parameters.csv, <id>_diagnostics.csv, <id>_qq.csv and
#                <flag>_residual_cells.csv of the six baseline models (step 07), config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/spatial_residual_cell_years.csv, residual_pair_correlations.csv,
#                FigS2_fleet_residuals.pdf/.png, FigS3_spatial_residuals.pdf/.png
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/12_check_residuals.R
# =====================================================================================================

source("analysis_walkthrough/plot_functions.R")
cat("Step 12: fitted count distributions and residuals (Table S2, Figs. S2 and S3)\n")
cat("第 12 步：拟合的尾数分布与残差（表 S2、图 S2 和 S3）\n")
dir <- file.path(out, "fits")
ids <- paste0(cfg$flags, "_all")
read_parts <- function(names) do.call(rbind, lapply(names, function(x) read.csv(file.path(dir, x))))

# ---- Table S2 / 表 S2 ----
# Dispersion (phi), shared range (km), SDs of the spatial and spatiotemporal fields, and proportions of zeros (%).
# 离散参数（phi）、共用范围参数（千米）、空间与时空随机场标准差，以及零值比例（%）。
parameters <- read_parts(paste0(ids, "_parameters.csv"))
diagnostics <- read_parts(paste0(ids, "_diagnostics.csv"))
table_s2 <- reshape(subset(parameters, group=="ran_pars", c(id, term, estimate)), idvar="id", timevar="term", direction="wide")
table_s2 <- merge(table_s2, diagnostics[c("id","observed_zero","expected_zero","residual_sd")], by="id")
table_s2[c("observed_zero","expected_zero")] <- 100*table_s2[c("observed_zero","expected_zero")]
print(table_s2[match(ids, table_s2$id),], digits=3, row.names=FALSE)

# ---- Residuals by cell and year (Section S4) / 按网格和年份的残差（补充材料 S4） ----
cell_year <- read_parts(paste0(cfg$flags, "_residual_cells.csv"))
stopifnot(!anyDuplicated(cell_year[c("flag","cell","Year")]), all(cell_year$records>0))
write.csv(cell_year, file.path(out, "spatial_residual_cell_years.csv"), row.names=FALSE)
# Each pair links the mean residuals of two observed cell-years; an edge between cells enters once in each orientation.
# 每一对连接两个有观测的网格年份的平均残差；网格之间的每条共边按两个方向各计入一次。
cell_year$mean_residual <- cell_year$sum_residual/cell_year$records
pair_rows <- function(d,dlon=0,dlat=0,dy=0) {
  keys <- c("LonC5","LatC5","Year")
  shifted <- d[c(keys,"mean_residual")]
  shifted$LonC5 <- shifted$LonC5-dlon
  shifted$LatC5 <- shifted$LatC5-dlat
  shifted$Year <- shifted$Year-dy
  merge(d[c(keys,"mean_residual")],shifted,by=keys,suffixes=c("_from","_to"))
}
correlations <- do.call(rbind,lapply(split(cell_year,cell_year$flag),function(d) {
  space <- rbind(pair_rows(d,dlon=5),pair_rows(d,dlat=5))
  time <- pair_rows(d,dy=1)
  a <- space$mean_residual_from; b <- space$mean_residual_to
  data.frame(flag=d$flag[1],adjacent_cell_pairs=nrow(space),adjacent_cell_r=cor(c(a,b),c(b,a)),
             consecutive_year_pairs=nrow(time),consecutive_year_r=cor(time$mean_residual_from,time$mean_residual_to))
}))
write.csv(correlations,file.path(out,"residual_pair_correlations.csv"),row.names=FALSE)
print(correlations, digits=2, row.names=FALSE)

# ---- Fig. S3: mean residual by fleet, cell and period / 图 S3：按船队、网格和时段的平均残差 ----
# Periods 2005-2009, 2010-2014, 2015-2019 and 2020-2023; each monthly record has equal weight, and cells with
# one to four records are grey.
# 时段为 2005-2009、2010-2014、2015-2019 和 2020-2023 年；每条月度记录权重相等，只有一至四条记录的网格以灰色表示。
blue <- "#2166AC"; red <- "#B2182B"; middle <- "#F7F7F2"
study <- cfg$years
starts <- seq(study[1],study[2],by=5L)
ends <- pmin(starts+4L,study[2])
cell_year$period <- cut(cell_year$Year,breaks=c(study[1]-1L,ends),labels=paste(starts,ends,sep="–"))
stopifnot(!anyNA(cell_year$period))
period <- aggregate(cbind(records,sum_residual)~flag+period+LonC5+LatC5,cell_year,sum)
period$fleet <- factor(period$flag,levels=names(fleets),labels=unname(fleets))
period$mean_residual <- period$sum_residual/period$records
cat("Fleet, cell and period combinations with at least five records / 至少五条记录的船队、网格与时段组合:",
    sum(period$records>=5), "| mean residual beyond -1 or 1 / 平均残差超出 -1 或 1:",
    sum(period$records>=5 & abs(period$mean_residual)>1), "\n")
residual_map <- ggplot(period,aes(LonC5,LatC5)) +
  geom_tile(data=subset(period,records<5),width=5,height=5,fill="grey82",colour="grey70",linewidth=.12) +
  geom_tile(data=subset(period,records>=5),aes(fill=mean_residual),width=5,height=5,colour="grey70",linewidth=.08) +
  map_land_layer + facet_grid(fleet~period,switch="y") +
  scale_fill_gradient2(name="Mean randomized quantile residual",low=blue,mid=middle,high=red,midpoint=0,
                       limits=c(-1,1),oob=scales::squish,breaks=c(-1,-.5,0,.5,1),
                       labels=c("≤−1","−0.5","0","0.5","≥1"),guide=guide_colourbar(barwidth=grid::unit(50,"mm"))) +
  coord_sf(default_crs=st_crs(4326),xlim=c(-151,-68),ylim=c(-46,11),expand=FALSE) +
  scale_x_continuous(breaks=c(-140,-100,-80),labels=lon_labels) +
  scale_y_continuous(breaks=c(-40,0),labels=lat_labels) +
  labs(x=NULL,y=NULL) + map_theme + theme(axis.text=element_text(size=7))
save_fig(residual_map,"FigS3_spatial_residuals",198)

# ---- Fig. S2: quantiles of the randomized quantile residuals / 图 S2：随机分位数残差的分位数 ----
theme_set(theme_bw(base_size=9,base_family="Arial") +
  theme(panel.grid=element_blank(),panel.border=element_rect(colour="grey25",linewidth=.50),
        strip.background=element_blank(),strip.text=element_text(hjust=0,face="bold"),
        legend.position="bottom",legend.title=element_blank(),legend.key.width=grid::unit(1.8,"lines")))
baselines <- read_parts(paste0(ids, "_qq.csv"))
baselines$panel <- factor(sub("_all$","",baselines$id),levels=names(fleets),labels=fleets)
axis_limits <- range(c(baselines$normal,baselines$residual))
p <- ggplot(baselines,aes(normal,residual)) +
  geom_abline(slope=1,intercept=0,colour="grey55",linewidth=.63) +
  geom_point(colour="black",size=.8) +
  facet_wrap(~panel,ncol=3) +
  coord_equal(xlim=axis_limits,ylim=axis_limits) +
  labs(x="Standard normal quantile",y="Randomized quantile residual")
save_fig(p,"FigS2_fleet_residuals",125)

# ---- Compare with the archived results / 与存档结果比较 ----
source("analysis_walkthrough/common.R")
compare_table("spatial_residual_cell_years.csv")
compare_table("residual_pair_correlations.csv")
for (name in c("FigS2_fleet_residuals.png","FigS3_spatial_residuals.png")) compare_figure(name)
