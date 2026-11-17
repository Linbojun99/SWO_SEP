# =====================================================================================================
# Step 08. Spain's index for the core-area cells in each assessment area (Section 2.1, Fig. S1)
# 第 08 步  各评估分区内核心作业区网格的西班牙指数（第 2.1 节，图 S1）
#
# Refits Spain's baseline model fitted to all eligible records (row 1 of tasks.csv) with the settings of step
# 07 and predicts it over the core-area cells of each assessment area of H2 (areas 2 to 5), with the cell
# weights, the 15-year mean and the interval method of the indices in Fig. 3. The refit must reproduce the
# indices of step 07 for the whole core area and east and west of 100W. A simulation check of the standard
# error from 1,000 joint draws is kept beside the delta-method standard error.
# 以第 07 步的设置重新拟合用全部合格记录拟合的西班牙基线模型（tasks.csv 第 1 行），并在 H2 各评估分区（分区 2 至 5）
# 的核心作业区网格上预测，网格权重、15 年均值和区间方法与图 3 的指数相同。重新拟合必须重现第 07 步在整个核心作业区
# 以及 100°W 以东和以西的指数。另保留由 1,000 次联合抽样得到的标准误，作为对 delta 法标准误的模拟核对。
#
# Reads / 读取:  analysis_walkthrough/output/input.rds (step 04), indices_baseline.csv (step 07),
#                config/directed_fishery.json, config/targeting_public_h2.json
# Writes / 写出: analysis_walkthrough/output/spain_area_indices.csv
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/08_spain_area_indices.R
# =====================================================================================================

source("analysis_walkthrough/model_functions.R")
out <- "analysis_walkthrough/output"
cat("Step 08: Spain's index by assessment area (Fig. S1) / 第 08 步：各评估分区的西班牙指数（图 S1）\n")
z <- readRDS(file.path(out, "input.rds"))
cfg <- jsonlite::fromJSON("config/directed_fishery.json")
h2 <- jsonlite::fromJSON("config/targeting_public_h2.json")$assessment_areas
assessment_area <- function(lon,lat) {
  stopifnot(!any(lon==h2$split_lon),!any(lat==h2$split_lat))
  ifelse(lat>=10,1L,ifelse(lon<h2$split_lon,ifelse(lat>h2$split_lat,2L,3L),ifelse(lat>h2$split_lat,4L,5L)))
}

# Refit Spain's baseline model with the data, mesh and settings of step 07.
# 以第 07 步的数据、三角网格和设置重新拟合西班牙基线模型。
i <- match("ESP_all",z$tasks$id); stopifnot(i==1L)
flag <- "ESP"; id <- "ESP_all"
d <- subset(z$d,Flag==flag)
years <- sort(unique(d$Year)); months <- sort(unique(d$Month))
d$year_f <- factor(d$Year,levels=years); d$month_f <- factor(d$Month,levels=months)
d$t <- match(d$Year,years)
g <- subset(z$grid,cell %in% d$cell)
geo <- fmesher::fm_mesh_2d(loc=as.matrix(g[c("x","y")]),cutoff=100,max.edge=c(850,1800),offset=c(650,1400))
mesh <- make_mesh(d,c("x","y"),mesh=geo)
fit <- sdmTMB(reformulate(c("year_f","month_f"),"SWOn"),data=d,mesh=mesh,time="t",
  family=nbinom2(),offset=log(d$Hooks/1000),spatial="on",spatiotemporal="iid",
  share_range=TRUE, control=sdmTMBcontrol(backend="tmb",newton_loops=2,get_joint_precision=TRUE))
gradient <- max(abs(fit$tmb_obj$gr(fit$model$par)))
stopifnot(fit$model$convergence==0,isTRUE(fit$sd_report$pdHess),gradient<cfg$gradient_limit)

# Predict every cell of the core area in every year for the earliest calendar month, and draw 1,000 joint
# realizations for the simulation check (seed cfg$seed + 1).
# 按最早的日历月份预测核心作业区每年的每个网格，并为模拟核对抽取 1,000 次联合实现（种子 cfg$seed + 1）。
masks <- subset(z$coverage,domain %in% cfg$flags & shared_years>=2)
pg <- subset(z$grid,cell %in% masks$cell)
p <- pg[rep(seq_len(nrow(pg)),times=length(years)),]
p$Year <- rep(years,each=nrow(pg)); p$t <- match(p$Year,years)
p$year_f <- factor(p$Year,levels=years); p$month_f <- factor(min(months),levels=months)
set.seed(cfg$seed+i)
sims <- exp(predict(fit,newdata=p,type="link",offset=rep(0,nrow(p)),nsim=cfg$draws))
stopifnot(all(dim(sims)==c(nrow(p),cfg$draws)),all(is.finite(sims)))
p$assessment_area <- assessment_area(p$LonC5,p$LatC5)

# Index of the whole core area, east and west of 100W, and of each assessment area, at each threshold.
# 在每个阈值下计算整个核心作业区、100°W 以东和以西以及各评估分区的指数。
series <- list()
for (threshold in 2:4) {
  target <- masks$cell[masks$domain=="ESP" & masks$shared_years>=threshold]
  for (region in c("all","east","west",paste0("area",2:5))) {
    keep <- p$cell %in% target
    if (region=="east") keep <- keep & p$LonC5>cfg$east_split
    if (region=="west") keep <- keep & p$LonC5<cfg$east_split
    if (startsWith(region,"area")) keep <- keep & p$assessment_area==as.integer(sub("area","",region))
    sub <- p[keep,]; target_cells <- unique(sub$cell)
    w <- sub$area/sum(pg$area[pg$cell %in% target_cells])
    bc <- bias_corrected_index(fit,sub,w)
    stopifnot(identical(bc$t,seq_along(years)))
    pt <- bc$est
    draws <- rowsum(sims[keep,,drop=FALSE]*w,sub$Year)
    ref <- years %in% z$refs
    idx <- pt/mean(pt[ref])
    ci <- cbind(bc$lower,bc$upper)/mean(pt[ref])
    rel <- sweep(draws,2,colMeans(draws[ref,,drop=FALSE]),"/")
    stopifnot(abs(mean(idx[ref])-1)<1e-10,max(abs(colMeans(rel[ref,,drop=FALSE])-1))<1e-10,
              max(abs((log(ci[,2])-log(ci[,1]))/(2*qnorm(.975))-bc$se))<1e-10)
    observed <- d[d$cell %in% target_cells,]
    count <- table(factor(observed$Year,levels=years))
    obs_cells <- unique(observed[c("cell","Year")])
    cell_count <- table(factor(obs_cells$Year,levels=years))
    series[[length(series)+1L]] <- data.frame(id,flag,model="all",annual_only=FALSE,domain="ESP",region,threshold,
      domain_cells=length(target_cells),year=years,records=as.integer(count),observed_cells=as.integer(cell_count),
      index=idx,lower=ci[,1],upper=ci[,2],cpue=pt,
      sd_log=(log(ci[,2])-log(ci[,1]))/(2*1.959964),sd_log_draws=apply(log(rel),1,sd),row.names=NULL)
  }
}
series <- do.call(rbind,series)

# The refit reproduces Spain's indices of step 07 for the whole core area and east and west of 100W.
# 重新拟合重现第 07 步西班牙在整个核心作业区以及 100°W 以东和以西的指数。
step07 <- subset(read.csv(file.path(out,"indices_baseline.csv")),id=="ESP_all" & domain=="ESP")
m <- merge(step07,series,by=c("region","threshold","year"),suffixes=c(".a",".r"))
stopifnot(nrow(m)==nrow(step07),nrow(step07)==3L*3L*length(years))
diffs <- sapply(c("index","lower","upper","cpue"),function(v) max(abs(m[[paste0(v,".a")]]-m[[paste0(v,".r")]])))
cat("Largest difference from step 07 / 与第 07 步的最大差值:", format(max(diffs), digits=2), "\n")
stopifnot(all(diffs<1e-6))

# Keep the four assessment areas, whose cells partition the core area at each threshold.
# 保留四个评估分区，其网格在每个阈值下构成核心作业区的划分。
areas <- subset(series,startsWith(region,"area"))
for (th in 2:4) stopifnot(sum(unique(areas[areas$threshold==th,c("region","domain_cells")])$domain_cells)==
  unique(series$domain_cells[series$region=="all" & series$threshold==th]))
write.csv(areas,file.path(out,"spain_area_indices.csv"),row.names=FALSE)
print(subset(areas,threshold==3 & year %in% c(2016,2020,2023))[c("region","domain_cells","year","index")],row.names=FALSE)

# ---- Compare with the archived results / 与存档结果比较 ----
source("analysis_walkthrough/common.R")
compare_table("spain_area_indices.csv")
