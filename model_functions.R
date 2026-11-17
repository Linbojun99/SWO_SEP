# =====================================================================================================
# Shared functions for fitting the fleet models to the public domain records
# 用于以公开记录拟合各船队模型的共用函数
#
# Sourced by steps 07, 08, 13, 16 and 17. fit_task() fits one row of tasks.csv: the spatiotemporal NB2 model of the
# retained swordfish count with year and month effects, the offset for the number of hooks, a spatial field
# and independent annual spatiotemporal fields (Eqs. 1 and 2), and, for the composition models, PC scores or
# clusters. It then calculates the bias-corrected annual indices over the core area or the areas of overlap.
# Optional parts write the parameter estimates, the randomized quantile residuals and the proportions of
# zeros (Section S3), and the residuals by cell and the relative CPUE in each cell (Sections S4 and S9).
# 由第 07、08、13、16 和 17 步读入。fit_task() 拟合 tasks.csv 的一行：剑鱼留存尾数的时空 NB2 模型，含年份效应、月份
# 效应、钩数偏置项、空间随机场和逐年独立的时空随机场（式 1 和式 2），组成模型另含 PC 得分或聚类。随后在核心作业区
# 或重叠区域上计算经偏差校正的年度指数。可选部分写出参数估计值、随机分位数残差和零值比例（补充材料 S3），以及按
# 网格汇总的残差和各网格的相对 CPUE（补充材料 S4 和 S9）。
#
# Reads / 读取:  nothing by itself / 本身不读取文件
# Writes / 写出: analysis_walkthrough/output/fits/<id>_fit.csv, <id>_series.csv and, when asked,
#                <id>_diagnostics.csv, <id>_parameters.csv, <id>_qq.csv, <flag>_residual_cells.csv,
#                <flag>_cpue_map.csv
# =====================================================================================================

suppressPackageStartupMessages(library(sdmTMB))
stopifnot(as.character(packageVersion("sdmTMB"))=="1.2.0")

# ---- Catch composition / 渔获组成 ----
# Square root of each species' proportion of the record's catch (Legendre and Gallagher, 2001).
# 每个鱼种占该记录渔获比例的平方根（Legendre and Gallagher, 2001）。
targeting_hellinger <- function(counts) {
  x <- as.matrix(counts); storage.mode(x) <- 'double'
  stopifnot(ncol(x)>=2,all(is.finite(x)),all(x>=0),all(rowSums(x)>0))
  sqrt(x/rowSums(x))
}
# PC scores of each record and the nearest cluster centre in the space of the retained scores.
# 每条记录的 PC 得分，以及在保留得分空间中最近的聚类中心。
targeting_apply <- function(model,counts) {
  h <- targeting_hellinger(as.matrix(counts)[,model$species,drop=FALSE])
  scores <- predict(model$pca,newdata=h)[,seq_len(model$n_pc),drop=FALSE]
  ds <- vapply(seq_len(model$k),function(j)rowSums(sweep(scores,2,model$centers[j,])^2),numeric(nrow(h)))
  ds <- matrix(ds,nrow=nrow(h),ncol=model$k)
  z <- as.data.frame(scores);names(z)<-paste0('PC',seq_len(ncol(z)))
  z$cluster <- factor(max.col(-ds,ties.method='first'),levels=seq_len(model$k));z
}
# Keep the fewest components that explain at least 80% of the variance (at most three), then choose two, three
# or four k-means clusters (20 random starts) by mean silhouette width on the same 1,000 records, among the
# solutions in which every cluster holds at least 2% of the records (Section S7).
# 保留累计解释至少 80% 方差的最少主成分（最多三个），再在每个聚类至少含 2% 记录的方案中，按同一批 1,000 条记录
# 的平均轮廓宽度，从两类、三类或四类 k 均值聚类（20 个随机起点）中选择（补充材料 S7）。
targeting_learn <- function(counts,seed=1,variance=.8,max_pc=3,k_candidates=2:4,silhouette_n=1000,min_share=.02) {
  h <- targeting_hellinger(counts)
  if(is.null(colnames(h)))colnames(h)<-paste0('species',seq_len(ncol(h)))
  if(nrow(h)<10||sum(apply(h,2,var))<1e-12)stop('Composition has insufficient variation or records')
  pca <- prcomp(h,center=TRUE,scale.=FALSE)
  ev <- pca$sdev^2/sum(pca$sdev^2)
  npc <- min(which(cumsum(ev)>=variance)[1],max_pc,sum(pca$sdev>1e-10))
  scores <- pca$x[,seq_len(npc),drop=FALSE]
  set.seed(seed);ii<-sort(sample.int(nrow(h),min(nrow(h),silhouette_n)))
  dist_sample <- dist(scores[ii,,drop=FALSE]);fits<-list();sel<-list()
  for(k in k_candidates) {
    set.seed(seed+k)
    km<-tryCatch(kmeans(scores,centers=k,nstart=20,iter.max=100),error=function(e)NULL)
    if(is.null(km))next
    share<-min(tabulate(km$cluster,nbins=k))/nrow(h)
    lab<-km$cluster[ii];sil<-if(length(unique(lab))>1)mean(cluster::silhouette(lab,dist_sample)[,'sil_width'])else NA_real_
    fits[[as.character(k)]]<-km;sel[[as.character(k)]]<-data.frame(k=k,silhouette=sil,min_share=share,eligible=is.finite(sil)&&share>=min_share)
  }
  selection<-do.call(rbind,sel)
  if(is.null(selection)||!any(selection$eligible))stop('No cluster candidate satisfies predeclared support rule')
  eligible<-selection[selection$eligible,];best<-eligible$k[which.max(eligible$silhouette)]
  km<-fits[[as.character(best)]];centers<-km$centers[order(km$centers[,1]),,drop=FALSE]
  list(species=colnames(h),pca=pca,n_pc=npc,k=best,centers=centers,
       selection=selection,explained_variance=ev,retained_variance=sum(ev[seq_len(npc)]),seed=seed)
}

# ---- Annual index / 年度指数 ----
# Bias-corrected, area-weighted annual index (epsilon estimator of Thorson and Kristensen, 2016) with its
# pointwise 95% interval exp(log(est) -/+ 1.96 se). The weights are cell areas that sum to one in each year.
# 经偏差校正、按面积加权的年度指数（Thorson and Kristensen, 2016 的 epsilon 估计）及其逐点 95% 区间
# exp(log(est) -/+ 1.96 se)。权重为网格面积，每年之和为 1。
bias_corrected_index <- function(fit, newdata, weights) {
  stopifnot(length(weights)==nrow(newdata), all(weights>0),
            all(abs(drop(rowsum(weights,newdata$t))-1)<1e-12))
  ix <- get_index(fit, newdata=newdata, offset=rep(0,nrow(newdata)), area=weights, bias_correct=TRUE)
  stopifnot(identical(as.integer(ix$t), sort(unique(as.integer(newdata$t)))),
            all(is.finite(ix$est)), all(ix$est>0), all(is.finite(ix$se)), all(ix$se>0),
            all(abs(ix$lwr/ix$est-exp(-qnorm(.975)*ix$se))<1e-10), all(abs(ix$upr/ix$est-exp(qnorm(.975)*ix$se))<1e-10))
  data.frame(t=as.integer(ix$t), est=ix$est, lower=ix$lwr, upper=ix$upr, se=ix$se)
}

# ---- Diagnostics (Section S3) / 诊断（补充材料 S3） ----
# Write the parameter estimates, AIC, the reported and fitted proportions of zeros and the randomized quantile
# residuals. The fixed parameters are held at their estimates and the random fields are drawn once from their
# approximate conditional distribution; the same draw gives the fitted probability of zero catch.
# 写出参数估计值、AIC、报告与拟合的零值比例以及随机分位数残差。固定参数取其估计值，随机场从其近似条件分布中
# 抽取一次；零渔获的拟合概率也由同一次抽样得到。
write_diagnostics <- function(fit, data, record, dir, seed) {
  parameters <- lapply(c("fixed", "ran_pars"), function(group) {
    estimates <- as.data.frame(generics::tidy(fit, effects=group, conf.int=TRUE))
    data.frame(id=record$id, group, estimates[c("term","estimate","conf.low","conf.high")])
  })
  write.csv(do.call(rbind, parameters), file.path(dir, paste0(record$id,"_parameters.csv")), row.names=FALSE)
  # Sample uniformly over the NB2 cumulative probability interval of each reported count; the upper-tail
  # form avoids rounding error for large counts.
  # 在每个报告尾数的 NB2 累积概率区间内均匀抽样；上尾写法可避免大尾数的舍入误差。
  zero_probability <- NULL
  quantile_residual <- function(object, y, mu, ...) {
    stopifnot(length(y)==nrow(data), identical(as.numeric(y),as.numeric(data$SWOn)),
              all(is.finite(mu)), all(mu>0))
    phi <- sigma(object)
    zero_probability <<- dnbinom(0, size=phi, mu=mu)
    upper <- y > mu
    a <- pnbinom(y-1, size=phi, mu=mu)
    b <- pnbinom(y, size=phi, mu=mu)
    a[upper] <- pnbinom(y[upper], size=phi, mu=mu[upper], lower.tail=FALSE)
    b[upper] <- pnbinom(y[upper]-1, size=phi, mu=mu[upper], lower.tail=FALSE)
    u <- runif(length(y), a, b)
    r <- qnorm(u)
    r[upper] <- qnorm(u[upper], lower.tail=FALSE)
    r
  }
  set.seed(seed)
  r <- residuals(fit, type="mle-mvn", qres_func=quantile_residual)
  stopifnot(length(r)==nrow(data), all(is.finite(r)))
  summary <- data.frame(record[c("id","flag","model","records","cells","optimizer","pdHess","gradient")],
    AIC=AIC(fit), observed_zero=mean(data$SWOn==0), expected_zero=mean(zero_probability),
    residual_mean=mean(r), residual_sd=sd(r), lower_tail=mean(r < -1.96),
    upper_tail=mean(r > 1.96), seed)
  write.csv(summary, file.path(dir, paste0(record$id,"_diagnostics.csv")), row.names=FALSE)
  prob <- c(.001,.002,.005,seq(.01,.99,.01),.995,.998,.999)
  qq <- data.frame(id=record$id, probability=prob, normal=qnorm(prob),
                   residual=unname(quantile(r, prob)))
  write.csv(qq, file.path(dir, paste0(record$id,"_qq.csv")), row.names=FALSE)
  invisible(r)
}

# ---- Residuals by cell and relative CPUE in each cell (Sections S4 and S9) / 按网格的残差与各网格相对 CPUE ----
# For a baseline model fitted to all eligible records, sum the residuals by cell and year, and predict CPUE
# in each cell of the fleet's area of overlap with Spain (for Spain, all five areas of overlap) for the
# earliest calendar month. Within each of 1,000 joint draws, each cell's annual prediction is divided by its
# mean over the 15 years; the median and 2.5th and 97.5th percentiles give the relative CPUE and its interval.
# 对用全部合格记录拟合的基线模型，按网格和年份汇总残差，并在该船队与西班牙的重叠区域（西班牙为全部五个重叠
# 区域）的每个网格上按最早的日历月份预测 CPUE。在 1,000 次联合抽样的每一次中，各网格的年度预测除以其 15 年
# 均值；中位数和第 2.5、97.5 百分位数给出相对 CPUE 及其区间。
write_spatial_outputs <- function(fit, d, residual, z, cfg, flag, task, dir) {
  stopifnot(z$tasks$model[task]=="all", length(residual)==nrow(d), all(is.finite(residual)))
  keys <- d[c("cell","LonC5","LatC5","Year")]
  by_cell <- aggregate(cbind(records=rep(1L,nrow(d)),sum_residual=residual),keys,sum)
  by_cell$flag <- flag
  write.csv(by_cell,file.path(dir,paste0(flag,"_residual_cells.csv")),row.names=FALSE)
  domains <- if (flag=="ESP") setdiff(cfg$flags,"ESP") else flag
  target <- with(z$coverage,cell[domain %in% domains & shared_years>=3])
  cells <- z$grid[z$grid$cell %in% target,]
  years <- sort(unique(d$Year))
  stopifnot(!anyDuplicated(cells$cell),length(z$refs)==15,all(z$refs %in% years))
  p <- cells[rep(seq_len(nrow(cells)),times=length(years)),]
  p$Year <- rep(years,each=nrow(cells))
  p$t <- match(p$Year,years)
  p$year_f <- factor(p$Year,levels=levels(d$year_f))
  p$month_f <- factor(min(d$Month),levels=levels(d$month_f))
  set.seed(cfg$seed+8000L+task)
  draws <- exp(predict(fit,newdata=p,type="link",offset=rep(0,nrow(p)),nsim=cfg$draws))
  stopifnot(all(dim(draws)==c(nrow(p),cfg$draws)),all(is.finite(draws)),all(draws>0))
  spatial <- lapply(cells$cell,function(cell) {
    idx <- which(p$cell==cell)
    reference <- colMeans(draws[idx[p$Year[idx] %in% z$refs],,drop=FALSE])
    ratios <- sweep(draws[idx,,drop=FALSE],2,reference,"/")
    stopifnot(max(abs(colMeans(ratios[p$Year[idx] %in% z$refs,,drop=FALSE])-1))<1e-10)
    qs <- t(apply(ratios,1,quantile,probs=c(.025,.5,.975),names=FALSE))
    colnames(qs) <- c("lower","ratio","upper")
    counts <- table(factor(d$Year[d$cell==cell],levels=years))
    data.frame(flag,year=years,p[idx,c("cell","LonC5","LatC5")],records=as.integer(counts),qs,row.names=NULL)
  })
  write.csv(do.call(rbind,spatial),file.path(dir,paste0(flag,"_cpue_map.csv")),row.names=FALSE)
}

# ---- One fit of tasks.csv / tasks.csv 中的一个拟合 ----
# Fit row i of tasks.csv and write its fit summary and annual indices to dir. The baseline model fitted to
# all eligible records uses every eligible record of the fleet; the other models use the records with at
# least one retained albacore, bigeye tuna, yellowfin tuna, blue marlin or striped marlin.
# 拟合 tasks.csv 的第 i 行，并把拟合汇总和年度指数写到 dir。用全部合格记录拟合的基线模型使用该船队的全部合格
# 记录；其他模型使用至少留存一尾长鳍金枪鱼、大眼金枪鱼、黄鳍金枪鱼、蓝枪鱼或条纹四鳍旗鱼的记录。
fit_task <- function(i, z, cfg, dir, diagnostics=FALSE, spatial=FALSE) {
  flag <- z$tasks$flag[i]; model <- z$tasks$model[i]; id <- z$tasks$id[i]
  fleet <- subset(z$d, Flag==flag)
  d <- if (model=="all") fleet else subset(fleet, common)
  # Year and month enter as factors; t numbers the years for the spatiotemporal fields.
  # 年份和月份作为因子进入模型；t 为时空随机场的年份编号。
  years <- sort(unique(d$Year)); months <- sort(unique(d$Month))
  d$year_f <- factor(d$Year,levels=years); d$month_f <- factor(d$Month,levels=months)
  d$t <- match(d$Year,years)
  rec <- data.frame(id,flag,model,annual_only=FALSE,records=nrow(d),cells=length(unique(d$cell)),
    version=as.character(packageVersion("sdmTMB")),
    optimizer=NA_integer_,pdHess=NA,gradient=NA_real_,k=NA_integer_,pcs=NA_integer_,
    status="error",detail="",stringsAsFactors=FALSE)
  run <- function() {
    stopifnot(all(z$refs %in% years))
    # For the composition models, add the PC scores or the cluster calculated from the five other species
    # (and swordfish when it is included in the composition); the reference cluster is cluster 1.
    # 组成模型加入由另外五个鱼种（组成中包含剑鱼时另加剑鱼）计算的 PC 得分或聚类；参照聚类为第 1 类。
    features <- character()
    if (grepl("cluster|pc",model)) {
      sp <- cfg$companion
      if (grepl("include",model)) sp <- c(sp,"SWOn")
      fm <- targeting_learn(d[sp],seed=cfg$seed+match(flag,cfg$flags))
      rec$k <<- fm$k; rec$pcs <<- fm$n_pc
      d <- cbind(d,targeting_apply(fm,d[sp]))
      features <- if (startsWith(model,"pc")) paste0("PC",seq_len(fm$n_pc)) else "cluster"
    }
    # The triangular mesh is built from the fleet's own cell centres, with coordinates in kilometres.
    # 三角网格由该船队自身的网格中心构建，坐标单位为千米。
    g <- subset(z$grid,cell %in% fleet$cell)
    geo <- fmesher::fm_mesh_2d(loc=as.matrix(g[c("x","y")]),cutoff=100,max.edge=c(850,1800),offset=c(650,1400))
    mesh <- make_mesh(d,c("x","y"),mesh=geo)
    # NB2 model with the offset for the number of hooks, log(hooks/1000), a spatial field and independent
    # annual spatiotemporal fields that share one Matern range.
    # NB2 模型，钩数偏置项为 log(钩数/1000)，空间随机场与逐年独立的时空随机场共用一个 Matern 范围参数。
    fit <- sdmTMB(reformulate(c("year_f","month_f",features),"SWOn"),data=d,mesh=mesh,time="t",
      family=nbinom2(),offset=log(d$Hooks/1000),spatial="on",spatiotemporal="iid",
      share_range=TRUE, control=sdmTMBcontrol(backend="tmb",newton_loops=2,get_joint_precision=TRUE))
    # A fit is accepted when the optimizer converged, the Hessian is positive definite and the maximum
    # absolute gradient is below 0.001.
    # 优化器收敛、Hessian 矩阵正定且最大绝对梯度小于 0.001 时，拟合被接受。
    rec$optimizer <<- fit$model$convergence; rec$pdHess <<- fit$sd_report$pdHess
    rec$gradient <<- max(abs(fit$tmb_obj$gr(fit$model$par)))
    if (!isTRUE(rec$optimizer==0 && rec$pdHess && rec$gradient<cfg$gradient_limit)) {
      rec$status <<- "fit_gate_failed"
      return()
    }
    # Predict for the earliest calendar month, at PC scores of zero or for the reference cluster. Spain is
    # predicted over its core area and all five areas of overlap, each other fleet over its own area of overlap.
    # 按最早的日历月份预测，PC 得分取零或取参照聚类。西班牙在核心作业区和全部五个重叠区域上预测，其他船队只在
    # 自身的重叠区域上预测。
    areas <- if (flag=="ESP") cfg$flags else flag
    masks <- subset(z$coverage,domain %in% areas & shared_years>=2)
    pg <- subset(z$grid,cell %in% masks$cell)
    p <- pg[rep(seq_len(nrow(pg)),times=length(years)),]
    p$Year <- rep(years,each=nrow(pg)); p$t <- match(p$Year,years)
    p$year_f <- factor(p$Year,levels=years); p$month_f <- factor(min(months),levels=months)
    if ("cluster" %in% features) p$cluster <- factor(1,levels=levels(d$cluster))
    for (k in grep("^PC",features,value=TRUE)) p[[k]] <- 0
    # One index for each area and for thresholds of two, three and four years; Spain's core area is also
    # split east and west of 100W.
    # 每个区域在两年、三年和四年阈值下各计算一个指数；西班牙核心作业区另分 100°W 以东和以西。
    series <- list()
    for (area in areas) for (threshold in 2:4) {
      target <- masks$cell[masks$domain==area & masks$shared_years>=threshold]
      regions <- if (flag=="ESP" && area=="ESP") c("all","east","west") else "all"
      for (region in regions) {
        keep <- p$cell %in% target
        if (region=="east") keep <- keep & p$LonC5>cfg$east_split
        if (region=="west") keep <- keep & p$LonC5<cfg$east_split
        sub <- p[keep,]; target_cells <- unique(sub$cell)
        # Weight each cell by the spherical area of the full 5-degree cell.
        # 每个网格按完整 5° 网格的球面面积加权。
        w <- sub$area/sum(pg$area[pg$cell %in% target_cells])
        # A failed bias correction makes the whole fit not accepted.
        # 偏差校正失败时，整个拟合不被接受。
        bc <- tryCatch(bias_corrected_index(fit,sub,w),error=function(e) e)
        if (inherits(bc,"error")) {
          rec$status <<- "bias_correction_failed"
          rec$detail <<- paste(area,region,threshold,conditionMessage(bc))
          return()
        }
        stopifnot(identical(bc$t,seq_along(years)))
        # Divide the index and both interval limits by the 15-year mean of the point estimates.
        # 指数及区间上下限都除以点估计的 15 年均值。
        pt <- bc$est
        ref <- years %in% z$refs
        idx <- pt/mean(pt[ref])
        ci <- cbind(bc$lower,bc$upper)/mean(pt[ref])
        stopifnot(abs(mean(idx[ref])-1)<1e-10,all(ci[,1]<idx),all(idx<ci[,2]),
                  max(abs((log(ci[,2])-log(ci[,1]))/(2*qnorm(.975))-bc$se))<1e-10)
        # Count the records and the cells with records in each year; open symbols in Fig. 3 mark years
        # without records in the cells indexed.
        # 统计每年的记录数和有记录的网格数；图 3 的空心符号标出所计算网格内没有记录的年份。
        observed <- d[d$cell %in% target_cells,]
        count <- table(factor(observed$Year,levels=years))
        obs_cells <- unique(observed[c("cell","Year")])
        cell_count <- table(factor(obs_cells$Year,levels=years))
        series[[length(series)+1L]] <- data.frame(id,flag,model,annual_only=FALSE,domain=area,region,threshold,
          domain_cells=length(target_cells),year=years,records=as.integer(count),observed_cells=as.integer(cell_count),
          index=idx,lower=ci[,1],upper=ci[,2],cpue=pt,row.names=NULL)
      }
    }
    # Diagnostics use the random seed cfg$seed + 4000 + row number of tasks.csv.
    # 诊断所用的随机种子为 cfg$seed + 4000 + tasks.csv 的行号。
    if (diagnostics) {
      r <- write_diagnostics(fit,d,rec,dir,cfg$seed+4000L+i)
      if (spatial) write_spatial_outputs(fit,d,r,z,cfg,flag,i,dir)
    }
    write.csv(do.call(rbind,series),file.path(dir,paste0(id,"_series.csv")),row.names=FALSE)
    rec$status <<- "accepted"
  }
  tryCatch(run(),error=function(e) rec$detail <<- conditionMessage(e))
  write.csv(rec,file.path(dir,paste0(id,"_fit.csv")),row.names=FALSE)
  # No archived fit needed the refit without the spatial field, so the walkthrough stops instead.
  # 没有存档拟合需要去掉空间随机场重新拟合，因此本流程在此停止。
  if (rec$status!="accepted")
    stop(id, " was not accepted (", rec$status, " ", rec$detail, "); the production analysis would refit it ",
         "without the spatial field (R/fit_directed_fishery.R, argument annual).")
  rec
}
