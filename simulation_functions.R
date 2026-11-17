# =====================================================================================================
# Shared set-up of the simulation of swordfish catches with known abundance (Section 2.4)
# 已知丰度的剑鱼渔获模拟的共用设定（第 2.4 节）
#
# Sourced by steps 20, 21, 24, 26 and 29. Defines the 141 simulated cells (92 in Spain's core area and 49 outside
# it), the numbers of hooks of each fleet by cell and year in the eight years 2009-2023 at two-year intervals,
# the three nested sampling designs of Table 2, the abundance scenarios of Table 1 and experiment(), which
# simulates one replicate and fits the estimation model to each sampling design. These are simulation
# results, not estimates from the fishery data.
# 由第 20、21、24、26 和 29 步读入。定义 141 个模拟网格（西班牙核心作业区 92 个、核心作业区以外 49 个）、2009-2023 年
# 每隔两年共八个年份各船队按网格和年份的钩数、表 2 的三种嵌套取样设计、表 1 的丰度情景，以及 experiment()：模拟
# 一次重复，并以估计模型拟合每种取样设计。这些是模拟结果，不是渔业数据的估计。
#
# Reads / 读取:  analysis_walkthrough/output/input.rds (step 04), config/directed_fishery.json,
#                data/raw/longline/PublicLLTunaBillfishNum.csv
# Writes / 写出: analysis_walkthrough/output/replicates/<run>/rep_<replicate>.rds (one file per replicate,
#                so an interrupted run resumes / 每次重复一个文件，中断后可接着运行)
# =====================================================================================================

suppressPackageStartupMessages(library(sdmTMB))
stopifnot(as.character(packageVersion("sdmTMB"))=="1.2.0")
out <- "analysis_walkthrough/output"
cfg <- jsonlite::fromJSON("config/directed_fishery.json")
z <- readRDS(file.path(out, "input.rds"))

# ---- Number of replicates and workers / 重复次数与并行进程数 ----
# The number of replicates is read from the config (100); a smaller number can be given for a quick run with
# --replicates N. Replicates run in parallel on --workers N processes (default 8). Each replicate sets its own
# random seed, so the results do not depend on the number of workers.
# 重复次数从配置文件读取（100）；快速运行时可用 --replicates N 指定较小的次数。各重复在 --workers N 个进程上并行
# 运行（默认 8 个）。每次重复设定自己的随机种子，因此结果与进程数无关。
argument <- function(name, default) {
  argv <- commandArgs(trailingOnly=TRUE)
  at <- match(name, argv)
  if (is.na(at)) default else as.integer(argv[at + 1L])
}
n_rep <- argument("--replicates", as.integer(cfg$coverage_replicates))
n_workers <- argument("--workers", 8L)
stopifnot(length(n_rep)==1L, n_rep>=1L, n_rep<=cfg$coverage_replicates, n_workers>=1L)
cat(sprintf("Replicates / 重复次数: %d of %d | workers / 并行进程: %d\n", n_rep, cfg$coverage_replicates, n_workers))

# ---- Simulated cells / 模拟网格 ----
# Keep the cells in which at least one fleet reported effort in three or more years from 2005 to 2023,
# counted from all public domain records in H2 with positive hooks and no catch in weight.
# 保留 2005-2023 年间至少有一支船队在三年或以上报告努力量的网格；按 H2 内钩数为正、无重量渔获的全部公开记录统计。
raw <- read.csv(cfg$num)
raw <- subset(raw,Flag %in% cfg$flags & Year>=min(cfg$years) & Year<=cfg$support_end &
              Hooks>0 & LatC5<cfg$north_boundary & LonC5>cfg$west_boundary)
raw <- raw[rowSums(raw[grepl("mt$",names(raw))]>0)==0,]
raw$cell <- paste(raw$LonC5,raw$LatC5,sep="_")
support <- aggregate(Year~cell+Flag,unique(raw[c("Year","Flag","cell")]),length)
g <- subset(z$grid,cell %in% support$cell[support$Year>=cfg$minimum_distinct_years])
# Mark the cells of Spain's core area (Spanish effort in at least three years).
# 标出西班牙核心作业区的网格（西班牙至少三年有努力量）。
g$spanish <- g$cell %in% subset(z$coverage,domain=="ESP" & shared_years>=3)$cell
years <- seq(2009,2023,by=2); nt <- length(years); n <- nrow(g)
stopifnot(n==141L, sum(g$spanish)==92L)

# ---- Reported hooks / 报告的钩数 ----
# Sum the reported hooks of the eligible records by cell, year and fleet. One simulated record is one fleet
# in one cell and year.
# 按网格、年份和船队汇总合格记录报告的钩数。一条模拟记录代表一支船队在某年某网格的作业。
d <- aggregate(Hooks~cell+Year+Flag,subset(z$d,Year %in% years & cell %in% g$cell),sum)
d <- merge(d,g,by="cell")
d$t <- match(d$Year,years)
d$year_f <- factor(d$Year)
d$fleet_f <- factor(d$Flag,levels=cfg$flags)

# ---- Sampling designs (Table 2) / 取样设计（表 2） ----
# Spain only; Spanish records plus the tuna-targeting longline fleets' records from the core area (+ core area);
# Spanish records plus their records from both inside and outside the core area (+ outside core area).
# 仅西班牙；西班牙记录加金枪鱼延绳钓船队在核心作业区的记录（+ 核心作业区）；西班牙记录加这些船队在核心作业区
# 内外的记录（+ 核心作业区以外）。
design_rows <- list(Spain=which(d$Flag=="ESP"),Shared=which(d$Flag=="ESP" | d$spanish),
                    Extended=seq_len(nrow(d)))

# ---- Prediction grid, fields and mesh / 预测网格、随机场与三角网格 ----
# Predict every cell in every year at Spain's catchability.
# 以西班牙的可捕性预测每年的每个网格。
p <- g[rep(seq_len(n),nt),]
p$t <- rep(seq_len(nt),each=n)
p$year_f <- factor(years[p$t],levels=years)
p$fleet_f <- factor("ESP",levels=cfg$flags)
# The generating fields have an exponential correlation function with correlation 0.37 at 900 km; the
# estimation model uses one triangular mesh for all designs.
# 生成渔获的随机场采用指数相关函数，900 千米处相关为 0.37；估计模型的所有设计使用同一个三角网格。
xy <- as.matrix(g[c("x","y")])
L <- t(chol(exp(-as.matrix(dist(xy))/900)))
geo <- fmesher::fm_mesh_2d(loc=xy,cutoff=500,max.edge=c(1100,2000),offset=c(700,1400))
# Catchability of each fleet relative to Spain's (Table S7).
# 各船队相对于西班牙的可捕性（表 S7）。
q <- c(ESP=1,JPN=.65,KOR=.50,CHN=.70,TWN=.65,VUT=.55)
# Change over the study period: the mean of 2021 and 2023 divided by the mean of 2009 and 2011.
# 研究期内的变化：2021 和 2023 年的均值除以 2009 和 2011 年的均值。
ratio <- function(a) mean(tail(a,2))/mean(head(a,2))

# ---- Abundance scenarios (Table 1) / 丰度情景（表 1） ----
# Multipliers of mean density outside the core area: constant (A1), a log-linear decline from 1 in 2009 to 0.5
# in 2023 (A2), and one cosine cycle from 1.2 in 2009 to 0.8 in 2016 and back to 1.2 in 2023 (A3), after the
# fluctuating population of Forrestal et al. (2019).
# 核心作业区以外平均密度的乘数：恒定（A1）、由 2009 年的 1 对数线性下降到 2023 年的 0.5（A2），以及一个余弦周期，
# 由 2009 年的 1.2 降到 2016 年的 0.8 再回到 2023 年的 1.2（A3），参照 Forrestal et al.（2019）的波动种群。
scenarios <- c("stable","periphery_decline","periphery_fluctuation")
multiplier <- list(stable=rep(1,nt),periphery_decline=.5^seq(0,1,length.out=nt),
                   periphery_fluctuation=1+.2*cos(2*pi*(years-min(years))/(max(years)-min(years))))
# Weight each cell by its area divided by the summed area of the 141 cells.
# 每个网格按其面积除以 141 个网格的总面积加权。
w <- p$area/sum(g$area)
stopifnot(all(abs(tapply(w,p$t,sum)-1)<1e-12))

# ---- Cache of replicates / 重复结果的缓存 ----
# Run fun(i) for each replicate in parallel and keep each result as one file, so an interrupted run resumes.
# Delete the folder after any change to the code, the config or the input data.
# 并行对每次重复运行 fun(i)，每个结果保存为一个文件，中断后可接着运行。代码、配置或输入数据改动后需删除该文件夹。
cached_replicates <- function(tag, ids, fun) {
  dir <- file.path(out,"replicates",tag)
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  path <- function(i) file.path(dir,sprintf("rep_%03d.rds",i))
  todo <- ids[!file.exists(path(ids))]
  cat(sprintf("[%s] %d of %d replicates already done, %d to run / 已完成 %d 次，待运行 %d 次\n",
              tag, length(ids)-length(todo), length(ids), length(todo), length(ids)-length(todo), length(todo)))
  if (length(todo)) {
    done <- parallel::mclapply(todo,function(i) {
      ans <- fun(i)
      part <- paste0(path(i),".part")
      saveRDS(ans,part)
      file.rename(part,path(i))
    },mc.cores=n_workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
    failed <- todo[!vapply(done,isTRUE,logical(1))]
    if (length(failed)) stop(sprintf("[%s] replicates failed: %s",tag,paste(failed,collapse=", ")))
  }
  lapply(ids,function(i) readRDS(path(i)))
}

# ---- One replicate / 一次重复 ----
# Simulate density in every cell and year, generate counts from the reported hooks and the fleets'
# catchability, and fit the estimation model to each sampling design. The scenarios of one replicate share
# the same random fields and the same uniform draws for the counts. q_trend is the annual increase in the
# catchability of the five tuna-targeting longline fleets (scenarios Q1 to Q6); 0 leaves it constant.
# 模拟每个网格每年的密度，由报告的钩数和各船队的可捕性生成尾数，并以估计模型拟合每种取样设计。同一次重复的各
# 情景共用相同的随机场和生成尾数的均匀随机数。q_trend 为五支金枪鱼延绳钓船队可捕性的年增幅（情景 Q1 至 Q6）；
# 为 0 时可捕性恒定。
experiment <- function(rep_id, run=scenarios, sampling_rows=design_rows, q_trend=0) {
  # One random seed per replicate: cfg$seed + 2000 + replicate number.
  # 每次重复一个随机种子：cfg$seed + 2000 + 重复编号。
  set.seed(cfg$seed+2000+rep_id)
  # Spatial field with SD 0.4 on the log scale, scaled to an area-weighted mean of one.
  # 对数尺度标准差 0.4 的空间随机场，缩放到面积加权均值为 1。
  baseline <- exp(.4*drop(L%*%rnorm(n)))
  baseline <- baseline/weighted.mean(baseline,g$area)
  # Independent annual spatiotemporal fields with SD 0.22 on the log scale.
  # 对数尺度标准差 0.22 的逐年独立时空随机场。
  latent <- exp(.22*(L%*%matrix(rnorm(n*nt),n))) * baseline
  # Rescale the density in each year so that the area mean inside and outside the core area keeps its initial value.
  # 逐年重新缩放密度，使核心作业区内外的面积均值保持初始值。
  for (s in c(TRUE,FALSE)) {
    mask <- g$spanish==s
    marginal <- colSums(latent[mask,,drop=FALSE]*g$area[mask])/sum(g$area[mask])
    latent[mask,] <- sweep(latent[mask,,drop=FALSE],2,marginal,"/")*weighted.mean(baseline[mask],g$area[mask])
  }
  # Uniform draws that turn expected catches into NB2 counts, shared by all scenarios.
  # 把期望渔获转为 NB2 尾数的均匀随机数，各情景共用。
  noise <- runif(nrow(d))
  rows <- list()
  paths <- list()
  counts <- list()
  area_truth <- list()
  # The abundance lost outside the core area in A2 is placed in the core area (L1) or spread over both areas
  # by the same proportion (L2), so that regional abundance is the same in every year (Section S11).
  # A2 中核心作业区以外损失的丰度在 L1 中改为发生在核心作业区，在 L2 中按相同比例分摊到两个区域，因此各年的区域
  # 丰度相同（补充材料 S11）。
  initial_core <- sum(baseline[g$spanish] * g$area[g$spanish])
  initial_outside <- sum(baseline[!g$spanish] * g$area[!g$spanish])
  loss <- initial_outside * (1 - multiplier$periphery_decline)
  for (scenario in run) {
    dens <- latent
    if (scenario=="core_decline") {
      core_factor <- 1 - loss / initial_core
      stopifnot(all(core_factor > 0))
      dens[g$spanish,] <- sweep(dens[g$spanish,,drop=FALSE],2,core_factor,"*")
    } else if (scenario=="uniform_decline") {
      dens <- sweep(dens,2,1 - loss / (initial_core + initial_outside),"*")
    } else if (scenario!="stable") {
      stopifnot(scenario %in% names(multiplier))
      dens[!g$spanish,] <- sweep(dens[!g$spanish,,drop=FALSE],2,multiplier[[scenario]],"*")
    }
    # True regional abundance is the area-weighted mean density over the 141 cells.
    # 真实区域丰度为 141 个网格上按面积加权的平均密度。
    known <- colSums(dens*g$area)
    if (scenario %in% c("periphery_decline","core_decline","uniform_decline")) {
      target <- initial_core + initial_outside - loss
      stopifnot(max(abs(known / target - 1)) < 1e-12)
      area_truth[[scenario]] <- data.frame(rep_id,scenario,year=years,regional=known,
        core=colSums(dens[g$spanish,,drop=FALSE] * g$area[g$spanish]),
        outside=colSums(dens[!g$spanish,,drop=FALSE] * g$area[!g$spanish]))
    }
    truth <- ratio(known/sum(g$area))
    # Expected catch = density x catchability x hooks/1000; the tuna-targeting longline fleets' catchability is
    # multiplied by (1 + q_trend)^(year - 2009). Counts follow the NB2 distribution with dispersion 5.
    # 期望渔获 = 密度 × 可捕性 × 钩数/1000；金枪鱼延绳钓船队的可捕性乘以 (1 + q_trend)^(年份 - 2009)。尾数服从
    # 离散参数为 5 的 NB2 分布。
    mu <- dens[cbind(match(d$cell,g$cell),d$t)] * q[d$Flag] * d$Hooks/1000
    if (q_trend!=0) mu <- mu * ifelse(d$Flag=="ESP",1,(1+q_trend)^(d$Year-min(years)))
    d$count <- qnbinom(noise,size=5,mu=mu)
    counts[[scenario]] <- d$count
    for (design in names(sampling_rows)) {
      train <- d[sampling_rows[[design]],]
      rec <- data.frame(rep_id,scenario,design,truth,records=nrow(train),target_cells=n,
        outside_spanish=sum(!g$spanish),status="error",estimate=NA_real_,gradient=NA_real_,detail="",
        uncorrected=NA_real_)
      tryCatch({
        # Estimation model: year effects, a spatial field and independent annual spatiotemporal fields,
        # the offset for the number of hooks and, for records from more than one fleet, a constant fleet effect.
        # 估计模型：年份效应、空间随机场和逐年独立的时空随机场、钩数偏置项；记录来自多支船队时另含恒定的船队效应。
        mesh <- make_mesh(train,c("x","y"),mesh=geo)
        form <- if (design=="Spain") count~year_f else count~year_f+fleet_f
        fit <- sdmTMB(form,data=train,time="t",mesh=mesh,family=nbinom2(),offset=log(train$Hooks/1000),
          spatial="on",spatiotemporal="iid",share_range=TRUE,
          control=sdmTMBcontrol(backend="tmb",newton_loops=2))
        rec$gradient <- max(abs(fit$tmb_obj$gr(fit$model$par)))
        # Accepted fits: converged, positive definite Hessian, maximum absolute gradient below 0.001,
        # and finite indices and standard errors.
        # 被接受的拟合：收敛、Hessian 矩阵正定、最大绝对梯度小于 0.001，且指数和标准误有限。
        if (fit$model$convergence==0 && isTRUE(fit$sd_report$pdHess) && rec$gradient<cfg$gradient_limit) {
          # Index without bias correction, kept only as a check of the earlier runs.
          # 未经偏差校正的指数，只作为与早期运行的核对而保留。
          plain <- predict(fit,newdata=p,offset=rep(0,nrow(p)),type="response")$est
          rec$uncorrected <- ratio(drop(rowsum(plain*p$area,p$t)))
          # Bias-corrected regional index at Spain's catchability, as in Section 2.2.
          # 以西班牙可捕性计算经偏差校正的区域指数，同第 2.2 节。
          idx <- tryCatch({
            pr <- predict(fit,newdata=p,offset=rep(0,nrow(p)),return_tmb_object=TRUE)
            ix <- get_index(pr,area=w,bias_correct=TRUE)
            stopifnot(nrow(ix)==nt,identical(as.integer(ix$t),seq_len(nt)),
                      all(is.finite(ix$est)),all(ix$est>0),all(is.finite(ix$se)))
            ix
          },error=function(e) e)
          if (inherits(idx,"error")) {
            rec$status <- "bias_correction_failed"
            rec$detail <- conditionMessage(idx)
          } else {
            # Each estimated and true series is scaled to its early mean (2009 and 2011).
            # 每条估计序列和真实序列都除以其早期均值（2009 和 2011 年）。
            rec$estimate <- ratio(idx$est)
            stopifnot(is.finite(rec$estimate),rec$estimate>0)
            ref <- mean(head(idx$est,2))
            paths[[length(paths)+1L]] <- data.frame(rep_id,scenario,design,year=years,
              estimated_index=idx$est/ref,true_index=known/mean(head(known,2)),
              lower=exp(idx$log_est-1.96*idx$se)/ref,upper=exp(idx$log_est+1.96*idx$se)/ref,
              trend=rec$estimate,true_trend=truth,trend_change_lwr=NA_real_,trend_change_upr=NA_real_)
            rec$status <- "accepted"
          }
        } else rec$status <- "fit_gate_failed"
      },error=function(e) rec$detail <<- conditionMessage(e))
      rows[[length(rows)+1]] <- rec
    }
  }
  list(fits=do.call(rbind,rows),annual=do.call(rbind,paths),area_truth=do.call(rbind,area_truth),
       counts=counts)
}

# ---- Two designs with the same number of records (Section S8) / 记录数相同的两种设计（补充材料 S8） ----
# In random order within each fleet and year, replace each core-area record of a tuna-targeting longline fleet
# with an unused record from outside the core area for the same fleet and year whose number of hooks differs
# by at most a factor of 1.05; keep the record when there is none. Samples are drawn in replicate order from
# one random stream, so replicates 1 to k are the same for any number of replicates of at least k.
# 在每支船队和年份内按随机顺序，把金枪鱼延绳钓船队在核心作业区的每条记录替换为同一船队同一年份、核心作业区以外
# 尚未使用且钩数相差不超过 1.05 倍的记录；没有这样的记录时保留原记录。各重复的样本按重复顺序从同一个随机数流中
# 抽取，因此只要重复次数不少于 k，第 1 至 k 次重复的样本都相同。
equal_record_samples <- function(d, seed, n=n_rep) {
  core <- which(d$Flag == "ESP" | d$spanish)
  bycatch <- which(d$Flag != "ESP")
  strata <- split(bycatch, interaction(d$Flag[bycatch], d$Year[bycatch], drop=TRUE))
  set.seed(seed)
  lapply(seq_len(n), function(rep_id) {
    replacement <- seq_len(nrow(d))
    for (stratum in strata) {
      inside <- stratum[d$spanish[stratum]]
      outside <- stratum[!d$spanish[stratum]]
      for (i in inside[sample.int(length(inside))]) {
        ratio <- d$Hooks[outside] / d$Hooks[i]
        candidates <- outside[ratio >= 1 / 1.05 & ratio <= 1.05]
        if (length(candidates) == 0L) next
        chosen <- candidates[sample.int(length(candidates), 1L)]
        replacement[i] <- chosen
        outside <- outside[outside != chosen]
      }
    }
    mixed <- replacement[core]
    stopifnot(!anyDuplicated(mixed), all(d$Flag[core] == d$Flag[mixed]),
              all(d$Year[core] == d$Year[mixed]),
              identical(core[d$Flag[core]=="ESP"], mixed[d$Flag[mixed]=="ESP"]),
              all(pmax(d$Hooks[core], d$Hooks[mixed]) <= 1.05 * pmin(d$Hooks[core], d$Hooks[mixed])))
    list(Core=core, Mixed=mixed)
  })
}
