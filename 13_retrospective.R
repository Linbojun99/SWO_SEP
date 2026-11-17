# =====================================================================================================
# Step 13. Retrospective analysis of the fleet indices (Section 2.2, Section S6, Fig. S4, Table S4)
# 第 13 步  各船队指数的回溯分析（第 2.2 节，补充材料 S6，图 S4，表 S4）
#
# For each fleet, omits the last one to five calendar years and refits the baseline model, with the same
# prediction area (the fleet's own cells at the three-year threshold), mesh and calendar month for prediction.
# Each truncated time series with an estimate in its terminal year gives the proportional difference between
# its terminal-year CPUE (swordfish per 1,000 hooks, before index scaling) and the estimate for the same year
# from the full time series. Mohn's rho is the mean of these differences (Mohn, 1999; Hurtado-Ferro et al., 2015).
# 对每支船队，去掉最后一至五个日历年并重新拟合基线模型，预测区域（该船队自身在三年阈值下的网格）、三角网格和预测
# 所用的日历月份保持不变。每条在末年有估计的截短时间序列，给出其末年 CPUE（每千钩剑鱼尾数，指数缩放之前）与完整
# 时间序列同年估计之间的比例差。Mohn's ρ 为这些比例差的均值（Mohn, 1999；Hurtado-Ferro et al., 2015）。
#
# Reads / 读取:  analysis_walkthrough/output/input.rds (step 04), indices_baseline.csv (step 07),
#                config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/retrospective_fits.csv, retrospective_indices.csv,
#                retrospective_rho.csv (Table S4 / 表 S4), FigS4_retrospective.pdf/.png
#
# Run from the project folder; the first line makes R use one core per process and read the Chinese folder
# names. The 36 fits run on --workers N processes (default 5).
# 在项目文件夹中运行；第一行让 R 每个进程只用一个核心，并能读取中文文件夹名。36 个拟合在 --workers N 个进程上运行
# （默认 5 个）。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/13_retrospective.R --workers 5
# =====================================================================================================

source("analysis_walkthrough/model_functions.R")
out <- "analysis_walkthrough/output"
cat("Step 13: retrospective analysis (Fig. S4, Table S4) / 第 13 步：回溯分析（图 S4、表 S4）\n")
cfg <- jsonlite::fromJSON("config/directed_fishery.json")
z <- readRDS(file.path(out, "input.rds"))
argv <- commandArgs(trailingOnly=TRUE)
workers <- if ("--workers" %in% argv) as.integer(argv[match("--workers", argv) + 1L]) else 5L

# Unscaled CPUE of each fleet's baseline index over its own cells at the three-year threshold (step 07).
# 第 07 步中各船队基线指数在其自身三年阈值网格上的未缩放 CPUE。
indices <- read.csv(file.path(out, "indices_baseline.csv"))
baseline <- subset(indices, region=="all" & domain==flag & threshold==3)

# Fit one truncated time series: peel 0 is the full time series and must reproduce step 07.
# 拟合一条截短时间序列：peel 0 为完整时间序列，必须重现第 07 步。
fit_peel <- function(flag,peel) {
  fleet <- subset(z$d,Flag==flag)
  cutoff <- max(fleet$Year)-peel
  d <- fleet[fleet$Year<=cutoff,]
  audit <- data.frame(flag,peel,cutoff,status="missing_terminal_year",n=nrow(d),optimizer=NA_integer_,
                      pdHess=NA,gradient=NA_real_,message="")
  # A truncated series without records in its terminal year has no terminal-year estimate.
  # 末年没有记录的截短序列没有末年估计。
  if (!cutoff %in% d$Year) return(list(audit=audit))
  years <- sort(unique(d$Year))
  month <- min(fleet$Month)
  g <- subset(z$grid,cell %in% fleet$cell)
  pg <- subset(z$grid,cell %in% z$coverage$cell[z$coverage$domain==flag & z$coverage$shared_years>=3])
  d$year_f <- factor(d$Year); d$month_f <- factor(d$Month); d$t <- match(d$Year,years)
  stopifnot(month %in% d$Month)
  geo <- fmesher::fm_mesh_2d(loc=as.matrix(g[c("x","y")]),cutoff=100,max.edge=c(850,1800),offset=c(650,1400))
  mesh <- make_mesh(d,c("x","y"),mesh=geo)
  fit <- sdmTMB(SWOn~year_f+month_f,data=d,mesh=mesh,time="t",family=nbinom2(),offset=log(d$Hooks/1000),
    spatial="on",spatiotemporal="iid",share_range=TRUE,
    control=sdmTMBcontrol(backend="tmb",newton_loops=2,get_joint_precision=TRUE))
  audit$optimizer <- fit$model$convergence; audit$pdHess <- fit$sd_report$pdHess
  audit$gradient <- max(abs(fit$tmb_obj$gr(fit$model$par)))
  if (!isTRUE(audit$optimizer==0 && audit$pdHess && audit$gradient<cfg$gradient_limit)) {
    audit$status <- "fit_gate_failed"
    return(list(audit=audit))
  }
  p <- pg[rep(seq_len(nrow(pg)),length(years)),]
  p$t <- rep(seq_along(years),each=nrow(pg))
  p$year_f <- factor(years[p$t],levels=years); p$month_f <- factor(month,levels=levels(d$month_f))
  bc <- tryCatch(bias_corrected_index(fit,p,p$area/sum(pg$area)),error=function(e) e)
  if (inherits(bc,"error")) {
    audit$status <- "bias_correction_failed"; audit$message <- conditionMessage(bc)
    return(list(audit=audit))
  }
  stopifnot(identical(bc$t,seq_along(years)))
  point <- bc$est
  ci <- cbind(bc$lower,bc$upper)
  reference <- baseline[baseline$flag==flag,]
  full <- reference$cpue[match(years,reference$year)]
  stopifnot(all(is.finite(full)),all(point>0))
  if (peel==0) stopifnot(max(abs(point/full-1))<1e-6)
  audit$status <- "accepted"
  list(audit=audit,series=data.frame(flag,peel,cutoff,year=years,cpue=point,lower=ci[,1],upper=ci[,2],
                                   full_cpue=full,relative_revision=point/full-1))
}

# Fit the six fleets with peels 0-5; the fits do not use random numbers, so the number of workers does not matter.
# 拟合六支船队的 peel 0-5；这些拟合不使用随机数，因此结果与进程数无关。
tasks <- expand.grid(flag=cfg$flags,peel=0:5,stringsAsFactors=FALSE)
ans <- parallel::mclapply(seq_len(nrow(tasks)),function(i) fit_peel(tasks$flag[i],tasks$peel[i]),
                          mc.cores=workers,mc.preschedule=FALSE)
stopifnot(all(vapply(ans,function(x)is.list(x)&&is.data.frame(x$audit),logical(1))))
audit <- do.call(rbind,lapply(ans,"[[","audit"))
series <- do.call(rbind,lapply(ans,"[[","series"))

# Mohn's rho: mean over the truncated series of the terminal-year estimate divided by the full-series estimate, minus 1.
# Mohn's ρ：各截短序列末年估计除以完整序列同年估计再减 1 的均值。
terminal <- subset(series,peel>0 & year==cutoff)
rho <- aggregate(relative_revision~flag,terminal,mean)
names(rho)[2] <- "mohn_rho"
rho$accepted_peels <- as.integer(table(terminal$flag)[rho$flag])
print(rho, digits=3, row.names=FALSE)
write.csv(audit,file.path(out,"retrospective_fits.csv"),row.names=FALSE)
write.csv(series,file.path(out,"retrospective_indices.csv"),row.names=FALSE)
write.csv(rho,file.path(out,"retrospective_rho.csv"),row.names=FALSE)

# ---- Fig. S4 / 图 S4 ----
# Full time series in black with its interval from 2015; truncated series in green with the terminal-year estimate and interval.
# 自 2015 年起绘制完整时间序列（黑色）及其区间；截短序列为绿色，并标出末年估计及其区间。
suppressPackageStartupMessages(library(ggplot2))
npg_green <- "#00A087"
theme_set(theme_bw(base_size=9,base_family="Arial") +
  theme(panel.grid=element_blank(),panel.border=element_rect(colour="grey25",linewidth=.50),
        strip.background=element_blank(),strip.text=element_text(hjust=0,face="bold"),
        legend.position="bottom",legend.title=element_blank(),legend.key.width=grid::unit(1.8,"lines")))
fleets <- c(ESP="Spain",JPN="Japan",KOR="Korea",CHN="China",TWN="Chinese Taipei",VUT="Vanuatu")
retro <- read.csv(file.path(out,"retrospective_indices.csv"))
rho <- read.csv(file.path(out,"retrospective_rho.csv"))
labs <- sprintf("%s   (ρ = %+.3f)",fleets,rho$mohn_rho[match(names(fleets),rho$flag)])
retro$panel <- factor(labs[match(retro$flag,names(fleets))],levels=labs)
retro <- subset(retro,year>=2015)
retro <- retro[order(retro$flag,retro$peel,retro$year),]
retro$group <- interaction(retro$flag,retro$peel,ave(retro$year,retro$flag,retro$peel,
                                                     FUN=function(y)cumsum(c(1,diff(y)!=1))))
full <- subset(retro,peel==0)
early <- subset(retro,peel>0)
endpoints <- subset(early,year==cutoff)
p <- ggplot(retro,aes(year,cpue,group=group)) + facet_wrap(~panel,ncol=2,scales="free_y") +
  geom_ribbon(data=full,aes(ymin=lower,ymax=upper),fill="black",alpha=.16,colour=NA) +
  geom_line(data=early,aes(colour="Truncated time series"),linewidth=.62,alpha=.82) +
  geom_line(data=full,aes(colour="Full time series"),linewidth=.82) +
  geom_linerange(data=endpoints,aes(ymin=lower,ymax=upper),colour=npg_green,linewidth=.70) +
  geom_point(data=endpoints,colour=npg_green,size=1.65) +
  scale_colour_manual(name="Series",values=c("Full time series"="black","Truncated time series"=npg_green)) +
  scale_x_continuous(breaks=unique(c(seq(min(retro$year),max(retro$year),by=3),max(retro$year)))) +
  labs(x="Year",y="Standardized swordfish catch per 1,000 hooks")
for (ext in c("pdf","png"))
  ggsave(file.path(out,paste0("FigS4_retrospective.",ext)),
         p + theme(legend.box.spacing=grid::unit(1.5,"mm"),legend.margin=margin(0,0,0,0)),width=190,height=155,
         units="mm",device=if (ext=="pdf") cairo_pdf else NULL,dpi=600,bg="white")

# ---- Compare with the archived results / 与存档结果比较 ----
source("analysis_walkthrough/common.R")
for (name in c("retrospective_fits.csv","retrospective_indices.csv","retrospective_rho.csv")) compare_table(name)
compare_figure("FigS4_retrospective.png")
