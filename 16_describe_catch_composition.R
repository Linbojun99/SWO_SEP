# =====================================================================================================
# Step 16. Catch composition of the tuna-targeting longline fleets (Section 2.3, Section S7, Table S5, Figs. S5-S8)
# 第 16 步  金枪鱼延绳钓船队的渔获组成（第 2.3 节，补充材料 S7，表 S5，图 S5-S8）
#
# For each tuna-targeting longline fleet, uses the eligible records with at least one retained albacore, bigeye
# tuna, yellowfin tuna, blue marlin or striped marlin. The square-root proportions of the five species are
# centred on the fleet means and summarized by principal components; the fewest components explaining at
# least 80% of the variance (at most three) are retained, and two, three or four k-means clusters of those
# scores are compared by mean silhouette width (Table S5). Writes the species loadings, the mean composition
# of each cluster and the annual share of records in each cluster, and draws Figs. S5 and S6. The sensitivity
# analysis repeats all of this with swordfish included in the composition (Figs. S7 and S8). No index is fitted.
# 对每支金枪鱼延绳钓船队，使用至少留存一尾长鳍金枪鱼、大眼金枪鱼、黄鳍金枪鱼、蓝枪鱼或条纹四鳍旗鱼的合格记录。五个
# 鱼种比例的平方根以船队均值中心化后进行主成分分析；保留累计解释至少 80% 方差的最少主成分（最多三个），并按平均
# 轮廓宽度比较这些得分的两类、三类或四类 k 均值聚类（表 S5）。写出物种载荷、各聚类的平均组成和各聚类记录的年度
# 占比，并绘制图 S5 和 S6。敏感性分析在组成中加入剑鱼后重复上述全部步骤（图 S7 和 S8）。本步不拟合指数。
#
# Reads / 读取:  analysis_walkthrough/output/input.rds (step 04), config/directed_fishery.json
# Writes / 写出: analysis_walkthrough/output/composition_selection.csv (Table S5 / 表 S5),
#                china_cluster_swordfish.csv, composition_loadings.csv, cluster_profiles.csv,
#                cluster_annual_shares.csv and the same three with the suffix _inclusive,
#                FigS5_fleet_pca, FigS6_fleet_clusters, FigS7_fleet_pca_inclusive, FigS8_fleet_clusters_inclusive (.pdf/.png)
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/16_describe_catch_composition.R
# =====================================================================================================

suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
source("analysis_walkthrough/model_functions.R")
out <- "analysis_walkthrough/output"
cat("Step 16: catch composition (Table S5, Figs. S5-S8) / 第 16 步：渔获组成（表 S5、图 S5-S8）\n")
cfg <- jsonlite::fromJSON("config/directed_fishery.json")
d <- readRDS(file.path(out, "input.rds"))$d
flags <- c(JPN="Japan", KOR="Korea", CHN="China", TWN="Chinese Taipei", VUT="Vanuatu")
five <- c(ALBn="Albacore", BETn="Bigeye\ntuna", YFTn="Yellowfin\ntuna", BUMn="Blue\nmarlin", MLSn="Striped\nmarlin")
stopifnot(identical(cfg$companion, names(five)))
selection <- list(); swordfish <- list()

for (inclusive in c(FALSE, TRUE)) {
  variant <- if (inclusive) "include" else "exclude"
  suffix <- if (inclusive) "_inclusive" else ""
  species <- if (inclusive) c(five, SWOn="Swordfish") else five
  cat(sprintf("Swordfish %s in the composition / 组成中%s剑鱼\n", if (inclusive) "included" else "excluded",
              if (inclusive) "包含" else "排除"))

  # ---- Scores, loadings and clusters of each fleet / 各船队的得分、载荷与聚类 ----
  # The components and clusters are those of the composition models of step 17 (seed cfg$seed + fleet number).
  # 主成分和聚类与第 17 步组成模型所用的相同（种子 cfg$seed + 船队编号）。
  scores <- loadings <- profiles <- annual <- list()
  for (flag in names(flags)) {
    fleet <- d[d$Flag==flag & d$common, ]
    counts <- fleet[names(species)]
    fm <- targeting_learn(counts, seed=cfg$seed+match(flag, cfg$flags))
    pc <- seq_len(fm$n_pc)
    feat <- targeting_apply(fm, counts)
    stopifnot(isTRUE(all.equal(as.matrix(feat[paste0("PC",pc)]), fm$pca$x[,pc], check.attributes=FALSE)))
    selection[[paste(variant, flag)]] <- data.frame(variant, flag, records=nrow(fleet), retained_pcs=fm$n_pc,
      retained_variance=fm$retained_variance, fm$selection, selected=fm$selection$k==fm$k)
    scores[[flag]] <- data.frame(flag, feat[c("PC1","PC2")])
    w <- as.data.frame(as.table(fm$pca$rotation[,pc,drop=FALSE]), responseName="loading")
    names(w)[1:2] <- c("species","component")
    w$variance_percent <- 100*fm$explained_variance[match(w$component,paste0("PC",pc))]
    loadings[[flag]] <- data.frame(flag,w)
    # Mean number proportions in each cluster, each record with equal weight, and the annual share of records.
    # 各聚类的平均尾数比例（每条记录权重相等），以及各聚类记录的年度占比。
    props <- as.matrix(counts)/rowSums(counts)
    sizes <- table(feat$cluster)
    means <- rowsum(props,feat$cluster)/as.numeric(sizes)
    profile <- as.data.frame(as.table(means),responseName="mean_share")
    names(profile)[1:2] <- c("cluster","species")
    profiles[[flag]] <- data.frame(flag, profile, records=as.integer(sizes[as.character(profile$cluster)]))
    yearly <- as.data.frame(table(year=fleet$Year,cluster=feat$cluster),responseName="records")
    yearly$year <- as.integer(as.character(yearly$year))
    yearly$share <- yearly$records/ave(yearly$records,yearly$year,FUN=sum)
    stopifnot(all(abs(rowSums(means)-1)<1e-12),sum(yearly$records)==nrow(fleet))
    annual[[flag]] <- data.frame(flag,yearly)
    # For China, the mean share of swordfish in the six-species catch by number in each cluster.
    # 对中国，计算各聚类中剑鱼占六种鱼尾数的平均比例。
    if (flag=="CHN") {
      share <- fleet$SWOn/rowSums(fleet[c(cfg$companion,"SWOn")])
      swordfish[[variant]] <- data.frame(variant, cluster=rownames(means), records=as.integer(sizes),
        leading_companion=cfg$companion[max.col(means[,cfg$companion],ties.method="first")],
        mean_percent_six_species_by_number=100*as.numeric(tapply(share,feat$cluster,mean)))
    }
  }
  stack_fleets <- function(x) {
    x <- do.call(rbind,x)
    x$flag <- factor(x$flag,levels=names(flags))
    rownames(x) <- NULL
    x
  }
  scores <- stack_fleets(scores); loadings <- stack_fleets(loadings)
  profiles <- stack_fleets(profiles); annual <- stack_fleets(annual)
  write.csv(loadings,file.path(out,paste0("composition_loadings",suffix,".csv")),row.names=FALSE)
  write.csv(profiles,file.path(out,paste0("cluster_profiles",suffix,".csv")),row.names=FALSE)
  write.csv(annual,file.path(out,paste0("cluster_annual_shares",suffix,".csv")),row.names=FALSE)

  # ---- Figs. S5 and S7: scores and loadings / 图 S5 和 S7：得分与载荷 ----
  fleet_labels <- flags
  fleet_labels["TWN"] <- "Chinese\nTaipei"
  fleet_rows <- function(free_y=FALSE) {
    facet_grid(rows=vars(flag), scales=if (free_y) "free_y" else "fixed",
               switch="y", labeller=labeller(flag=fleet_labels))
  }
  style <- theme_bw(base_size=9,base_family="Arial") +
    theme(panel.grid=element_blank(),panel.border=element_rect(colour="grey30",linewidth=.48),
          strip.background=element_blank(),strip.placement="outside",
          strip.text.y.left=element_text(angle=0,face="bold",hjust=0),
          axis.text=element_text(colour="grey15"),axis.ticks=element_line(linewidth=.4),
          panel.spacing.y=grid::unit(3.6,"mm"),plot.title=element_text(face="bold",size=9),
          legend.position="bottom",legend.title=element_text(size=8),
          legend.key.height=grid::unit(3,"mm"))
  hidden_strip <- theme(strip.text.y.left=element_blank())
  species_x <- scale_x_discrete(limits=names(species),labels=species,expand=c(0,0))
  cell_y <- scale_y_discrete(expand=c(0,0))
  loadings$axis <- with(loadings, sprintf("%s (%.1f%%)",component,variance_percent))
  loadings$axis <- factor(loadings$axis,levels=rev(unique(loadings$axis[order(loadings$component)])))
  score_plot <- ggplot(scores,aes(PC1,PC2)) +
    geom_bin_2d(binwidth=c(.05,.05)) + fleet_rows() +
    scale_fill_gradient(low="#bdceda",high="#153c5c",transform="log10",
                        breaks=c(1,10,100,1000),name="Records per bin") +
    scale_y_continuous(breaks=c(-.5,0,.5)) +
    labs(x="PC1 score",y="PC2 score",title="(a) Component scores") + style
  loadings$weight_label <- sub("-0.00", "0.00", sprintf("%.2f",loadings$loading), fixed=TRUE)
  loading_plot <- ggplot(loadings,aes(species,axis,fill=loading)) +
    geom_tile(colour="white",linewidth=.4) +
    geom_text(aes(label=weight_label,colour=abs(loading)>.58),size=3.05) +
    scale_colour_manual(values=c("FALSE"="black","TRUE"="white"),guide="none") +
    scale_fill_gradient2(low="#984328",mid="white",high="#214e72",limits=c(-1,1),
                         midpoint=0,breaks=c(-1,0,1),name="Loading") +
    species_x + cell_y + fleet_rows(TRUE) +
    labs(x=NULL,y=NULL,title="(b) Species loadings") + style + hidden_strip
  pca <- (score_plot | loading_plot) + plot_layout(widths=c(1,1.45),guides="keep") &
    theme(legend.position="bottom",legend.title.position="top",
          legend.title=element_text(size=8,hjust=.5),legend.key.width=grid::unit(7,"mm"))

  # ---- Figs. S6 and S8: cluster composition and annual shares / 图 S6 和 S8：聚类组成与年度占比 ----
  profiles$label <- with(profiles,sprintf("C%s (%s)",cluster,format(records,big.mark=",",trim=TRUE)))
  profiles$label <- factor(profiles$label,levels=rev(unique(profiles$label[order(profiles$cluster)])))
  profiles$percent <- 100*profiles$mean_share
  profiles$percent_label <- ifelse(profiles$percent>0 & profiles$percent<.5,"<1",sprintf("%.0f",profiles$percent))
  if (inclusive) {
    swo <- profiles$species=="SWOn"
    profiles$percent_label[swo] <- sprintf("%.1f",profiles$percent[swo])
  }
  profile_plot <- ggplot(profiles,aes(species,label,fill=percent)) +
    geom_tile(colour="white",linewidth=.4) +
    geom_text(aes(label=percent_label,colour=percent>55),size=3.05) +
    scale_colour_manual(values=c("FALSE"="black","TRUE"="white"),guide="none") +
    scale_fill_gradient(low="white",high="#214e72",limits=c(0,100),guide="none") +
    species_x + cell_y + fleet_rows(TRUE) +
    labs(x=NULL,y=NULL,title="(a) Mean catch composition (%)") + style
  annual_plot <- ggplot(annual,aes(year,100*share,fill=cluster)) +
    geom_col(width=.82,position=position_stack(reverse=TRUE)) + fleet_rows() +
    scale_x_continuous(limits=cfg$years+c(-.5,.5),breaks=unique(c(seq(cfg$years[1],cfg$years[2],by=5),cfg$years[2])),expand=c(0,0)) +
    scale_y_continuous(breaks=c(0,50,100),expand=c(0,0)) +
    scale_fill_manual(values=c("1"="#294f6e","2"="#d39655","3"="#6ba79a"),
                      labels=function(x) paste0("C",x),name="Fleet-specific cluster") +
    labs(x="Year",y="Monthly records (%)",title="(b) Records by cluster and year") +
    style + hidden_strip
  clusters <- (profile_plot | annual_plot) + plot_layout(widths=c(if (inclusive) 1.4 else 1.2,1),guides="collect") &
    theme(legend.position="bottom")
  save_pair <- function(stem,p) {
    ggsave(file.path(out,paste0(stem,".pdf")),p,width=190,height=193,units="mm",device=grDevices::cairo_pdf,bg="white")
    ggsave(file.path(out,paste0(stem,".png")),p,width=190,height=193,units="mm",dpi=600,bg="white")
  }
  save_pair(if (inclusive) "FigS7_fleet_pca_inclusive" else "FigS5_fleet_pca",pca)
  save_pair(if (inclusive) "FigS8_fleet_clusters_inclusive" else "FigS6_fleet_clusters",clusters)
}
selection <- do.call(rbind, selection); rownames(selection) <- NULL
swordfish <- do.call(rbind, swordfish)
write.csv(selection, file.path(out, "composition_selection.csv"), row.names=FALSE)
write.csv(swordfish, file.path(out, "china_cluster_swordfish.csv"), row.names=FALSE)
print(subset(selection, variant=="exclude")[c("flag","retained_pcs","retained_variance","k","silhouette","selected")],
      digits=3, row.names=FALSE)
print(swordfish, digits=3, row.names=FALSE)

# ---- Compare with the archived results / 与存档结果比较 ----
# The archived candidate solutions and Chinese swordfish shares are kept in the two composition receipts.
# 存档的候选方案和中国剑鱼比例保存在两份组成记录文件中。
source("analysis_walkthrough/common.R")
receipt <- function(variant) jsonlite::fromJSON(file.path(archive, paste0("composition_figure_receipt",
  if (variant=="include") "_inclusive" else "", ".json")))$checks
old_selection <- do.call(rbind, lapply(c("exclude","include"), function(v) do.call(rbind, lapply(names(flags), function(f) {
  x <- receipt(v)[[f]]
  data.frame(variant=v, flag=f, records=x$records, retained_pcs=x$retained_pcs, retained_variance=x$retained_variance,
             x$candidates, selected=x$candidates$k==as.integer(names(x$cluster_counts))[length(x$cluster_counts)])
}))))
compare_table("composition_selection.csv", old=old_selection)
compare_table("china_cluster_swordfish.csv", old=do.call(rbind, lapply(c("exclude","include"), function(v)
  data.frame(variant=v, receipt(v)$CHN$swordfish_by_cluster))))
for (name in c("composition_loadings","cluster_profiles","cluster_annual_shares"))
  for (suffix in c("","_inclusive")) compare_table(paste0(name, suffix, ".csv"))
for (name in c("FigS5_fleet_pca","FigS6_fleet_clusters","FigS7_fleet_pca_inclusive","FigS8_fleet_clusters_inclusive"))
  compare_figure(paste0(name, ".png"))
