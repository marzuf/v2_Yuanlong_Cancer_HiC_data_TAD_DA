
# Rscript conserv_signif_obsShuffle.R

require(ggplot2) 
require(ggsci)
require(reshape2)
require(foreach)
require(doMC)

registerDoMC(40)

atMin <- 2
nPermut <- 1000

plotType <- "svg"
myHeight <- 5
myHeightGG <- 5
myWidth <- 6
myWidthGG <- 7

plotCex <- 1.4

fontFamily <- "Hershey"

mycols <- c("obs."="midnightblue","shuffle"="darkolivegreen")

inFolder_obs <- "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2"
inFolder_perm <- "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_SHUFFLE"

outFolder <- "CONSERV_SIGNIF_OBSSHUFFLE"
dir.create(outFolder, recursive = TRUE)

obsCons_dt <- get(load(file.path(inFolder_obs, 
                                 "plot_matching_dt_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
m_obsCons_dt <- melt(t(obsCons_dt))
colnames(m_obsCons_dt) <- c("region", "dataset", "conserved")
count_obsCons_dt <- aggregate(conserved ~ region, data=m_obsCons_dt, FUN=sum)
count_obsCons_dt <- count_obsCons_dt[order(count_obsCons_dt$conserved),]
count_obsCons_dt$region_rank <- 1:nrow(count_obsCons_dt)

m_obsCons_dt2 <- m_obsCons_dt[m_obsCons_dt$conserved == 1,]
colnames(m_obsCons_dt2) <- c("region", "dataset", "conserved")
count_obsCons_dt2 <- aggregate(dataset ~ region, data=m_obsCons_dt2, FUN=function(x)length(unique(x)))
colnames(count_obsCons_dt2)[colnames(count_obsCons_dt2) == "dataset"] <- "conserved"
count_obsCons_dt2 <- count_obsCons_dt2[order(count_obsCons_dt2$conserved),]
count_obsCons_dt2$region_rank <- 1:nrow(count_obsCons_dt2)
stopifnot(all.equal(count_obsCons_dt, count_obsCons_dt2, check.attributes=FALSE))

all_perm_dt <- foreach(i_perm=1:nPermut, .combine='rbind') %dopar% {
  cons_dt <- get(load(file.path(inFolder_perm,
                                i_perm,
                                "plot_matching_dt_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
  m_cons_dt <- melt(t(cons_dt))
  colnames(m_cons_dt) <- c("region", "dataset", "conserved")
  count_cons_dt <- aggregate(conserved ~ region, data=m_cons_dt, FUN=sum)
  count_cons_dt <- count_cons_dt[order(count_cons_dt$conserved),]
  m_cons_dt2 <- m_cons_dt[m_cons_dt$conserved == 1,]
  colnames(m_cons_dt2) <- c("region", "dataset", "conserved")
  count_cons_dt2 <- aggregate(dataset ~ region, data=m_cons_dt2, FUN=function(x)length(unique(x)))
  colnames(count_cons_dt2)[colnames(count_cons_dt2) == "dataset"] <- "conserved"
  count_cons_dt2 <- count_cons_dt2[order(count_cons_dt2$conserved),]
  stopifnot(all.equal(count_cons_dt, count_cons_dt2, check.attributes=FALSE))
  stopifnot(count_cons_dt2$conserved >= atMin)
  stopifnot(count_cons_dt$conserved >= atMin)
  # m_cons_dt3 <- m_cons_dt2
  # m_cons_dt3$exprds <- dirname(as.character(m_cons_dt3$dataset))
  # count_cons_dt3_exprds <- aggregate(exprds ~ region, data=m_cons_dt3, FUN=function(x)length(unique(x)))
  # colnames(count_cons_dt3_exprds)[colnames(count_cons_dt3_exprds) == "exprds"] <- "conserved"
  # count_cons_dt3_exprds <- count_cons_dt3_exprds[order(count_cons_dt3_exprds$conserved),]
  # count_cons_dt3_exprds <- count_cons_dt3_exprds[count_cons_dt3_exprds$conserved >= atMin,]
  count_cons_dt$region_rank <- 1:nrow(count_cons_dt)
  count_cons_dt$region_rank_resc <- count_cons_dt$region_rank/max(count_cons_dt$region_rank)
  count_cons_dt$permut <- i_perm
  count_cons_dt
}
  

outFile <- file.path(outFolder, "all_perm_dt.Rdata")
save(all_perm_dt, file = outFile)


plotTit <- "Conserved signif. regions"
subTit <- paste0("# obs. cons. regions = ", nrow(count_obsCons_dt), "; # permut = ", nPermut)
  
outFile <- file.path(outFolder, paste0("obsCons_permCons_cons_line.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")
plot(NULL,
     xlim = c(0,1),
     ylim = range(c(count_obsCons_dt$conserved, all_perm_dt$conserved)),
     xlab = "Ranked regions (rescaled)",
     ylab = "Conserved in \"x\" datasets",
     cex.main=plotCex,
     cex.axis=plotCex,
     cex.lab=plotCex,
axes=F,
     main=plotTit)
axis(2, lwd=0, lwd.ticks=1)
axis(1, labels=F, tick=F)
box(type="L")
lines(x=count_obsCons_dt$region_rank/max(count_obsCons_dt$region_rank), 
      y  = count_obsCons_dt$conserved,
      col = mycols["obs."]
      )
mtext(side=3, text=subTit, font=3)

for(i in 1:nPermut) {
  tmp_perm_dt <- all_perm_dt[all_perm_dt$permut == i,]
  lines(x=tmp_perm_dt$region_rank_resc, 
        y  = tmp_perm_dt$conserved, 
        col=mycols["shuffle"])
}

legend("topleft",
       legend= c("obs.", paste0("shuffle (# permut. = ", nPermut, ")")),
       lty=1,
       col=c(mycols["obs."],mycols["shuffle"]),
       bty="n")


foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

stop("-ok")

nConsByPermut_dt <- aggregate(region ~ conserved + permut, FUN=function(x) length(unique(x)),data=all_perm_dt)
colnames(nConsByPermut_dt)[colnames(nConsByPermut_dt) == "region"] <- "nRegions"

nConsByPermut_dt$conserved <- factor(as.character(nConsByPermut_dt$conserved), 
                                     levels=as.character(2:max(c(count_obsCons_dt$conserved, all_perm_dt$conserved))))

agg_count_obsCons_dt <- aggregate(region ~ conserved, FUN=function(x)length(unique(x)), data=count_obsCons_dt)
colnames(agg_count_obsCons_dt)[colnames(agg_count_obsCons_dt) == "region"] <- "nRegions"

agg_count_obsCons_dt$conserved <- factor(as.character(agg_count_obsCons_dt$conserved), 
                                         levels=as.character(2:max(c(count_obsCons_dt$conserved, all_perm_dt$conserved))))

agg_count_obsCons_dt$type <- "obs."
agg_count_obsCons_dt$type <- factor(agg_count_obsCons_dt$type, levels=c("obs.", "shuffle"))

plotTit <- "Conserved signif. regions"
subTit <- paste0("# obs. cons. regions = ", nrow(count_obsCons_dt), "; # permut = ", nPermut)

stopifnot(nrow(count_obsCons_dt) == sum(agg_count_obsCons_dt$nRegions))

barbox_p <- ggplot(data=agg_count_obsCons_dt, 
         aes(x=conserved, y = nRegions, color=type,fill=type))+ 
  ggtitle(plotTit, subtitle = subTit)+
  geom_bar(stat="identity", alpha=0.2) + 
  labs(fill="", color="", x="Conserved in \"x\" datasets", y="# of regions")+
  scale_color_manual(values=mycols, breaks = names(mycols), labels=names(mycols), drop=F)+
  # scale_y_continuous(expand=c(0,5))+
  scale_fill_manual(values=mycols, drop=F)+
  geom_boxplot(data=nConsByPermut_dt, aes(x=conserved, y = nRegions), 
               color = mycols["shuffle"],
               inherit.aes = FALSE, notch = TRUE)+
  theme(
  text = element_text(family=fontFamily),
  panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
  panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
  panel.background = element_rect(fill = "transparent"),
  panel.grid.major.x =  element_blank(),
  panel.grid.minor.x =  element_blank(),
  axis.line = element_line(),
  axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
  axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
  axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
  axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
  plot.title = element_text(hjust=0.5, size = 16, face="bold"),
  plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
  legend.title = element_text(face="bold"),
  legend.text = element_text(size=12)
) 

outFile <- file.path(outFolder, paste0("obsCons_permCons_nCons_barplot_boxplot.", plotType))
ggsave(barbox_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))




all_perm_dt$consType <- ifelse(all_perm_dt$conserved == 2, "2",
                               ifelse(all_perm_dt$conserved > 2 & all_perm_dt$conserved <= 5, ">2,<=5",
                                      ifelse(all_perm_dt$conserved > 5 & all_perm_dt$conserved <= 10, ">5,<=10",
                                             ifelse(all_perm_dt$conserved > 10, ">10", NA))))
stopifnot(!is.na(all_perm_dt$consType))
cat_perm_dt <- aggregate(region~consType, data=all_perm_dt, function(x) length(x)/nPermut)
cat_perm_dt$ds_type <- "shuffle (mean)"

count_obsCons_dt$consType <- ifelse(count_obsCons_dt$conserved == 2, "2",
                               ifelse(count_obsCons_dt$conserved > 2 & count_obsCons_dt$conserved <= 5, ">2,<=5",
                                      ifelse(count_obsCons_dt$conserved > 5 & count_obsCons_dt$conserved <= 10, ">5,<=10",
                                             ifelse(count_obsCons_dt$conserved > 10, ">10", NA))))
stopifnot(!is.na(count_obsCons_dt$consType))
cat_obs_dt <- aggregate(region~consType, data=count_obsCons_dt, function(x) length(x))
cat_obs_dt$ds_type <- "obs."

plot_dt <- rbind(cat_obs_dt, cat_perm_dt)

plot_dt$consType <- factor(plot_dt$consType, levels=rev(c("2", ">2,<=5", ">5,<=10", ">10")))

ggsci_pal <- "uchicago"
ggsci_subpal <- ""


plotTit <- "Conserved signif. regions"
subTit <- paste0("# permut = ", nPermut)

p_cat <- ggplot(data=plot_dt, 
                   aes(x=ds_type, y = region, color=consType,fill=consType))+ 
  ggtitle(plotTit, subtitle = subTit)+
  geom_bar(stat="identity", position="stack") + 
  labs(fill="Conserved in\n\"x\" datasets", color="Conserved in\n\"x\" datasets", y="# of regions")+
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")")))+
  # scale_y_continuous(expand=c(0,5))+
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
  theme(
    text = element_text(family=fontFamily),
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x =  element_blank(),
    axis.line = element_line(),
    # axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    axis.text.x = element_text(size=14, hjust=0.5, vjust=0.5),
    plot.title = element_text(hjust=0.5, size = 16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.title = element_text(face="bold"),
    legend.text = element_text(size=12)
  ) 

outFile <- file.path(outFolder, paste0("obsCons_permCons_byCat_barplot.", plotType))
ggsave(p_cat, filename = outFile, height=myHeightGG, width=myWidthGG*0.8)
cat(paste0("... written: ", outFile, "\n"))



########## ~TRASH

# min_max_permDT <- do.call(rbind, by(all_perm_dt, all_perm_dt$region_rank, function(x){
#   data.frame(
#             region_rank = unique(x$region_rank),
#              minCons = min(x$conserved),
#              maxCons = max(x$conserved), stringsAsFactors = FALSE)
#   }))

# plot(NULL,
#      xlim = range(c(min_max_permDT$region_rank, count_obsCons_dt$region_rank)),
#      ylim = range(c(count_obsCons_dt$conserved, all_perm_dt$conserved)),
#      cex.main=plotCex,
#      cex.main=plotCex,
#      cex.main=plotCex,
#      main="")
#   
# lines(x=count_obsCons_dt$region_rank, y  = count_obsCons_dt$conserved)
# 
# for(i in 1:nPermut) {
#   
#   tmp_perm_dt <- all_perm_dt[all_perm_dt$permut == i,]
#   lines(x=tmp_perm_dt$region_rank, y  = tmp_perm_dt$conserved, col="grey")
# }

# polygon(c(min_max_permDT$region_rank, rev(min_max_permDT$region_rank)),
#         c(min_max_permDT$minCons, rev(min_max_permDT$maxCons)),
#         col="grey"
#         )  


# geom_boxplot(data=nConsByPermut_dt, aes(x=conserved, y = nRegions), inherit.aes = FALSE, notch = TRUE, outlier.shape=NA)+
# geom_boxplot(notch = TRUE, outlier.shape=NA)+
# geom_point(position=position_jitterdodge(),  alpha=0.5) +
# geom_jitter(data=nConsByPermut_dt, aes(x=conserved, y = nRegions), inherit.aes = FALSE,alpha=0.5)+
