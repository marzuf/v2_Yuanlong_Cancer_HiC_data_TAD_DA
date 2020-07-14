
# Rscript conserv_signif_obsShuffle_purityFilter.R

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


args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  purity_ds <- args[1]  
  purity_plot_name <- "EPIC"
} else{
  purity_ds <- ""
  purity_plot_name <- "aran"
}

### HARD-CODED - MAIN SETTINGS
corMet <- "pearson"
transfExpr <- "log10"

plotCex <- 1.4

fontFamily <- "Hershey"

mycols <- c("obs."="midnightblue","shuffle"="darkolivegreen")


outFolder <- file.path("CONSERV_SIGNIF_OBSSHUFFLE_PURITYFILTER", purity_ds, transfExpr)
dir.create(outFolder, recursive = TRUE)

mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[!grepl("PERMUT", all_hicds) & !grepl("RANDOM", all_hicds)]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

inFolder_perm <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_SHUFFLE_PURITYFILTER/", purity_plot_name, transfExpr)
inFolder_obs <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER", purity_ds, transfExpr)
                                                      #_PURITYFILTER/EPIC/log10/conserved_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata"))
# obsCons_dt <- get(load(file.path(inFolder_obs, 
#                                  "plot_matching_dt_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))


obs_data <- get(load(file.path(inFolder_obs,
                                  "conserved_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))

all_nCons <- lengths(obs_data)
dsByReg_dt <- do.call(rbind, lapply(obs_data, function(x) as.numeric(all_datasets %in% unique(dirname(x)))))
colnames(dsByReg_dt) <- all_datasets
stopifnot(sum(dsByReg_dt) == sum(all_nCons))
stopifnot(max(rowSums(dsByReg_dt)) == max(all_nCons))
obsCons_dt <- dsByReg_dt

m_obsCons_dt <- melt(obsCons_dt)
colnames(m_obsCons_dt) <- c("region", "dataset", "conserved")
stopifnot(grepl("region", m_obsCons_dt$region))
count_obsCons_dt <- aggregate(conserved ~ region, data=m_obsCons_dt, FUN=sum)
count_obsCons_dt <- count_obsCons_dt[order(count_obsCons_dt$conserved),]
count_obsCons_dt$region_rank <- 1:nrow(count_obsCons_dt)

m_obsCons_dt2 <- m_obsCons_dt[m_obsCons_dt$conserved == 1,]
colnames(m_obsCons_dt2) <- c("region", "dataset", "conserved")
stopifnot(grepl("region", m_obsCons_dt2$region))
count_obsCons_dt2 <- aggregate(dataset ~ region, data=m_obsCons_dt2, FUN=function(x)length(unique(x)))
colnames(count_obsCons_dt2)[colnames(count_obsCons_dt2) == "dataset"] <- "conserved"
count_obsCons_dt2 <- count_obsCons_dt2[order(count_obsCons_dt2$conserved),]
count_obsCons_dt2$region_rank <- 1:nrow(count_obsCons_dt2)
stopifnot(all.equal(count_obsCons_dt, count_obsCons_dt2, check.attributes=FALSE))

all_perm_dt <- foreach(i_perm=1:nPermut, .combine='rbind') %dopar% {
  # cons_dt <- get(load(file.path(inFolder_perm,
  #                               i_perm,
  #                               "plot_matching_dt_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
  perm_data <- get(load(file.path(inFolder_perm,
                                  i_perm,
                                 "conserved_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
  
  #lapply(perm_data, function(x) duplicated(dirname(x)))
  # !! corrected -> cannot be lengths(perm_data) because might be that several TADs of 1 dataset match a conserved region
  all_nCons <- unlist(lapply(perm_data, function(x) length(unique(dirname(x)))))
  dsByReg_dt <- do.call(rbind, lapply(perm_data, function(x) as.numeric(all_datasets %in% unique(dirname(x)))))
  colnames(dsByReg_dt) <- all_datasets
  stopifnot(sum(dsByReg_dt) == sum(all_nCons))
  stopifnot(max(rowSums(dsByReg_dt)) == max(all_nCons))
  cons_dt <- dsByReg_dt

  m_cons_dt <- melt(cons_dt)
  colnames(m_cons_dt) <- c("region", "dataset", "conserved")
  stopifnot(grepl("region", m_cons_dt$region))
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

obs_cons_dt <- count_obsCons_dt
saveFile <- file.path(outFolder, paste0("supp_fig4A_obs_cons_dt.Rdata"))
save(obs_cons_dt, file=saveFile, version=2)
cat(paste0("... written:" , saveFile, "\n"))

permut_cons_dt <- all_perm_dt
saveFile <- file.path(outFolder, paste0("supp_fig4A_permut_cons_dt.Rdata"))
save(permut_cons_dt, file=saveFile, version=2)
cat(paste0("... written:" , saveFile, "\n"))


# stop("-ok")


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
