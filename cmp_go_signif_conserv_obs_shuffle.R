# collect all

# Rscript cmp_go_signif_conserv_obs_shuffle.R

require(foreach)
require(doMC)
registerDoMC(40)
require(ggplot2)
fontFamily <- "Hershey"

plotType <- "svg"
myHeightGG <- 9
myWidthGG <- 9

outFolder <- "CMP_GO_SIGNIF_CONSERV_OBS_SHUFFLE"
dir.create(outFolder, recursive = TRUE)

nPerm <- 1000

signif_col <- "p.adjust"
signif_thresh <- 0.05

my_box_theme <- theme(
  text = element_text(family=fontFamily),
  panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
  panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
  panel.background = element_rect(fill = "transparent"),
  panel.grid.major.x =  element_blank(),
  panel.grid.minor.x =  element_blank(),
  axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
  axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
  axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
  axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
  plot.title = element_text(hjust=0.5, size = 16, face="bold"),
  plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
  legend.title = element_text(face="bold")
) 


all_perm_signifGO <- foreach(i_perm=1:nPerm) %dopar% {
  perm_go_dt <- get(load(file.path("GO_SIGNIF_ACROSS_HICDS_v2_SHUFFLE", i_perm, "conserved_signif_enrich_resultDT.Rdata")))
  perm_go_dt$ID[perm_go_dt[,paste0(signif_col)] <= signif_thresh]
}

obs_go_dt <- get(load(file.path("GO_SIGNIF_ACROSS_HICDS_v2", "conserved_signif_enrich_resultDT.Rdata")))
obs_signifGO <- obs_go_dt$ID[obs_go_dt[,paste0(signif_col)] <= signif_thresh]

nPermGO <- lengths(all_perm_signifGO)

plot_dt <- rbind(
  data.frame(
  GO_type = c("obs."),
  nGO = length(obs_signifGO),
  stringsAsFactors = FALSE
  ),
  data.frame(
    GO_type = c("shuffle"),
    nGO = lengths(all_perm_signifGO),
    stringsAsFactors = FALSE
  )
)


obsGO_perm <- unlist(lapply(obs_signifGO, function(curr_go) sum(unlist(lapply(all_perm_signifGO, function(sublist) curr_go %in% sublist)))))
names(obsGO_perm) <- obs_signifGO

plot_dt2 <- obs_go_dt[,c("ID", signif_col)]
plot_dt2 <- plot_dt2[plot_dt2[,paste0(signif_col)] <= signif_thresh,]
plot_dt2$ID_permCount <- obsGO_perm[plot_dt2$ID]
stopifnot(!is.na(plot_dt2))
plot_dt2 <- plot_dt2[order(plot_dt2[,c(signif_col)]),]
plot_dt2$ID <- factor(plot_dt2$ID, levels = plot_dt2$ID)

# scale_x_discrete(labels=function(x){sub("\\s", "\n", x)})
save(all_perm_signifGO, file=file.path(outFolder, "all_perm_signifGO.Rdata"))
save(plot_dt, file=file.path(outFolder, "plot_dt.Rdata"))
save(plot_dt2, file=file.path(outFolder, "plot_dt2.Rdata"))
# 
# load("CMP_GO_SIGNIF_CONSERV_OBS_SHUFFLE/plot_dt2.Rdata")
# load("CMP_GO_SIGNIF_CONSERV_OBS_SHUFFLE/plot_dt.Rdata")

# load("CMP_GO_SIGNIF_CONSERV_OBS_SHUFFLE/all_perm_signifGO.Rdata")




plotTop <- 10

strwdth <- 25

plotTit <- paste0("GO occurences in permut. data")
subTit <- paste0("# permut = ", nPerm, "; GO ", signif_col , " <= ", signif_thresh)

occ_p <- ggplot(plot_dt2[1:min(plotTop, nrow(plot_dt2)),], aes(x=ID, y=ID_permCount)) +
  ggtitle(plotTit, subtitle=subTit)+
  geom_bar(stat="identity") + 
  scale_x_discrete(labels=function(x) {
    unlist(lapply(strwrap(gsub("_", " ", x), width = strwdth, simplify=FALSE), function(x) paste0(x, collapse="\n")))})+
  labs(x=paste0("Ranked signif. GOs (top ", plotTop, ")") , y="# occ. in permut. data") + 
  my_box_theme+
  theme(
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    axis.line=element_line()
  )

outFile <- file.path(outFolder, paste0("nOcc_obsSignifGO_top", plotTop, ".", plotType))
ggsave(occ_p, file=outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

occ_p <- ggplot(plot_dt2, aes(x=ID, y=ID_permCount)) +
  ggtitle(plotTit, subtitle=subTit)+
  geom_bar(stat="identity") + 
  scale_x_discrete(labels=function(x) {
    unlist(lapply(strwrap(gsub("_", " ", x), width = strwdth, simplify=FALSE), function(x) paste0(x, collapse="\n")))})+
  labs(x=paste0("Ranked signif. GOs (top ", plotTop, ")") , y="# occ. in permut. data") + 
  my_box_theme+
  theme(
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    axis.line=element_line()
  )

outFile <- file.path(outFolder, paste0("nOcc_obsSignifGO_allSignif.", plotType))
ggsave(occ_p, file=outFile, height=myHeightGG, width=myWidthGG*2)
cat(paste0("... written: ", outFile, "\n"))



plotTit <- paste0("# signif. GOs")
subTit <- paste0("# permut = ", nPerm, "; GO ", signif_col , " <= ", signif_thresh)
  
nsig_p <- ggplot(data=plot_dt, aes(x=GO_type, y=nGO))+
  ggtitle(plotTit, subtitle=subTit)+
  geom_boxplot()+
  labs(x="Data", y="# signif. GOs") + 
  my_box_theme +
  theme(
    axis.text.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.line=element_line()
  )


outFile <- file.path(outFolder, paste0("nSignif_obs_permut.", plotType))
ggsave(nsig_p, file=outFile, height=myHeightGG*0.8, width=myWidthGG*0.8)
cat(paste0("... written: ", outFile, "\n"))


all_permGOs <- table(unlist(all_perm_signifGO))

perm_dt <- data.frame(
  GO=names(all_permGOs),
  nSignif=as.numeric(all_permGOs),
  stringsAsFactors = FALSE
)

perm_dt <- perm_dt[order(perm_dt$nSignif, decreasing = TRUE),]

perm_dt$GO <- factor(perm_dt$GO, levels=perm_dt$GO)

plotTit <- paste0("Permut GO occurences")
subTit <- paste0("# permut = ", nPerm, "; GO ", signif_col , " <= ", signif_thresh)

perm_p <- ggplot(perm_dt[1:min(plotTop, nrow(perm_dt)),], aes(x=GO, y=nSignif)) +
  ggtitle(plotTit, subtitle=subTit)+
  geom_bar(stat="identity") + 
  scale_x_discrete(labels=function(x) {
    unlist(lapply(strwrap(gsub("_", " ", x), width = strwdth, simplify=FALSE), function(x) paste0(x, collapse="\n")))})+
  labs(x=paste0("Most frequent signif. GOs (top ", plotTop, ")") , y="# occ. in permut. data") + 
  my_box_theme+
  theme(
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    axis.line=element_line()
  )

outFile <- file.path(outFolder, paste0("nOcc_permutSignifGO_top", plotTop, ".", plotType))
ggsave(perm_p, file=outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))





