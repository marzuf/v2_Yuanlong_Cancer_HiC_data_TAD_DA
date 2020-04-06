# Rscript cmp_nSignif_sampleSize.R

setDir <- "/media/electron"
setDir <- ""

require(foreach)
require(doMC)
registerDoMC(40)
require(ggplot2)
require(ggpubr)
require(reshape2)

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 9

outFolder <- "CMP_NSIGNIF_SAMPLESIZE"
dir.create(outFolder, recursive = TRUE)

tadSignifThresh <- 0.01
geneSignifThresh <- 0.01

init_hicds <- "ENCSR312KHQ_SK-MEL-5"
exprds <- "TCGAskcm_lowInf_highInf"

script1_name <- "1_runGeneDE"
script11_name <- "11sameNbr_runEmpPvalCombined"

nsub=""
nsub=20
all_tadSignif_dt <- foreach(nsub = c("", seq(from=20, to=100, by=20)), .combine='rbind') %dopar% {
  
  if(nsub == "") {
    hicds <- paste0(init_hicds,"_", nsub, "40kb")
  } else{
    hicds <- paste0(init_hicds,"_RANDOMSUB", nsub, "_40kb")
  }
  source(file.path("PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R")))  
  s1 <- get(load(file.path(setDir, sample1_file)))
  s2 <- get(load(file.path(setDir, sample2_file)))
  stopifnot(length(s1) == length(s2))
  if(nsub != "") stopifnot(length(s1) == nsub)
  
  pval <- get(load(file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, script11_name, "emp_pval_combined.Rdata")))
  adj_pval <- p.adjust(pval, method="BH")
  
  data.frame(
    hicds = hicds,
    exprds = exprds,
    nSamp1 = length(s1),
    nSamp2 = length(s2),
    region = names(adj_pval),
    adjCombPval = as.numeric(adj_pval),
    stringsAsFactors = FALSE
  )
  
  
}
outFile <- file.path(outFolder, "all_tadSignif_dt.Rdata")
save(all_tadSignif_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))


all_tadSignif_dt$nSamp <- all_tadSignif_dt$nSamp1 +  all_tadSignif_dt$nSamp2
all_tadSignif_dt <- all_tadSignif_dt[order(all_tadSignif_dt$nSamp),]
all_tadSignif_dt$sampLab <- paste0(all_tadSignif_dt$nSamp1, "+", all_tadSignif_dt$nSamp2 )

labLevels <- unique(as.character(all_tadSignif_dt$sampLab))


tadSignif_agg_dt <- aggregate(adjCombPval~hicds+exprds+sampLab,data=all_tadSignif_dt,FUN=function(x) mean(x <= tadSignifThresh))
colnames(tadSignif_agg_dt)[colnames(tadSignif_agg_dt) == "adjCombPval"] <- "ratioSignifTADs"

outFile <- file.path(outFolder, "tadSignif_agg_dt.Rdata")
save(tadSignif_agg_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))


all_geneSignif_dt <- foreach(nsub = c("", seq(from=20, to=100, by=20)), .combine='rbind') %dopar% {
  
  if(nsub == "") {
    hicds <- paste0("ENCSR312KHQ_SK-MEL-5_", nsub, "40kb")
  } else{
    hicds <- paste0("ENCSR312KHQ_SK-MEL-5_RANDOMSUB", nsub, "_40kb")
  }
  source(file.path("PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R")))  
  s1 <- get(load(file.path(setDir, sample1_file)))
  s2 <- get(load(file.path(setDir, sample2_file)))
  stopifnot(length(s1) == length(s2))
  if(nsub != "") stopifnot(length(s1) == nsub)
  
  de_dt <- get(load(file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, script1_name, "DE_topTable.Rdata")))
  
  data.frame(
    hicds = hicds,
    exprds = exprds,
    nSamp1 = length(s1),
    nSamp2 = length(s2),
    gene = as.character(de_dt$genes) ,
    adjPval = de_dt$adj.P.Val,
    stringsAsFactors = FALSE
  )
  
  
}
outFile <- file.path(outFolder, "all_geneSignif_dt.Rdata")
save(all_geneSignif_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))


all_geneSignif_dt$nSamp <- all_geneSignif_dt$nSamp1 +  all_geneSignif_dt$nSamp2
all_geneSignif_dt$nSamp <- all_geneSignif_dt[order(all_geneSignif_dt$nSamp),]
all_geneSignif_dt$sampLab <- paste0(all_geneSignif_dt$nSamp1, "+", all_geneSignif_dt$nSamp2 )


geneSignif_agg_dt <- aggregate(adjPval~hicds+exprds+sampLab,data=all_geneSignif_dt,FUN=function(x) mean(x <= geneSignifThresh))
colnames(geneSignif_agg_dt)[colnames(geneSignif_agg_dt) == "adjPval"] <- "ratioSignifGenes"

outFile <- file.path(outFolder, "geneSignif_agg_dt.Rdata")
save(geneSignif_agg_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

load("CMP_NSIGNIF_SAMPLESIZE/geneSignif_agg_dt.Rdata")
load("CMP_NSIGNIF_SAMPLESIZE/tadSignif_agg_dt.Rdata")

allSignif_agg_dt <- merge(geneSignif_agg_dt, tadSignif_agg_dt, by=c("hicds", "exprds", "sampLab"), all=TRUE)
stopifnot(!is.na(allSignif_agg_dt))


plot_dt <- melt(allSignif_agg_dt, id.vars = c("hicds", "exprds", "sampLab"))


plot_dt$sampLab <- factor(plot_dt$sampLab, levels=labLevels)

subTit <- "(variable sample size)"

col1 <- "darkblue"
col2 <- "darkorange"

signif_p <- ggplot(
  data = plot_dt,
 aes_string(x="sampLab",y="value", fill = "variable"))+
  geom_bar(stat="identity", position="dodge")+
  scale_x_discrete(name= "# of samples")+
  scale_y_continuous(name="Ratio signif. features")+
  scale_fill_manual(values = c(col1, col2)) +
  labs(fill  = paste0("") )
       

outFile <- file.path(outFolder, paste0(init_hicds, "_", exprds, "_signifSubSamples_boxplot.", plotType))
ggsave(plot = signif_p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

       
# nSignif_barplot <-  ggplot(plot_dt, aes(x = nSamp, y = value, fill = variable)) + 
#   geom_bar(position="dodge", stat="identity")
# 
# ggtitle(paste0(hicds, " - ", exprds), subtitle = paste0(subTit))+
#   scale_x_discrete(name=my_xlab)+
#   scale_y_continuous(name=paste0(my_ylab),
#                      breaks = scales::pretty_breaks(n = 20))+
#   
#   scale_shape_manual(
#     values = c(15,8),
#     breaks = c("noMut", "withMut"),
#     labels = c("not mut.", "mut.")
#   )+
#   
#   scale_color_manual(values=c(col1, col2))+
#   scale_fill_manual(values=c(col1, col2))+
#   
#   labs(fill  = paste0("Cond."), color=paste0("Cond."), shape=paste0("KEAP1|NFE2L2")) +
#   theme( 
#     plot.title = element_text(hjust = 0.5, face = "bold", size=16),
#     plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
#     panel.grid = element_blank(),
#     panel.grid.major.y = element_line(colour = "grey"),
#     panel.grid.minor.y = element_line(colour = "grey"),
#     axis.line.x= element_line(size = .2, color = "black"),
#     axis.line.y = element_line(size = .2, color = "black"),
#     axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
#     axis.text.x =element_text(color="black", hjust=0.5,vjust = 0.5, size=12, face="bold"),
#     # axis.ticks.x = element_blank(),
#     axis.title.y = element_text(color="black", size=14),
#     axis.title.x = element_text(color="black", size=14),
#     panel.border = element_blank(),
#     panel.background = element_rect(fill = "transparent"),
#     legend.background =  element_rect(),
#     legend.text = element_text(size=12),
#     legend.key = element_blank(),
#     legend.title = element_text(face="bold", size=12)
#   )

# outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_plot, "_allSamples_exprValues_boxplot_vShape.", plotType))
# ggsave(plot = p_var_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
# cat(paste0("... written: ", outFile, "\n"))





