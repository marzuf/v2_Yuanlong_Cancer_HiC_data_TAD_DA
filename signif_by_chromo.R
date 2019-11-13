require(ggpubr)

plotType <- "png"
myHeight <- 7
myWidth <- 10

outFolder <- file.path("SIGNIF_BY_CHROMO")
dir.create(outFolder, recursive = TRUE)

all_chromo <- paste0("chr", 1:22)

final_DT <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

final_DT$chromo <- gsub("(chr.+)_.+", "\\1", final_DT$region)
final_DT$chromo <- factor(final_DT$chromo, levels=all_chromo)
stopifnot(!is.na(final_DT$chromo))

final_DT$adjPvalComb_log10 <- -log10(final_DT$adjPvalComb)

p_pvalByChromo <- ggboxplot(final_DT, x = "chromo", y = "adjPvalComb_log10", xlab="", ylab="adj. comb. p-val [-log10]",
                            title = paste0("all datasets, all TADs")
                            ) +
                theme(
                  axis.title.x = element_blank(),
                  axis.text.x = element_text(angle = 90)
                )

outFile <- file.path(outFolder, paste0("adjPval_allTADs_allDS_by_chromo.", plotType))
ggsave(p_pvalByChromo, filename = outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))

p_pvalSignifByChromo <- ggboxplot(final_DT[final_DT$adjPvalComb <= 0.01,], x = "chromo", y = "adjPvalComb_log10",
                                  title = paste0("all datasets, only signif. TADs")
                              ) +
                                theme(
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_text(angle = 90)
                                )
outFile <- file.path(outFolder, paste0("adjPval_signifTADs_allDS_by_chromo.", plotType))
ggsave(p_pvalSignifByChromo, filename = outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))



countSignif_DT <- aggregate(adjPvalComb ~ hicds+exprds+chromo, data=final_DT, function(x) sum(x<=0.01))

p_countSignifChromo <- ggboxplot(x="chromo", y="adjPvalComb",
                                 ylab=paste0("# signif. TADs"),
                                 data=countSignif_DT) +
                        theme(
                          axis.title.x = element_blank(),
                          axis.text.x = element_text(angle = 90)
                        )
                        
outFile <- file.path(outFolder, paste0("nbr_signifTADs_allDS_by_chromo.", plotType))
ggsave(p_countSignifChromo, filename = outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))




