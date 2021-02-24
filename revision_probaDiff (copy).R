
tadSignifThresh <- 0.01

# 1st part => can be done for all data
# pval vs. inter/intra 
final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))
final_table_DT$regionID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
final_table_DT$signif_lab <- ifelse(final_table_DT$adjPvalComb <= 0.01, "signif.", "not signif.")

all_inter_intra1_dt <- get(load("REVISION_INTER_INTRA_PROBA/all_inter_intra_dt.Rdata"))
all_inter_intra2_dt <- get(load("REVISION_INTER_INTRA_PROBA2/all_inter_intra_dt.Rdata"))
stopifnot(! all_inter_intra1_dt$hicds %in% all_inter_intra2_dt$hicds)
all_inter_intra_dt <- rbind(all_inter_intra1_dt, all_inter_intra2_dt)
stopifnot(final_table_DT$hicds %in% all_inter_intra_dt$hicds)

# compute inter/intra for the next
all_inter_intra_dt$interOverIntraNorm_next <- all_inter_intra_dt$mean_inter_nextNorm/all_inter_intra_dt$mean_intraNorm
# compute inter/intra for the prev
all_inter_intra_dt$interOverIntraNorm_prev <- all_inter_intra_dt$mean_inter_prevNorm/all_inter_intra_dt$mean_intraNorm
# compute the average next and prev inter/intra
all_inter_intra_dt$interOverIntraNorm_mean <- sapply(1:nrow(all_inter_intra_dt), function(x) {
  # if there was no next TAD -> keep only prev ratio
  if(is.na(all_inter_intra_dt$mean_inter_nextNorm[x])) {
    return(all_inter_intra_dt$interOverIntraNorm_prev[x])
  }
  # if there was no prev TAD -> keep only next ratio
  if(is.na(all_inter_intra_dt$mean_inter_prevNorm[x])) {
    return(all_inter_intra_dt$interOverIntraNorm_next[x])
  }
  return(0.5*(all_inter_intra_dt$interOverIntraNorm_next[x]+all_inter_intra_dt$interOverIntraNorm_prev[x]))
})

###### for the not norm
# compute inter/intra for the next
all_inter_intra_dt$interOverIntra_next <- all_inter_intra_dt$mean_inter_next/all_inter_intra_dt$mean_intra
# compute inter/intra for the prev
all_inter_intra_dt$interOverIntra_prev <- all_inter_intra_dt$mean_inter_prev/all_inter_intra_dt$mean_intra
# compute the average next and prev inter/intra
all_inter_intra_dt$interOverIntra_mean <- sapply(1:nrow(all_inter_intra_dt), function(x) {
  # if there was no next TAD -> keep only prev ratio
  if(is.na(all_inter_intra_dt$mean_inter_next[x])) {
    return(all_inter_intra_dt$interOverIntra_prev[x])
  }
  # if there was no prev TAD -> keep only next ratio
  if(is.na(all_inter_intra_dt$mean_inter_prev[x])) {
    return(all_inter_intra_dt$interOverIntra_next[x])
  }
  return(0.5*(all_inter_intra_dt$interOverIntra_next[x]+all_inter_intra_dt$interOverIntra_prev[x]))
})




all_inter_intra_dt$hicds_regionID <- file.path(all_inter_intra_dt$hicds, all_inter_intra_dt$region)
final_table_DT$hicds_regionID <- file.path(final_table_DT$hicds, final_table_DT$region)
stopifnot(final_table_DT$hicds_regionID %in% all_inter_intra_dt$hicds_regionID)

merge_dt <- merge(final_table_DT, all_inter_intra_dt, all=FALSE, by=c("hicds", "region", "hicds_regionID"))

qts <- quantile(na.omit(merge_dt$interOverIntraNorm_mean), probs=c(0.05, 0.95))


sub_dt <- merge_dt[merge_dt$interOverIntraNorm_mean >= qts[1] & merge_dt$interOverIntraNorm_mean <= qts[2] ,]

p3 <- ggdensity(sub_dt,
                x = paste0("interOverIntraNorm_mean"),
                y = "..density..",
                # combine = TRUE,                  # Combine the 3 plots
                xlab = paste0("interOverIntraNorm_mean"),
                # add = "median",                  # Add median line.
                rug = FALSE,                      # Add marginal rug
                color = "signif_lab",
                fill = "signif_lab",
                palette = "jco"
) +
  ggtitle(plotTit, subtitle = mySub)+
  scale_color_manual(values=my_cols)+
  scale_fill_manual(values=my_cols)  +
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
  guides(color=FALSE)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  mytheme

outFile <- file.path(outFolder, paste0("tad_rankDiff_signif_notsignif_density.", plotType))
ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

qts <- quantile(na.omit(merge_dt$interOverIntra_mean), probs=c(0.05, 0.95))

sub_dt <- merge_dt[merge_dt$interOverIntra_mean >= qts[1] & merge_dt$interOverIntra_mean <= qts[2] ,]
ggdensity(sub_dt,
                x = paste0("interOverIntra_mean"),
                y = "..density..",
                # combine = TRUE,                  # Combine the 3 plots
                xlab = paste0("interOverIntra_mean"),
                # add = "median",                  # Add median line.
                rug = FALSE,                      # Add marginal rug
                color = "signif_lab",
                fill = "signif_lab",
                palette = "jco"
) 





# compute inter-intra for the next
all_inter_intra_dt$interMinusIntraNorm_next <- all_inter_intra_dt$mean_inter_nextNorm-all_inter_intra_dt$mean_intraNorm
# compute inter/intra for the prev
all_inter_intra_dt$interMinusIntraNorm_prev <- all_inter_intra_dt$mean_inter_prevNorm-all_inter_intra_dt$mean_intraNorm
# compute the average next and prev inter/intra
all_inter_intra_dt$interMinusIntraNorm_mean <- sapply(1:nrow(all_inter_intra_dt), function(x) {
  # if there was no next TAD -> keep only prev ratio
  if(is.na(all_inter_intra_dt$mean_inter_nextNorm[x])) {
    return(all_inter_intra_dt$interMinusIntraNorm_prev[x])
  }
  # if there was no prev TAD -> keep only next ratio
  if(is.na(all_inter_intra_dt$mean_inter_prevNorm[x])) {
    return(all_inter_intra_dt$interMinusIntraNorm_next[x])
  }
  return(0.5*(all_inter_intra_dt$interMinusIntraNorm_next[x]+all_inter_intra_dt$interMinusIntraNorm_prev[x]))
})


all_inter_intra_dt$hicds_regionID <- file.path(all_inter_intra_dt$hicds, all_inter_intra_dt$region)
final_table_DT$hicds_regionID <- file.path(final_table_DT$hicds, final_table_DT$region)
stopifnot(final_table_DT$hicds_regionID %in% all_inter_intra_dt$hicds_regionID)

merge_dt <- merge(final_table_DT, all_inter_intra_dt, all=FALSE, by=c("hicds", "region", "hicds_regionID"))

qts <- quantile(na.omit(merge_dt$interMinusIntraNorm_mean), probs=c(0.05, 0.95))


sub_dt <- merge_dt[merge_dt$interMinusIntraNorm_mean >= qts[1] & merge_dt$interMinusIntraNorm_mean <= qts[2] ,]

ggdensity(sub_dt,
                x = paste0("interMinusIntraNorm_mean"),
                y = "..density..",
                # combine = TRUE,                  # Combine the 3 plots
                xlab = paste0("interMinusIntraNorm_mean"),
                # add = "median",                  # Add median line.
                rug = FALSE,                      # Add marginal rug
                color = "signif_lab",
                fill = "signif_lab",
                palette = "jco"
) 
  





# 2nd part  => for matching pairs
refID signif. or not signif.

similar to normFC and tumorFC -> normInterIntraProba tumor

diff. inter/intra vs. signif/not signif. norm tumor
diff. inter/intra vs. norm tumor FC