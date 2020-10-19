
require(reshape2)
require(ggpubr)
require(ggplot2)
require(ggsci)

outFolder <- "CMP_SAMEFAMINSAMETAD_REGIONOCC_WITHRANDOM"
dir.create(outFolder)

plotType <- "svg"

myHeight <- 6
myWidth <- 7

cat_labels <- c(
  cat1 = "1",
  cat2 = "2-3",
  cat3 = "4-5",
  cat4 = ">5"
)
nRandom <- 100

########################## PREPARE OBSERVED

obs_data <- get(load("SAMEFAMINSAMETAD_REGIONOCC/all_cat_nbrUniqueTADs.Rdata"))
obs_plot_dt <- melt(obs_data)
obs_plot_dt$nbr_cat <- gsub(".+_(.+)", "\\1", obs_plot_dt$variable)
obs_plot_dt$data <- gsub("(.+)_.+", "\\1", obs_plot_dt$variable)
obs_plot_dt$data <- ifelse(obs_plot_dt$data == "obs", "observed", 
                       ifelse(obs_plot_dt$data == "rd", "random", NA))
obs_plot_dt$nbr_cat_label <- cat_labels[obs_plot_dt$nbr_cat]
stopifnot(!is.na(obs_plot_dt))
tot_dt <- aggregate(value ~ data , FUN=sum, data=obs_plot_dt)
tot_by_data <- setNames(tot_dt$value, tot_dt$data)
obs_agg_dt <- aggregate(value ~ data + nbr_cat_label, FUN=sum, data=obs_plot_dt)
obs_agg_dt$value_ratio <- obs_agg_dt$value/tot_by_data[obs_agg_dt$data]
obs_agg_dt$value_log10 <- log10(obs_agg_dt$value)
obs_agg_dt$nbr_cat_label <- factor(obs_agg_dt$nbr_cat_label, levels=cat_labels)
stopifnot(!is.na(obs_agg_dt$nbr_cat_label))
obs_agg_dt$nbr_cat_label <- factor(obs_agg_dt$nbr_cat_label, levels=cat_labels)
stopifnot(!is.na(obs_agg_dt$nbr_cat_label))

# ggbarplot(obs_agg_dt, x="data", y="value_ratio", fill="nbr_cat_label",
#           xlab = "", ylab = "ratio unique TADs")+
#   scale_fill_nejm()+
#   labs(fill="") + 
#   ggtitle(plotTit, subtitle=subTit)+
#   theme(
#     plot.title = element_text(size=16, face = "bold", hjust=0.5),
#     plot.subtitle = element_text(size=14, face = "italic", hjust=0.5)
#   )


########################## PREPARE RANDOMMIDPOSDISC

randomdisc_data <- get(load("SAMEFAMINSAMETAD_REGIONOCC_RANDOMMIDPOSDISC//all_cat_nbrUniqueTADs.Rdata"))
randomdisc_plot_dt <- melt(randomdisc_data)
randomdisc_plot_dt$nbr_cat <- gsub(".+_(.+)", "\\1", randomdisc_plot_dt$variable)
randomdisc_plot_dt$data <- gsub("(.+)_.+", "\\1", randomdisc_plot_dt$variable)
randomdisc_plot_dt$data <- ifelse(randomdisc_plot_dt$data == "obs", "RANDOMMIDPOSDISC", 
                           ifelse(randomdisc_plot_dt$data == "rd", "random", NA))
randomdisc_plot_dt$nbr_cat_label <- cat_labels[randomdisc_plot_dt$nbr_cat]
stopifnot(!is.na(randomdisc_plot_dt))
tot_dt <- aggregate(value ~ data , FUN=sum, data=randomdisc_plot_dt)
tot_by_data <- setNames(tot_dt$value, tot_dt$data)
randomdisc_agg_dt <- aggregate(value ~ data + nbr_cat_label, FUN=sum, data=randomdisc_plot_dt)
randomdisc_agg_dt$value_ratio <- randomdisc_agg_dt$value/tot_by_data[randomdisc_agg_dt$data]
randomdisc_agg_dt$value_log10 <- log10(randomdisc_agg_dt$value)
randomdisc_agg_dt$nbr_cat_label <- factor(randomdisc_agg_dt$nbr_cat_label, levels=cat_labels)
stopifnot(!is.na(randomdisc_agg_dt$nbr_cat_label))
randomdisc_agg_dt$nbr_cat_label <- factor(randomdisc_agg_dt$nbr_cat_label, levels=cat_labels)
stopifnot(!is.na(randomdisc_agg_dt$nbr_cat_label))

# ggbarplot(randomdisc_agg_dt, x="data", y="value_ratio", fill="nbr_cat_label",
#           xlab = "", ylab = "ratio unique TADs")+
#   scale_fill_nejm()+
#   labs(fill="") + 
#   ggtitle(plotTit, subtitle=subTit)+
#   theme(
#     plot.title = element_text(size=16, face = "bold", hjust=0.5),
#     plot.subtitle = element_text(size=14, face = "italic", hjust=0.5)
#   )


########################## PREPARE RANDOMMIDPOSSTRICT

randomstrict_data <- get(load("SAMEFAMINSAMETAD_REGIONOCC_RANDOM//all_cat_nbrUniqueTADs.Rdata"))
randomstrict_plot_dt <- melt(randomstrict_data)
randomstrict_plot_dt$nbr_cat <- gsub(".+_(.+)", "\\1", randomstrict_plot_dt$variable)
randomstrict_plot_dt$data <- gsub("(.+)_.+", "\\1", randomstrict_plot_dt$variable)
randomstrict_plot_dt$data <- ifelse(randomstrict_plot_dt$data == "obs", "RANDOMMIDPOSSTRICT", 
                                  ifelse(randomstrict_plot_dt$data == "rd", "random", NA))
randomstrict_plot_dt$nbr_cat_label <- cat_labels[randomstrict_plot_dt$nbr_cat]
stopifnot(!is.na(randomstrict_plot_dt))
tot_dt <- aggregate(value ~ data , FUN=sum, data=randomstrict_plot_dt)
tot_by_data <- setNames(tot_dt$value, tot_dt$data)
randomstrict_agg_dt <- aggregate(value ~ data + nbr_cat_label, FUN=sum, data=randomstrict_plot_dt)
randomstrict_agg_dt$value_ratio <- randomstrict_agg_dt$value/tot_by_data[randomstrict_agg_dt$data]
randomstrict_agg_dt$value_log10 <- log10(randomstrict_agg_dt$value)
randomstrict_agg_dt$nbr_cat_label <- factor(randomstrict_agg_dt$nbr_cat_label, levels=cat_labels)
stopifnot(!is.na(randomstrict_agg_dt$nbr_cat_label))
randomstrict_agg_dt$nbr_cat_label <- factor(randomstrict_agg_dt$nbr_cat_label, levels=cat_labels)
stopifnot(!is.na(randomstrict_agg_dt$nbr_cat_label))

# ggbarplot(randomstrict_agg_dt, x="data", y="value_ratio", fill="nbr_cat_label",
#           xlab = "", ylab = "ratio unique TADs")+
#   scale_fill_nejm()+
#   labs(fill="") + 
#   ggtitle(plotTit, subtitle=subTit)+
#   theme(
#     plot.title = element_text(size=16, face = "bold", hjust=0.5),
#     plot.subtitle = element_text(size=14, face = "italic", hjust=0.5)
#   )


stopifnot(nrow(randomstrict_data) == nrow(randomdisc_data))
stopifnot(nrow(randomstrict_data) == nrow(obs_data))

all_agg_dt <- rbind(randomstrict_agg_dt, rbind(randomdisc_agg_dt, obs_agg_dt))
all_agg_dt <- all_agg_dt[all_agg_dt$data != "random",]
all_agg_dt$data <- factor(all_agg_dt$data, levels=c("observed", "RANDOMMIDPOSDISC", "RANDOMMIDPOSSTRICT"))
stopifnot(!is.na(all_agg_dt$data))


# subTit <- paste0("# DS = ", length(all_ds_results), "; # permut = ", nRandom)
subTit <- ""
plotTit <- "ratio unique TADs by family (by occurence)"

out_p <- ggbarplot(all_agg_dt, x="data", y="value_ratio", fill="nbr_cat_label",
          xlab = "", ylab = "ratio unique TADs")+
  scale_fill_nejm()+
  labs(fill="") + 
  ggtitle(plotTit, subtitle=subTit)+
  theme(
    plot.title = element_text(size=16, face = "bold", hjust=0.5),
    plot.subtitle = element_text(size=14, face = "italic", hjust=0.5)
  ) +
  # geom_hline(yintercept=all_agg_dt$value_ratio[all_agg_dt$data=="observed" & 
  #                                                (as.character(all_agg_dt$nbr_cat_label) == "1" |as.character(all_agg_dt$nbr_cat_label) == ">5") ])
  geom_hline(yintercept=1- all_agg_dt$value_ratio[all_agg_dt$data=="observed" & as.character(all_agg_dt$nbr_cat_label) == "1"  ], color="red") + 
geom_hline(yintercept=all_agg_dt$value_ratio[all_agg_dt$data=="observed" & as.character(all_agg_dt$nbr_cat_label) == ">5" ], color="red", linetype=1)


outFile <- file.path(outFolder, paste0("nbrUniqueTADsByFam_ratio_by_cat_withRandom_barplot.", plotType))
ggsave(out_p, filename=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile,  "\n"))
