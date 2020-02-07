

options(scipen=100)

# Rscript toy_fcc_scores.R

script_name <- "toy_fcc_scores.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")


outFolder <- file.path("TOY_FCC_SCORES")
dir.create(outFolder, recursive = TRUE)

plotCex <- 1.4
plotType <- "svg"
myHeight <- 7
myWidth <- 9
# plotType <- "png"
# myHeight <- 400
# myWidth <- 500


get_fcc <- function(fc_vect) {
  (2* sum(fc_vect < 0)/length(fc_vect) -1) *  (2* sum(abs(fc_vect[fc_vect<0]))/sum(abs(fc_vect)) -1)
}

get_ratioDown <- function(fc_vect) {
  sum(fc_vect < 0)/length(fc_vect) 
}

get_ratioFC <- function(fc_vect) {
  sum(abs(fc_vect[fc_vect<0]))/sum(abs(fc_vect))
}


################################################
# CASE 1: very discordant
################################################

very_discordant1 <- c(-0.2,-0.1,-0.4,-0.3,-0.25, 5)

very_discordant2 <- c(-2,-1,-0.4,-0.3,0.2, 4, 5)

very_discordant3 <- c(-4, 3, -0.2, 0.1, -1.5, 2, 0.1, 0.1, -0.05)



################################################
# CASE 2: fifty-fifty
################################################

fifty_fifty <- c(-1, -1.5, -2, -2.5, -3, 
                1, 1.5, 2, 2.5, 3)


################################################
# CASE 3: down concordant
################################################

down_concordant <- c(-2,-2.5,-3,-3.5, -1.5,-4)


################################################
# CASE 4: up concordant
################################################

up_concordant <- abs(down_concordant)


################################################
# CASE 5: mostly down concordant
################################################

almost_down_concordant <- c(-2,-2.5,-3,-3.5, -1.5,-4, 0.2, 0.4, 0.3)



################################################
# CASE 6: mostly up concordant
################################################

almost_up_concordant <- -1*almost_down_concordant

################################################
# CASE 7: weird 0 cases
################################################

weird0_1 <- c(-0.1, -0.2, -0.3, -0.4, 2.5, 3, 3.5, 4)

weird0_2 <- c(-5, 1, 1.5, 1.5, 1)




plot_vect = "very_discordant"


for(plot_vect in c( "very_discordant1", "very_discordant2", "very_discordant3", "fifty_fifty", "down_concordant", "up_concordant", "almost_down_concordant", "almost_up_concordant", "weird0_1", "weird0_2")) {
  
  fc_vect <- sort(get(plot_vect))
  
  plotTit <-  paste0("FCC = ", sprintf("%.4f", get_fcc(fc_vect)))
  
  
  outFile <- file.path(outFolder, paste0(plot_vect, "_toy_fcc.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  
  plot(fc_vect, type="h", col = as.numeric(fc_vect < 0)+1, lwd=3, axes=F,
       ylim = c(-max(abs(fc_vect)), max(abs(fc_vect))),
       xlab = "",
       ylab = "log2FC",
       main=plotTit,
       cex.axis=plotCex,
       cex.lab=plotCex,
       cex.main=plotCex)
  axis(1, pos=0, labels=FALSE, lwd.ticks = -1, at=0:(length(fc_vect)+1))
  axis(2)
  
  
  legend(
    "topleft",
    legend = c(
      paste0("ratioDown = ", sprintf("%.4f", get_ratioDown(fc_vect))),
      paste0("ratioFC = ", sprintf("%.4f", get_ratioFC(fc_vect)))
    ),
    bty="n"
  )
  
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}








