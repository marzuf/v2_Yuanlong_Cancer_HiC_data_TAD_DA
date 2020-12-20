
# nSingleton_lab0 <- get_fract_lab0(vect_values=tad_sing_dt$nSingleton, range_levels = nSingleton_levels)
# stopifnot(nSingleton_lab0 == tad_sing_dt$nSingleton_lab)


get_fract_lab2 <- function(vect_values, range_levels) {
  stopifnot(is.numeric(vect_values))
  stopifnot(is.numeric(range_levels))
  range_levels <- sort(range_levels)
  my_cmd_first <- paste0("ifelse(vect_values <= ", range_levels[1], ",\"<=", range_levels[1], "\",")
  my_cmd_last <- paste0("ifelse(vect_values > ", range_levels[length(range_levels)], ",\">",range_levels[length(range_levels)],"\",NA")
  my_cmd_end <- paste0(rep(")", length(range_levels)+1), collapse="")
  my_cmd_mid <- paste0("ifelse(vect_values >",range_levels[1:(length(range_levels)-1)], " & vect_values <=", range_levels[2:(length(range_levels))],
                       ",\">",range_levels[1:(length(range_levels)-1)], " & <=", range_levels[2:(length(range_levels))],"\","
  )
  full_cmd <- paste0(my_cmd_first,
                     paste0(my_cmd_mid, collapse=""),
                     my_cmd_last,
                     my_cmd_end, collapse=",")
  return(eval(parse(text=full_cmd)))
}
# nSingleton_lab2 <- get_fract_lab2(vect_values=tad_sing_dt$nSingleton, range_levels = nSingleton_levels)
# stopifnot(nSingleton_lab2 == tad_sing_dt$nSingleton_lab)


# tad_sing_dt$nSingleton_lab <- ifelse(tad_sing_dt$nSingleton == 0, "<=0",
#                                      ifelse(tad_sing_dt$nSingleton > 0 & tad_sing_dt$nSingleton <= 5, ">0 & <=5",
#                                             ifelse(tad_sing_dt$nSingleton > 5 & tad_sing_dt$nSingleton <= 10, ">5 & <=10",
#                                                    ifelse(tad_sing_dt$nSingleton > 10 & tad_sing_dt$nSingleton <= 15, ">10 & <=15",
#                                                           ifelse(tad_sing_dt$nSingleton > 10, ">15", NA)))))
# 
# stopifnot(!is.na(tad_sing_dt$nSingleton_lab))
# 
# nSingleton_levels <- c(0,5,10,15)

get_level_labs <- function(range_levels) {
  range_levels <- sort(range_levels)
  range_labels <- c(paste0("<=", range_levels[1]), paste0(">",range_levels[1:(length(range_levels)-1)], " & <=", range_levels[2:(length(range_levels))]),
                    paste0(">", range_levels[length(range_levels)]))
  return(range_labels)
}

get_fract_lab0 <- function(vect_values, range_levels) {
  stopifnot(is.numeric(vect_values))
  stopifnot(is.numeric(range_levels))
  range_levels <- sort(range_levels)
  range_labels <- c(paste0("<=", range_levels[1]), paste0(">",range_levels[1:(length(range_levels)-1)], " & <=", range_levels[2:(length(range_levels))]),
                    paste0(">", range_levels[length(range_levels)]))
  
  vect_labels <- sapply(vect_values, function(x) {
    tmpx <- hist(x, breaks=c(-Inf,range_levels, Inf),plot=F)$counts
    stopifnot(sum(tmpx) == 1)
    stopifnot(length(tmpx) == length(range_labels))
    range_labels[which(tmpx == 1)]
  })
  vect_labels <- factor(vect_labels, levels = range_labels)
  stopifnot(!is.na(vect_labels))
  stopifnot(table(vect_labels) == hist(vect_values, breaks=c(-Inf, range_levels, Inf), plot=FALSE)$counts)
  return(as.character(vect_labels))
}

do_densplot_withCorr <- function(xvar, yvar, plot_dt) {
  my_x <- plot_dt[,paste0(xvar)]
  my_y <- plot_dt[,paste0(yvar)]
  densplot(
    x=my_x,
    y=my_y,
    xlab=paste0(xvar),
    ylab=paste0(yvar),
    cex.main=plotCex,
    cex.axis=plotCex,
    cex.lab = plotCex,
    pch=16
  )
  addCorr(x=my_x, y=my_y, bty="n")
}


plot_density <- function(p) {
  p2 <- p+  
    guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(
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
  return(p2)
}
