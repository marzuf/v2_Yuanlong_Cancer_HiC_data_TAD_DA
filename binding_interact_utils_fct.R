
#*******************************************************************************************
#*******************************************************************************************
#*******************************************************************************************
# generate a random partition from a set of domains preserving
# the size of the domains and (if wanted) the pattern of domain/non-domain regions
# adapted from Duggal et al. 2013

shuffle_chromoPartition <- function(domainDT, chrSize, preservePattern=TRUE, setSeed = F , seed = 42) {
  if(setSeed) set.seed(seed)
  stopifnot(all(colnames(domainDT) == c("chromo", "start", "end")))
  chromo <- domainDT$chromo[1]
  ### STEP1: label all regions and associate interval with region label
  # prepare domain data
  domainDT$regType <- rep("domain", nrow(domainDT))
  # prepare gaps data
  gaps_start <- domainDT$end + 1
  gaps_start <- gaps_start[-length(gaps_start)]
  gaps_end <- domainDT$start - 1
  gaps_end <- gaps_end[-1]
  stopifnot(length(gaps_start) == length(gaps_end))
  gapsDT <- data.frame(chromo = rep(chromo, length(gaps_end)),
                       start= gaps_start,
                       end = gaps_end,
                       regType=rep("non_domain", length(gaps_end)))
  gapsDT <- gapsDT[gapsDT$end > gapsDT$start,]
  chromoDT <- rbind(domainDT, gapsDT)
  if(domainDT$start[1] > 1) {
    startDT <- data.frame(chromo = chromo,
                       start= 1,
                       end = domainDT$start[1]-1,
                       regType = "non_domain")
    chromoDT <- rbind(chromoDT, startDT)
  }
  if(domainDT$end[nrow(domainDT)] < chrSize) {
    endDT <- data.frame(chromo = chromo,
                          start= domainDT$end[nrow(domainDT)] + 1,
                          end = chrSize,
                        regType = "non_domain")
    chromoDT <- rbind(chromoDT, endDT)
  }
  chromoDT <- chromoDT[order(chromoDT$start),]
  chromoDT$len <- chromoDT$end - chromoDT$start + 1
  ### STEP2: generate a pair for each interval, O_list the list with all pairs in the observed ordering
  O_list <- setNames(chromoDT$len, chromoDT$regType) 
  ### STEP3: from O_list create two separate list for domains and non-domains
  O_list_dom <- O_list[grep("^domain$", names(O_list))]
  O_list_non <- O_list[grep("^non_domain$", names(O_list))]
  ### STEP4: randomly shuffle both lists
  # WARNING - undesirable behaviour if x is numeric and of length 1 !!!
  #           -> it does not sample but draws an integer from 1:x
  if(length(O_list_dom) > 1) {
    O_list_dom_rand <- sample(x=O_list_dom, size=length(O_list_dom))  
  } else{
    O_list_dom_rand <- O_list_dom
  }
  if(length(O_list_non) > 1) {
    O_list_non_rand <- sample(x=O_list_non, size=length(O_list_non))  
  } else{
    O_list_non_rand <- O_list_non
  }
  ### STEP5: traverse the observed sequence 0 of domains and non-domains
  ### if preservePattern == TRUE -> output succession of domain/non-domains as in the observed data
  if(!preservePattern) {
    if(length(O_list) > 1){
      O_list <- sample(x=O_list, size = length(O_list))  
    } else{
      O_list <- O_list
    }
  }
  dom_index <- non_index <- 1
  # !!! DO NOT USE dopar WHEN MODFIYING "GLOBAL" VARIABLE !!!
  randomPartition <- foreach(reg = names(O_list), .combine='c') %do% {
    if(reg == "domain") {
     z <- setNames(O_list_dom_rand[dom_index], "domain")
     dom_index <- dom_index + 1
    } else if(reg == "non_domain") {
      z <- setNames(O_list_non_rand[non_index], "non_domain")
      non_index <- non_index + 1
    } else {
		stop("error")
	}
    z
  }
  stopifnot(all(names(randomPartition) == names(O_list)))
  stopifnot( (dom_index - 1) == length(O_list_dom))
  stopifnot( (non_index - 1) == length(O_list_non))
  ### STEP6: reconstruct the domain list
  cum_randomPart <- cumsum(randomPartition)
  newChromoDT <- data.frame(chromo = rep(chromo, length(cum_randomPart)), 
                        start = c(1, cum_randomPart[-length(cum_randomPart)]+1 ),
                        end = as.numeric(cum_randomPart),
                        regType = names(cum_randomPart))
  stopifnot(nrow(newChromoDT) == nrow(chromoDT))
  newDomainDT <- newChromoDT[newChromoDT$regType == "domain",c("chromo", "start", "end")]
  stopifnot(nrow(newDomainDT) == nrow(domainDT))
 
  obsSum <- sum(domainDT$end-domainDT$start+1)
  shuffSum <- sum(newDomainDT$end-newDomainDT$start+1)
  if( abs(obsSum-shuffSum) > 1e-10) {
	cat(paste0( round(obsSum, 2) , " == ", round(shuffSum, 2), "\n"))
	stop("abs(obsSum-shuffSum) > 1e-10\n")
  }
  return(newDomainDT)
}

#**************************************************************************************
#**************************************************************************************
#**************************************************************************************
fill_DT <- function(dt, chr_len) {
    stopifnot(ncol(dt) == 3)
    chromo <- as.character(dt[1,1])
    # do not take the 1st column with the "chr6"    
    dt <- dt[,2:3]
    colnames(dt) <- c("start", "end")
    # add a column filled with 0
    dt$region <- "TAD"
    # test that nrow dt is bigger than 1 otherwise the 2:... will create NA
    if(nrow(dt) > 1) {
      # create data frame that will hold the gap for the 1st dataset
      # start of the gap = end of the end of the TAD + 1 (do not take the last row)    
      # end of the gap = start of the TAD - 1 (do not take the first row)
      dt_gaps <- data.frame( start = (dt$end[1:(nrow(dt)-1)] + 1),
                              end = (dt$start[2:nrow(dt)] -1))
      stopifnot(is.numeric(dt_gaps$start[1]))
      stopifnot(is.numeric(dt_gaps$end[1]))    
      # select only the row with end > start
      dt_gaps <- dt_gaps[dt_gaps$start < dt_gaps$end,]
    } else{
      dt_gaps <- data.frame(start = numeric(), end=numeric())
    }
    # ad gaps at the beginning until first TAD and at the end until end of chromo size
    # CHANGE MZ: FIRST DOMAIN APPENDED SHOULD START WITH 1 NOT WITH 0
    #pgaps1 = pgaps1.append(pd.DataFrame([[0, p1.iloc[0,0]-1], [p1.iloc[p1.shape[0]-1,1]+1, chr_len]], columns=['Start', 'End']), ignore_index=True)
    # THERE WAS A PROBLEM IF THE LAST TAD WAS UNTIL THE LAST CHROMO IT WHOULD HAD LAST ROW WITH CHR_LEN+1 CHR_LEN
    #dt_gaps = dt_gaps.append(pd.DataFrame([[1, dt.iloc[0,0]-1], [dt.iloc[dt.shape[0]-1,1]+1, chr_len]], columns=['Start', 'End']), ignore_index=True)
    # if needed, add gap before 1st TAD
    if(dt$start[1] > 1) {
        tmpDT <- data.frame(start = 1, end = dt$start[1]-1)
        dt_gaps <- rbind(tmpDT, dt_gaps)
    }
    # if needed, add gap until chromosome end    
    if(dt$end[nrow(dt)] < chr_len) {
        tmpDT <- data.frame(start = dt$end[nrow(dt)] + 1, end = chr_len)
        dt_gaps <- rbind(dt_gaps, tmpDT)        
    }
    # add a column to indicate there are gaps
    if(nrow(dt_gaps) > 0) {
        dt_gaps$region <- "gap"
        dt_final <- rbind(dt, dt_gaps)
    } else{
        dt_final <- dt
    }
    dt_final <- dt_final[order(dt_final$start),]
    rownames(dt_final) <- NULL
    dt_final$chromo <- chromo
    return(dt_final)
}
#**************************************************************************************
#**************************************************************************************
#**************************************************************************************


get_binding_interact_with_tad <- function(dt_para, dt_g2t, check_entrez=NULL) {
  
  
  bindinginteract_g2t_dt_tmp <- inner_join(dt_para, dt_g2t[,c("entrezID", "region")], by = c("entrezID_a" = "entrezID"))
  colnames(bindinginteract_g2t_dt_tmp)[colnames(bindinginteract_g2t_dt_tmp) == "region"] <- "region_a"
  bindinginteract_g2t_dt <- inner_join(bindinginteract_g2t_dt_tmp, dt_g2t[,c("entrezID", "region")], by = c("entrezID_b" = "entrezID"))
  colnames(bindinginteract_g2t_dt)[colnames(bindinginteract_g2t_dt) == "region"] <- "region_b"
  nrow(bindinginteract_g2t_dt)
  # 158274
  stopifnot(is.character(bindinginteract_g2t_dt$entrezID_a))
  stopifnot(is.character(bindinginteract_g2t_dt$entrezID_b))
  
  if(!is.null(check_entrez))
    stopifnot(setequal(check_entrez, unique(c(bindinginteract_g2t_dt$entrezID_b, bindinginteract_g2t_dt$entrezID_b))))
  
  
  sameTAD_ratio <- sum(bindinginteract_g2t_dt$region_a == bindinginteract_g2t_dt$region_b)/nrow(bindinginteract_g2t_dt)
  
  return(sameTAD_ratio)
  
}
