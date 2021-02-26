setDir <- "/media/electron"
setDir <- ""
cell_line <- "Compendium_LG1"

chromo <- "chr21"

matfile <- file.path(setDir, "/mnt/ndata/Yuanlong/2.Results/1.Juicer",
                     cell_line, "contact_mat", paste0("mat_", chromo, "_40kb_ob.txt.gz"))

binSize <- 40000

mat_dt <- data.frame(
  coordA = binSize*c( rep(0, 7), rep(1,6), rep(2,5), rep(3,4), rep(4,3), rep(5,2), rep(6,1) ), 
  coordB = c( 0:6*binSize, 1:6*binSize, 2:6*binSize, 3:6*binSize, 4:6*binSize,  5:6*binSize, 6:6*binSize ),
  count = c(c(10,2,2,1,3,0,1), c(9,3,4,5,7,8), c(12,6,6,4,2), c(14,5,5,1), c(8,5,4), c(15,3), 11)
)

mat_dt$binA <- mat_dt$coordA/binSize
mat_dt$binB <- mat_dt$coordB/binSize
stopifnot(mat_dt$binA %% 1 == 0)
stopifnot(mat_dt$binB %% 1 == 0)
stopifnot(mat_dt$binA <= mat_dt$binB)

mat_dt$diagoDist <- mat_dt$binB-mat_dt$binA

mainDiagoSize <- max(mat_dt$binB)+1
mat_dt$diagoSize <- mainDiagoSize - mat_dt$diagoDist

by(mat_dt, mat_dt$diagoSize, function(x) nrow(x))



init_nrow <- nrow(mat_dt)
mat_dt <- na.omit(mat_dt)
cat(paste0("... after discarding NA values: ", nrow(mat_dt), "/", init_nrow, "\n"))


cat(paste0("... performing z-score normalization...\n"))

matNorm_dt <- do.call(rbind, by(mat_dt, mat_dt$diagoDist, function(x) {x$normCount <- as.numeric(scale(x$count)); x}))
# can produce Na if not enough value at one diagodist
matNorm_dt <- na.omit(matNorm_dt)

cat(paste0("... after discarding NA values: ", nrow(matNorm_dt), "/", nrow(mat_dt), "\n"))
  
# RAO BIN / BIN SIZE -> THIS GIVES THE 0-BASED COORD
require(Matrix)
mat_mat <- sparseMatrix(i = mat_dt$binA, 
             j = mat_dt$binB,
             x = mat_dt$count, index1=FALSE)
  
# 7 x 7 sparse Matrix of class "dgCMatrix"
# 
# [1,] 10 2  2  1 3  0  1
# [2,]  . 9  3  4 5  7  8
# [3,]  . . 12  6 6  4  2
# [4,]  . .  . 14 5  5  1
# [5,]  . .  .  . 8  5  4
# [6,]  . .  .  . . 15  3
# [7,]  . .  .  . .  . 11


tad_dt <- data.frame(
  region=c("tad1", "tad2", "tad3"),
  start = c(1,120001, 200001),
  end = c(120000,200000, 280000)
)
# convert to 0-based bin
tad_dt$startBin <- (tad_dt$start-1)/binSize
tad_dt$endBin <- (tad_dt$end)/binSize-1
stopifnot(tad_dt$startBin %% 1 == 0)
stopifnot(tad_dt$endBin %% 1 == 0)
stopifnot(tad_dt$startBin <= tad_dt$endBin)

stopifnot(nrow(tad_dt) > 2)

i=1
for(i in 1:c(nrow(tad_dt)-1)) {

  t_start <- tad_dt$startBin[i]
  t_end <- tad_dt$endBin[i]
  
  next_start <- tad_dt$startBin[i+1]
  next_end <- tad_dt$endBin[i+1]

  prev_start <- tad_dt$startBin[i-1]
  prev_end <- tad_dt$endBin[i-1]
  
  intra_values <- mat_dt$count[(t_start <= mat_dt$binA  & t_end >= mat_dt$binA) &
                                 (t_start <= mat_dt$binB  & t_end >= mat_dt$binB)]
  intra_values
  
  # store upper right corner, so binA in current, binB in next
  next_inter_values <- mat_dt$count[(t_start <= mat_dt$binA  & t_end >= mat_dt$binA) &
                                (next_start <= mat_dt$binB  & next_end >= mat_dt$binB)]
  next_inter_values
  
  # store upper right corner, so binA in previous, binB in current
  prev_inter_values <- mat_dt$count[(prev_start <= mat_dt$binA  & prev_end >= mat_dt$binA) &
                                      (t_start <= mat_dt$binB  & t_end >= mat_dt$binB)]
  prev_inter_values
  
  
  
}



  # read.csv(gzfile(matfile,'rt')  ,header=FALSE, sep="\t", col.names=c("coordA", "coordB", "count")) 

# # data from Rao et al. indicate bins by genome coordinates, they need to be turned into indeces
#   chr.data$binsA = chr.data$binsA/bin.size
#   chr.data$binsB = chr.data$binsB/bin.size
# }
# 
# chr.data$binsA = chr.data$binsA+1
# chr.data$binsB = chr.data$binsB+1