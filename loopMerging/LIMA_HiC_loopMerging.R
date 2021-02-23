
# INITIALIZE   -----------------------------------------------------------------------------

#install.packages('dbscan')
library(dbscan)
library(readr)



# FUNCTIONS   -----------------------------------------------------------------------------

# General function for rounding a number to a variable base
mround <- function(x,base, mode="std")
{ 
  if(mode == "std"){
    out = base*round(x/base) 
  }
  
  if(mode == "floor"){
    out = base*floor(x/base) 
  }
  
  if(mode == "ceiling"){
    out = base*ceiling(x/base) 
  }
  
  return(out)
}

# For surveying loop call quality (single hic/loop file)
makeHicPlot <- function(hic, loopList=diffLoops, type="box",
                        chr, str, end, res=10000,
                        X, Y, Z, W=2, H=2, justification=c("left", "top")){
  tempPLT = bb_plotHic(hic, chrom=chr, chromstart=str, chromend=end, zrange=c(0,Z), resolution=res,
                       x=unit(X,"inches"), y=unit(Y,"inches"), just=justification, 
                       width=unit(W,"inches"), height=unit(H,"inches"))
  
  if(!is.null(loopList) & nrow(loopList) != 0){
    bb_annotateLoops(hic=tempPLT, loops=loopList, type=type)
  }
  
  bb_labelGenome(tempPLT,x=unit(X,units="inches"), y=unit(Y+H,units="inches"), just=justification, scale="Mb")
}

# For converting SIP coordinates from 5kb --> 10kb resolution
convertCoords <- function(loopDF, newres=10000){
  newLoopDF = loopDF
  
  newLoopDF[,2] = mround(newLoopDF[,2], base=newres, mode="floor")
  newLoopDF[,3] = mround(newLoopDF[,3], base=newres, mode="ceiling")
  
  newLoopDF[,5] = mround(newLoopDF[,5], base=newres, mode="floor")
  newLoopDF[,6] = mround(newLoopDF[,6], base=newres, mode="ceiling")
  
  return(newLoopDF)
}

# For finding an "average" value based on mode, then mean
modeFinder <- function(x){
  ux <- unique(x)
  tabx <- tabulate(match(x,ux))
  
}

# For reporting an "average" loop based on position mode; assumes  BEDPE structure for col1-6
loopPicker <- function(loopDF, res=10000, strongestQuery="APScoreAvg") {
  # pull x1/y1 coordinates
  x <- loopDF[,2] 
  y <- loopDF[,5] 
  
  # unique x/y coords
  ux <- unique(x) 
  uy <- unique(y)  
  
  # counts per unique entry
  tabx <- tabulate(match(x, ux))
  taby <- tabulate(match(y, uy)) 
  
  # find mode, OR mean of modes (rounded down)
  valx <- mround(mean(ux[tabx == max(tabx)]), base=res, mode="floor") 
  valy <- mround(mean(uy[taby == max(taby)]), base=res, mode="floor")
  
  # find the strongest loop, based on strongestQuery, for filling in other loop info
  strongest = loopDF[which.max(loopDF[,strongestQuery]), ]
  
  # sum the loop counts (LIMA-specific!!!)
  callCounts = colSums(loopDF[,16:24])
  
  # build a "representative loop" from the strongest loop, adjusting coordinates to average
  repLoop <- strongest
  repLoop[,2:3] <- c(valx, valx+res)
  repLoop[,5:6] <- c(valy, valy+res)
  repLoop[,16:24] <- callCounts
  
  return(repLoop)
}



# READ IN     -----------------------------------------------------------------------------

# Read in loop calls for each TP (for merging)
loops0000 = read.table("./data/input/loops/perTP/SIP_newOpt2_MAPQ0/LIMA_THP1_WT_LPIF_0000_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g2_t2000.bedpe", header=T)
loops0030 = read.table("./data/input/loops/perTP/SIP_newOpt2_MAPQ0/LIMA_THP1_WT_LPIF_0030_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g2_t2000.bedpe", header=T)
loops0060 = read.table("./data/input/loops/perTP/SIP_newOpt2_MAPQ0/LIMA_THP1_WT_LPIF_0060_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g2_t2000.bedpe", header=T)
loops0090 = read.table("./data/input/loops/perTP/SIP_newOpt2_MAPQ0/LIMA_THP1_WT_LPIF_0090_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g2_t2000.bedpe", header=T)
loops0120 = read.table("./data/input/loops/perTP/SIP_newOpt2_MAPQ0/LIMA_THP1_WT_LPIF_0120_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g2_t2000.bedpe", header=T)
loops0240 = read.table("./data/input/loops/perTP/SIP_newOpt2_MAPQ0/LIMA_THP1_WT_LPIF_0240_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g2_t2000.bedpe", header=T)
loops0360 = read.table("./data/input/loops/perTP/SIP_newOpt2_MAPQ0/LIMA_THP1_WT_LPIF_0360_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g2_t2000.bedpe", header=T)
loops1440 = read.table("./data/input/loops/perTP/SIP_newOpt2_MAPQ0/LIMA_THP1_WT_LPIF_1440_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g2_t2000.bedpe", header=T)
loopOMEGA = read.table("./data/input/loops/omega_5kb_SIP_loops_firstOpt_MAPQ0.bedpe", header=T)

# loops0000 = read.table("./data/input/loops/perTP/SIP_newOpt1_MAPQ0/LIMA_THP1_WT_LPIF_0000_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loops0030 = read.table("./data/input/loops/perTP/SIP_newOpt1_MAPQ0/LIMA_THP1_WT_LPIF_0030_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loops0060 = read.table("./data/input/loops/perTP/SIP_newOpt1_MAPQ0/LIMA_THP1_WT_LPIF_0060_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loops0090 = read.table("./data/input/loops/perTP/SIP_newOpt1_MAPQ0/LIMA_THP1_WT_LPIF_0090_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loops0120 = read.table("./data/input/loops/perTP/SIP_newOpt1_MAPQ0/LIMA_THP1_WT_LPIF_0120_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loops0240 = read.table("./data/input/loops/perTP/SIP_newOpt1_MAPQ0/LIMA_THP1_WT_LPIF_0240_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loops0360 = read.table("./data/input/loops/perTP/SIP_newOpt1_MAPQ0/LIMA_THP1_WT_LPIF_0360_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loops1440 = read.table("./data/input/loops/perTP/SIP_newOpt1_MAPQ0/LIMA_THP1_WT_LPIF_1440_S_0.0.0_MAPQ0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loopOMEGA = read.table("./data/input/loops/omega_5kb_SIP_loops_newOpt1_MAPQ0.bedpe", header=T)

# loops0000 = read.table("./data/input/loops/perTP/SIP_newOpt1/LIMA_THP1_WT_LPIF_0000_S_0.0.0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loops0030 = read.table("./data/input/loops/perTP/SIP_newOpt1/LIMA_THP1_WT_LPIF_0030_S_0.0.0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loops0060 = read.table("./data/input/loops/perTP/SIP_newOpt1/LIMA_THP1_WT_LPIF_0060_S_0.0.0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loops0090 = read.table("./data/input/loops/perTP/SIP_newOpt1/LIMA_THP1_WT_LPIF_0090_S_0.0.0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loops0120 = read.table("./data/input/loops/perTP/SIP_newOpt1/LIMA_THP1_WT_LPIF_0120_S_0.0.0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loops0240 = read.table("./data/input/loops/perTP/SIP_newOpt1/LIMA_THP1_WT_LPIF_0240_S_0.0.0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loops0360 = read.table("./data/input/loops/perTP/SIP_newOpt1/LIMA_THP1_WT_LPIF_0360_S_0.0.0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loops1440 = read.table("./data/input/loops/perTP/SIP_newOpt1/LIMA_THP1_WT_LPIF_1440_S_0.0.0_SIP_loops_5kb_f05_g1_t2500.bedpe", header=T)
# loopOMEGA = read.table("./data/input/loops/omega_5kb_SIP_loops_newOpt1.bedpe", header=T)

# loops0000 = read.table("./data/input/loops/perTP/hiccups/LIMA_THP1_WT_LPIF_0000_S_0.0.0_megaMap_postprocessed_pixels_5000.bedpe", header=F)
# loops0030 = read.table("./data/input/loops/perTP/hiccups/LIMA_THP1_WT_LPIF_0000_S_0.0.0_megaMap_postprocessed_pixels_5000.bedpe", header=F)
# loops0060 = read.table("./data/input/loops/perTP/hiccups/LIMA_THP1_WT_LPIF_0000_S_0.0.0_megaMap_postprocessed_pixels_5000.bedpe", header=F)
# loops0090 = read.table("./data/input/loops/perTP/hiccups/LIMA_THP1_WT_LPIF_0000_S_0.0.0_megaMap_postprocessed_pixels_5000.bedpe", header=F)
# loops0120 = read.table("./data/input/loops/perTP/hiccups/LIMA_THP1_WT_LPIF_0000_S_0.0.0_megaMap_postprocessed_pixels_5000.bedpe", header=F)
# loops0240 = read.table("./data/input/loops/perTP/hiccups/LIMA_THP1_WT_LPIF_0000_S_0.0.0_megaMap_postprocessed_pixels_5000.bedpe", header=F)
# loops0360 = read.table("./data/input/loops/perTP/hiccups/LIMA_THP1_WT_LPIF_0000_S_0.0.0_megaMap_postprocessed_pixels_5000.bedpe", header=F)
# loops1440 = read.table("./data/input/loops/perTP/hiccups/LIMA_THP1_WT_LPIF_0000_S_0.0.0_megaMap_postprocessed_pixels_5000.bedpe", header=F)
# loopOMEGA = read.table("./data/input/loops/LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_30_withNorms_HiCCUPS_postprocessed_pixels_5000.bedpe", header=T)


# Convert to 10kb coordinates
loops0000 = convertCoords(loops0000)
loops0030 = convertCoords(loops0030)
loops0060 = convertCoords(loops0060)
loops0090 = convertCoords(loops0090)
loops0120 = convertCoords(loops0120)
loops0240 = convertCoords(loops0240)
loops0360 = convertCoords(loops0360)
loops1440 = convertCoords(loops1440)
loopOMEGA = convertCoords(loopOMEGA)

# Add in TP counter columns
counterRow = data.frame(
  "call0000" = 0, "call0030" = 0, "call0060" = 0, "call0090" = 0,
  "call0120" = 0, "call0240" = 0, "call0360" = 0, "call1440" = 0, "callOMGA" = 0
)
addCounter <- function(loopDF, i){
  rowN = nrow(loopDF)
  colN = ncol(loopDF)
  
  loopDF = cbind(loopDF, counterRow[rep(1, times=rowN),])
  loopDF[,colN+i] = 1
  return(loopDF)
}

loopDFlist = list(loops0000, loops0030, loops0060, loops0090,
                  loops0120, loops0240, loops0360, loops1440, loopOMEGA)
for (n in 1:length(loopDFlist)){
  loopDFlist[[n]] = addCounter(loopDFlist[[n]], n)
  n=n+1
}

# Read in hic files for each TP (for surveys)
hic0000="/Volumes/KatieEHD/Data/LIMA_share/hic/hic_maps/LIMA_THP1_WT_LPIF_0000_S_0.0.0_megaMap_inter_30.hic"
hic0030="/Volumes/KatieEHD/Data/LIMA_share/hic/hic_maps/LIMA_THP1_WT_LPIF_0030_S_0.0.0_megaMap_inter_30.hic"
hic0060="/Volumes/KatieEHD/Data/LIMA_share/hic/hic_maps/LIMA_THP1_WT_LPIF_0060_S_0.0.0_megaMap_inter_30.hic"
hic0090="/Volumes/KatieEHD/Data/LIMA_share/hic/hic_maps/LIMA_THP1_WT_LPIF_0090_S_0.0.0_megaMap_inter_30.hic"
hic0120="/Volumes/KatieEHD/Data/LIMA_share/hic/hic_maps/LIMA_THP1_WT_LPIF_0120_S_0.0.0_megaMap_inter_30.hic"
hic0240="/Volumes/KatieEHD/Data/LIMA_share/hic/hic_maps/LIMA_THP1_WT_LPIF_0240_S_0.0.0_megaMap_inter_30.hic"
hic0360="/Volumes/KatieEHD/Data/LIMA_share/hic/hic_maps/LIMA_THP1_WT_LPIF_0360_S_0.0.0_megaMap_inter_30.hic"
hic1440="/Volumes/KatieEHD/Data/LIMA_share/hic/hic_maps/LIMA_THP1_WT_LPIF_1440_S_0.0.0_megaMap_inter_30.hic"
hicOMGA="/Volumes/KatieEHD/Data/LIMA_share/hic/hic_maps/LIMA_THP1_WT_LPIF_omega_S_0.0.0_megaMap_inter_30.hic"

# Format loop calls with loop names
rownames(loopDFlist[[1]]) = paste0("loop0000_", 1:nrow(loopDFlist[[1]]))
rownames(loopDFlist[[2]]) = paste0("loop0030_", 1:nrow(loopDFlist[[2]]))
rownames(loopDFlist[[3]]) = paste0("loop0060_", 1:nrow(loopDFlist[[3]]))
rownames(loopDFlist[[4]]) = paste0("loop0090_", 1:nrow(loopDFlist[[4]]))
rownames(loopDFlist[[5]]) = paste0("loop0120_", 1:nrow(loopDFlist[[5]]))
rownames(loopDFlist[[6]]) = paste0("loop0240_", 1:nrow(loopDFlist[[6]]))
rownames(loopDFlist[[7]]) = paste0("loop0360_", 1:nrow(loopDFlist[[7]]))
rownames(loopDFlist[[8]]) = paste0("loop1440_", 1:nrow(loopDFlist[[8]]))
rownames(loopDFlist[[9]]) = paste0("loopOMGA_", 1:nrow(loopDFlist[[9]]))



# RUN    -----------------------------------------------------------------------------

# ** Concept ===============
# Combine all loops together (simple rbind, no merging; duplicates possible)
concatLoops = rbind(loopDFlist[[1]], loopDFlist[[2]], loopDFlist[[3]], loopDFlist[[4]], 
                    loopDFlist[[5]], loopDFlist[[6]], loopDFlist[[7]], loopDFlist[[8]], 
                    loopDFlist[[9]])

# Initialize
mergedLoops = c()
i=0

# For each chromosome....
for (chr in levels(concatLoops[,1])){
  i=i+1
  
  # ... pull all loops on that chromosome, extract 2D coords and cluster (NOTE: can change epsilon)
  chrLoops = concatLoops[concatLoops[,1] == chr,]
  coords = chrLoops[, c(2,5)]
  mdist = dist(coords, method="manhattan")
  clustr = dbscan(mdist, eps=20000, minPts = 2)
  
  # ... identify loops that were not clustered (no merging); should be few
  singleLoops = chrLoops[which(clustr$cluster == 0), ] # you could add on columns here and in the loop to note things like # over TP, TP of max, etc
  
  # ... for the loops that are clustered, find the strongest (based on observed counts [V12] here) and save their coords
  overlapLoops = c()
  if (max(clustr$cluster) != 0){ # ... if there are any clusters...
    for (c in 1: max(clustr$cluster)){  
      clustLoops = chrLoops[which(clustr$cluster == c), ]
      # strongest = clustLoops[which.max(clustLoops[,12]), ]
      # strongest = chrLoops[which(clustr$cluster == c)[which.max(chrLoops[which(clustr$cluster == c), 12])], ]  ## one-liner version
      
      overlapLoops[[c]] <- loopPicker(clustLoops)
    }
    # ... combine all the merged loops with the single loops from before 
    allOverlapLoops = do.call(rbind, overlapLoops)
    allChrLoops = rbind(singleLoops, allOverlapLoops)
  } else{
    allChrLoops = singleLoops # ... or don't, if there weren't any clusters (i.e. chrY)
  }
  
  # ... order by start position
  allChrLoops = allChrLoops[order(allChrLoops[,2]),]
  
  # ... add on to the list of merged loops for this chromosome (i)
  mergedLoops[[i]] = allChrLoops
}

# ... combine loops from every chromosome into one data frame
allMergedLoops = do.call(rbind, mergedLoops)

# Compare with original list
nrow(allMergedLoops)
nrow(concatLoops)

# Swap coordinates where y1 became smaller than x1
swapped = (allMergedLoops$x1 > allMergedLoops$y1) # n=26
allMergedLoops[swapped,1:6] = allMergedLoops[swapped,c(4:6, 1:3)]

# Write the loop list to a file for extracting counts --> run in DESeq
#write_tsv(allMergedLoops, "./data/input/loops/LIMA_SIP_10kbE_mergedLoops_eps20kMH.bedpe")
#write_tsv(allMergedLoops, "./data/input/loops/LIMA_SIP_MAPQ0_10kbE_mergedLoops_eps20kMH.bedpe")
#write_tsv(allMergedLoops, "./data/input/loops/LIMA_SIP_MAPQ0_10kbE_mergedLoops_conserv_eps20kMH.bedpe")

# Use this BEDPE file to pull counts on cluster!


# ** Notes/Testing ===============

# Understanding Epsilon
test = data.frame(x=rep(seq(100000, 125000, 5000),times=2), y=c(rep(200000, times=6,), seq(200000, 225000, 5000)))
dbscan(test, eps=20000, minPts=2) # 1 cluster (12)
dbscan(test, eps=10000, minPts=2) # 1 cluster (12)
dbscan(test, eps=7072, minPts=2) # 1 cluster (12)
dbscan(test, eps=5000, minPts=2)  # 1 cluster (8) + 4 outside

dist(test)

# You can use manhattan distance instead of euclidean if you make a dist object first!
mdist = dist(test, method="manhattan")
dbscan(mdist, eps=20000, minPts=2)  # 1 cluster (12)
dbscan(mdist, eps=10000, minPts=2)  # 1 cluster (12)
dbscan(mdist, eps=8000, minPts=2)  # 1 cluster (8) + 4 outside
dbscan(mdist, eps=5000, minPts=2)  # 1 cluster (8) + 4 outside
mdist


test = data.frame(x=rep(seq(100000, 125000, 5000), times=6), y=rep(seq(100000, 125000, 5000), each=6))
dbscan(test, eps=20000, minPts=2) # 1 cluster (36)
dbscan(test, eps=10000, minPts=2) # 1 cluster (36)
dbscan(test, eps=5000, minPts=2) # 1 cluster (36)


