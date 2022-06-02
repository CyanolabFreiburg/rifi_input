
##########################Synechococcus_chromosome#######################
setwd("/home/loub/microarrays_species/synechococcus/")

load("input_df.RData")

input_df$start <- input_df$position - 50

time <- as.numeric(colnames(input_df[1:6]))
counts <- as.matrix(input_df[,1:6])#60 time point is excluded due to background effect
rowRanges <- GRanges(seqnames = "chr",
                     IRanges(start = as.numeric(input_df[,"start"]), 
                             end = as.numeric(input_df[,"position"])),
                     strand = input_df[,"strand"]
)
colData <- DataFrame(timepoint = time, replicate = 1)
synechoco_SE <- SummarizedExperiment(assays = counts,
                                 rowRanges = rowRanges, colData = colData)

save(synechoco_SE, file="synechoco_SE.rda")
