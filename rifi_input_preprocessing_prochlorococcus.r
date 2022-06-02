##########################Prochlorococcus_chromosome#######################
setwd("/home/loub/microarrays_species/steglish_2010/fromClaudia/Prochlorococcus_single_probes")

load("y2.rdata")
setwd("/home/loub/microarrays_species/steglish_2010")
load("Data_pro_minus60.rdata")
setwd("/home/loub/microarrays_species/steglish_2010/proc_se/")


strandSelection <- function(data, Str){
    df <- data %>%
        filter(strand==Str)
    df <- df[order(df[,"position"], decreasing = F),] 
    return(df)
}

daf <- y2[,c(9:11)]

colnames(daf) <- c("start","position", "strand")
rownames(Data) <- NULL
daf$start <- convertTonumeric(daf, "start")
daf$position <- convertTonumeric(daf, "position")
df1 <- strandSelection(daf, "+")
df2 <- strandSelection(daf, "-")
df_joined <- rbind(df1, df2)#32068
df_joined <- right_join(df_joined, Data, by=c("position", "strand"))
df <- df_joined
convertTonumeric <- function(data, x){
    x <- as.numeric(as.character(data[,x]))
    return(x)
}

df$start <- convertTonumeric(df, "start")
df$position <- convertTonumeric(df, "position")

time <- as.numeric(c(0, 2.5, 5, 10, 20, 40))
df <- df_joined
counts <- as.matrix(df[, c(4:9)])#60 time point is excluded due to background effect
rowRanges <- GRanges(seqnames = "chr",
                     IRanges(start = (df[, "start"]), 
                             end = (df[, "position"])),
                     strand = df[,"strand"]
)
colData <- DataFrame(timepoint = time, replicate = 1)
# rowRanges$seqinfo <- data.frame(seqnames = "chr", seqlengths = 3573470,isCircular = TRUE,  genome = "Mitschke2011")

proch_SE <- SummarizedExperiment(assays = counts,
                                 rowRanges = rowRanges,
                                 colData = colData)

save(proch_SE, file="proch_SE.rda")