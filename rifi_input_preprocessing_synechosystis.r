##########################Synechosystis_chromosome#######################
setwd("/home/loub/Documents/Hess_lab/HESS_array/SC/")
load("fit1_SC.rdata")
time <- c(0, 2, 4, 8, 16, 32, 64)
chr <- which(fit1$genes$Replicon == "chromosome")
colnames(fit1$coefficients) <- as.numeric(time)
df <- cbind.data.frame(fit1$genes[chr, c("Replicon", "start", "end", 
                                         "orientation")], 
                       2^fit1$coefficients[chr,])
df <- df[order(df$start, decreasing = F),]
#change orientation as microarrays have the opposite strand orientation
str <- which(df$orientation == "+")
df$orientation[str] <- NA
df$orientation[-str] <- "+"
df$orientation[str] <- "-"
#eliminate those probes with same position at start or end
pos <- df %>%
    filter(orientation == "+") %>%
    distinct(end, .keep_all = TRUE)
neg <- df %>%
    filter(orientation == "-") %>%
    distinct(start, .keep_all = TRUE)
df <- rbind(pos, neg)
df$Replicon <- "chr"
df_rifi <- df

counts <- as.matrix(df[,5:ncol(df)])
rowRanges <- GRanges(seqnames = df[,"Replicon"],
                     IRanges(start = df[,"start"], 
                                     end = df[, "end"]),
                    strand = df[,"orientation"]
                    )
colData <- DataFrame(timepoint = time, replicate = 1,
                     row.names=time)
seqinfo <- data.frame(seqnames = "chr", seqlengths = 3573470,
                     isCircular = TRUE,  genome = "Mitschke2011")

syne_chr_SE <- SummarizedExperiment(assays=counts,
                           rowRanges=rowRanges, colData=colData)
setwd("/home/loub/Documents/Hess_lab/HESS_array/SC/se_syneco/")
save(syne_chr_SE, file="syne_chr_SE.rda")

##########################Synechosystis_chromosome#######################
setwd("/home/loub/Documents/Hess_lab/HESS_array/iron_depletion/")
load("fit1_iron.rdata")
time <- c(0, 2, 4, 8, 16, 32, 64)
chr <- which(fit1$genes$Replicon == "chromosome")
colnames(fit1$coefficients) <- as.numeric(time)
df <- cbind.data.frame(fit1$genes[chr, c("Replicon", "start", "end", 
                                         "orientation")], 
                       2^fit1$coefficients[chr,])
df <- df[order(df$start, decreasing = F),]
#change orientation as microarrays have the opposite strand orientation
str <- which(df$orientation == "+")
df$orientation[str] <- NA
df$orientation[-str] <- "+"
df$orientation[str] <- "-"
#eliminate those probes with same position at start or end
pos <- df %>%
    filter(orientation == "+") %>%
    distinct(end, .keep_all = TRUE)
neg <- df %>%
    filter(orientation == "-") %>%
    distinct(start, .keep_all = TRUE)
df <- rbind(pos, neg)
df$Replicon <- "chr"
df_rifi <- df
save(df_rifi, file="df_rifi.rda")
head(df_rifi)


counts <- as.matrix(df[,5:ncol(df)])
rowRanges <- GRanges(seqnames = df[,"Replicon"],
                     IRanges(start = df[,"start"], 
                             end = df[, "end"]),
                     strand = df[,"orientation"]
)
colData <- DataFrame(timepoint = time, replicate = 1,
                     row.names=time)
seqinfo <- data.frame(seqnames = "chr", seqlengths = 3573470,
                      isCircular = TRUE,  genome = "Mitschke2011")

syne_chr_SE_iron <- SummarizedExperiment(assays=counts,
                                    rowRanges=rowRanges, colData=colData)
setwd("/home/loub/Documents/Hess_lab/HESS_array/iron_depletion/se_iron")
save(syne_chr_SE_iron, file="syne_chr_SE_iron.rda")

