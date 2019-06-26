#!/usr/bin/env Rscript
# 
for(pkg in c("snpStats", "doParallel", "SNPRelate", "GenABEL")){
      if(!require(pkg, character.only = T)) {
            stop(pgk, " is required for this script. Please install it on your system.")
      }
}

source("GWASfunction.R")
load("PhenoGenoMap.RData")

# Use SNP call rate of 100%, MAF of 0.1 (very stringent)
maf <- 0.1
callRate <- 1
SNPstats <- col.summary(genData$SNP)

maf_call <- with(SNPstats, MAF > maf & Call.rate == callRate)
genData$SNP <- genData$SNP[,maf_call]
genData$MAP <- genData$MAP[maf_call,]
SNPstats <- SNPstats[maf_call,]

# Sample call rate & heterozygosity
callMat <- !is.na(genData$SNP)
Sampstats <- row.summary(genData$SNP)
hetExp <- callMat %*% (2 * SNPstats$MAF * (1 - SNPstats$MAF)) # Hardy-Weinberg heterozygosity (expected)
hetObs <- with(Sampstats, Heterozygosity * (ncol(genData$SNP)) * Call.rate)
Sampstats$hetF <- 1-(hetObs/hetExp)
# Use sample call rate of 100%, het threshold of 0.1 (very stringent)
het <- 0.1 # Set cutoff for inbreeding coefficient;
het_call <- with(Sampstats, abs(hetF) < het & Call.rate == 1)
genData$SNP <- genData$SNP[het_call,]
genData$LIP <- genData$LIP[het_call,]

# LD and kinship coeff
ld <- .2
kin <- .1
snpgdsBED2GDS(bed.fn = "convertGDS.bed", bim.fn = "convertGDS.bim",
              fam.fn = "convertGDS.fam", out.gdsfn = "myGDS",
              cvt.chr = "char")
genofile <- snpgdsOpen("myGDS", readonly = F)
gds.ids <- read.gdsn(index.gdsn(genofile,  "sample.id"))
gds.ids <- sub("-1", "", gds.ids)
add.gdsn(genofile, "sample.id", gds.ids, replace = T)
geno.sample.ids <- rownames(genData$SNP)
# First filter for LD
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld,
                          sample.id = geno.sample.ids,
                          snp.id = colnames(genData$SNP))
snpset.ibd <- unlist(snpSUB, use.names = F)
# And now filter for MoM
ibd <- snpgdsIBDMoM(genofile, kinship = T,
                    sample.id = geno.sample.ids,
                    snp.id = snpset.ibd,
                    num.thread = 1)
ibdcoef <- snpgdsIBDSelection(ibd)
ibdcoef <- ibdcoef[ibdcoef$kinship >= kin,]

# Filter samples out
related.samples <- NULL
while (nrow(ibdcoef) > 0) {
      # count the number of occurrences of each and take the top one
      sample.counts <- sort(table(c(ibdcoef$ID1, ibdcoef$ID2)), decreasing = T)
      rm.sample <- names(sample.counts)[1]
      cat("Removing sample", rm.sample, "too closely related to",
          sample.counts[1], "other samples.\n")
      
      # remove from ibdcoef and add to list
      ibdcoef <- ibdcoef[ibdcoef$ID1 != rm.sample & ibdcoef$ID2 != rm.sample,]
      related.samples <- c(as.character(rm.sample), related.samples)
}
genData$SNP <- genData$SNP[!(rownames(genData$SNP) %in% related.samples),]
genData$LIP <- genData$LIP[!(rownames(genData$LIP) %in% related.samples),]

# PCA
set.seed(100)
pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids, 
                 snp.id = snpset.ibd, num.thread = 1)
pctab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1],
                    PC2 = pca$eigenvect[,2],
                    stringsAsFactors = F)

# Subset and/or reorder origin accordingly
origin <- origin[match(pca$sample.id, origin$sample.id),]

pcaCol <- rep(rgb(0,0,0,.3), length(pca$sample.id)) # Set black for chinese
pcaCol[origin$Country == "I"] <- rgb(1,0,0,.3) # red for indian
pcaCol[origin$Country == "M"] <- rgb(0,.7,0,.3) # green for malay

png("PCApopulation.png", width = 500, height = 500)
plot(pctab$PC1, pctab$PC2, xlab = "PC1", ylab = "PC2", col = pcaCol, pch = 16)
abline(h = 0, v = 0, lty = 2, col = "grey")
legend("top", legend = c("Chinese", "Indian", "Malay"), col = 1:3, pch = 16, bty = "n")
dev.off()

# Choose trait for association analysis, use colnames(genData$LIP) for listing
# NOTE: Ignore the first column of genData$LIP (gender)
target <- "Cholesterol"

phenodata <- data.frame("id" = rownames(genData$LIP),
                        "phenotype" = scale(genData$LIP[,target]), stringsAsFactors = F)

# Conduct GWAS (will take a while)
start <- Sys.time()
GWAA(genodata = genData$SNP, phenodata = phenodata, filename = paste(target, ".txt", sep = ""))
Sys.time() - start # benchmark

# Manhattan plot
GWASout <- read.table(paste(target, ".txt", sep = ""), header = T, colClasses = c("character", rep("numeric",4)))
GWASout$type <- rep("typed", nrow(GWASout))
GWASout$Neg_logP <- -log10(GWASout$p.value)
GWASout <- merge(GWASout, genData$MAP[,c("SNP", "chr", "position")])
GWASout <- GWASout[order(GWASout$Neg_logP, decreasing = T),]

png(paste(target, ".png", sep = ""), height = 500,width = 1000)
GWAS_Manhattan(GWASout)
dev.off()

# QQ plot using GenABEL estlambda function
png(paste(target, "_QQplot.png", sep = ""), width = 500, height = 500)
lambda <- estlambda(GWASout$t.value**2, plot = T, method = "median")
dev.off()

writeLines(capture.output(sessionInfo()), "sessionInfo")