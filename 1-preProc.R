#!/usr/bin/env Rscript
# 
if(!require("snpStats")) {
      stop("snpStats is required for this script. Please install it on your system.")
}
load("conversionTable.RData")

pathM <- paste("public/Genomics/108Malay_2527458snps", c(".bed", ".bim", ".fam"), sep = "")
SNP_M <- read.plink(pathM[1], pathM[2], pathM[3])

pathI <- paste("public/Genomics/105Indian_2527458snps", c(".bed", ".bim", ".fam"), sep = "")
SNP_I <- read.plink(pathI[1], pathI[2], pathI[3])

pathC <- paste("public/Genomics/110Chinese_2527458snps", c(".bed", ".bim", ".fam"), sep = "")
SNP_C <- read.plink(pathC[1], pathC[2], pathC[3])

# Ensure == number of markers across the three populations
if(ncol(SNP_C$genotypes) != ncol(SNP_I$genotypes)){
        stop("Different number of columns in input files detected. This is not allowed.")
}
if(ncol(SNP_I$genotypes) != ncol(SNP_M$genotypes)){
        stop("Different number of columns in input files detected. This is not allowed.")
}

# Merge the three SNP datasets
SNP <- SNP_M
SNP$genotypes <- rbind(SNP_M$genotypes, SNP_I$genotypes, SNP_C$genotypes)
colnames(SNP$map) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2") # same for all three
SNP$fam<- rbind(SNP_M$fam, SNP_I$fam, SNP_C$fam)

# Rename SNPs present in the conversion table into rs IDs
mappedSNPs <- intersect(SNP$map$SNP, names(conversionTable))
newIDs <- conversionTable[match(SNP$map$SNP[SNP$map$SNP %in% mappedSNPs], names(conversionTable))]
SNP$map$SNP[rownames(SNP$map) %in% mappedSNPs] <- newIDs

# Load lipid datasets & match SNP-Lipidomics samples
lipidsMalay <- read.delim("public/Lipidomic/117Malay_282lipids.txt", row.names = 1)
lipidsIndian <- read.delim("public/Lipidomic/120Indian_282lipids.txt", row.names = 1)
lipidsChinese <- read.delim("public/Lipidomic/122Chinese_282lipids.txt", row.names = 1)

all(Reduce(intersect, list(colnames(lipidsMalay),
                           colnames(lipidsIndian),
                           colnames(lipidsChinese))) == colnames(lipidsMalay)) # TRUE
lip <- rbind(lipidsMalay, lipidsIndian, lipidsChinese)

# Country
country <- sapply(list(SNP_M, SNP_I, SNP_C), function(k){
        nrow(k$genotypes)
})
origin <- data.frame(sample.id = rownames(SNP$genotypes),
                     Country = factor(rep(c("M", "I", "C"), country)))

matchingSamples <- intersect(rownames(lip), rownames(SNP$genotypes))
SNP$genotypes <- SNP$genotypes[matchingSamples,]
lip <- lip[matchingSamples,]
origin <- origin[match(matchingSamples, origin$sample.id),]
# Combine SNP and Lipidomics
genData <- list(SNP = SNP$genotype, MAP = SNP$map, LIP = lip)

# Write processed omics and GDS
save(genData, origin, file = "PhenoGenoMap.RData")
write.plink("convertGDS", snps = SNP$genotypes)

# Clear memory
rm(list = ls())
