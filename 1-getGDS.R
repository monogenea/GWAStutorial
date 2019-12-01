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

# Merge the three SNP datasets
SNP_M$genotypes <- rbind(SNP_M$genotypes, SNP_I$genotypes, SNP_C$genotypes)
# Take one bim map (all 3 maps are based on the same ordered set of SNPs)
map <- SNP_M$map
colnames(map) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
SNP_M$fam<- rbind(SNP_M$fam, SNP_I$fam, SNP_C$fam)

# Rename SNPs present in the conversion table into rs IDs
mappedSNPs <- intersect(SNP_M$map$SNP, names(conversionTable))
newIDs <- conversionTable[match(SNP_M$map$SNP[SNP_M$map$SNP %in% mappedSNPs], names(conversionTable))]
SNP_M$map$SNP[rownames(SNP_M$map) %in% mappedSNPs] <- newIDs

write.plink("convertGDS", snps = SNP_M$genotypes)
