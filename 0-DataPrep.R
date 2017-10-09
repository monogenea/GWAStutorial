

library(snpStats)
load("conversionTable.RData")

pathM <- paste("Genomics/108Malay_2527458snps", c(".bed", ".bim", ".fam"), sep = "")
SNP_M <- read.plink(pathM[1], pathM[2], pathM[3])

pathI <- paste("Genomics/105Indian_2527458snps", c(".bed", ".bim", ".fam"), sep = "")
SNP_I <- read.plink(pathI[1], pathI[2], pathI[3])

pathC <- paste("Genomics/110Chinese_2527458snps", c(".bed", ".bim", ".fam"), sep = "")
SNP_C <- read.plink(pathC[1], pathC[2], pathC[3])

# Merge the three SNP datasets
SNP <- rbind(SNP_M$genotypes, SNP_I$genotypes, SNP_C$genotypes)

# Take one bim map (all 3 maps are based on the same ordered set of SNPs)
map <- SNP_M$map
colnames(map) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")

# Rename SNPs present in the conversion table into rs IDs
mappedSNPs <- intersect(map$SNP, names(conversionTable))
newIDs <- conversionTable[match(map$SNP[map$SNP %in% mappedSNPs], names(conversionTable))]
map$SNP[rownames(map) %in% mappedSNPs] <- newIDs

# Load lipid datasets & match SNP-Lipidomics samples
lipidsMalay <- read.delim("Lipidomic/117Malay_282lipids.txt", row.names = 1)
lipidsIndian <- read.delim("Lipidomic/120Indian_282lipids.txt", row.names = 1)
lipidsChinese <- read.delim("Lipidomic/122Chinese_282lipids.txt", row.names = 1)
all(Reduce(intersect, list(colnames(lipidsMalay),
                           colnames(lipidsIndian),
                           colnames(lipidsChinese))) == colnames(lipidsMalay)) # TRUE
lip <- rbind(lipidsMalay, lipidsIndian, lipidsChinese)

matchingSamples <- intersect(rownames(lip), rownames(SNP))
SNP <- SNP[matchingSamples,]
lip <- lip[matchingSamples,]

genData <- list(SNP = SNP, MAP = map, LIP = lip)
save(genData, file = "PhenoGenoMap.RData")

# Clear memory
rm(list = ls())
