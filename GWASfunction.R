GWAA <- function(genodata=genotypes, phenodata=phenotypes, family = gaussian, filename=NULL, append=FALSE, workers=getOption("mc.cores",2L), flip=TRUE,
                 select.snps=NULL, hosts=NULL, nSplits=10){
      if (!require(doParallel)) { stop("Missing doParallel package") }
      #Check that a filename was specified
      if(is.null(filename)) stop("Must specify a filename for output.")
      #Check that the genotype data is of class 'SnpMatrix'
      if( class(genodata)!="SnpMatrix") stop("Genotype data must of class 'SnpMatrix'.")
      #Check that there is a variable named 'phenotype' in phenodata table
      if( !"phenotype" %in% colnames(phenodata)) stop("Phenotype data must have column named 'phenotype'")
      #Check that there is a variable named 'id' in phenodata table
      if( !"id" %in% colnames(phenodata)) stop("Phenotype data must have column named 'id'.")
      #If a vector of SNPs is given, subset genotype data for these SNPs
      if(!is.null(select.snps)) genodata<-genodata[,which(colnames(genodata)%in%select.snps)]
      #Check that there are still SNPs in 'SnpMatrix' object
      if(ncol(genodata)==0) stop("There are no SNPs in the 'SnpMatrix' object.")
      #Print the number of SNPs to be checked
      cat(paste(ncol(genodata), " SNPs included in analysis.\n"))
      #If append=FALSE than we will overwrite file with column names
      if(!isTRUE(append)) {
      columns<-c("SNP", "Estimate", "Std.Error", "t-value", "p-value")
      write.table(t(columns), filename, row.names=FALSE, col.names=FALSE, quote=FALSE) }
      # Check sample counts
      if (nrow(phenodata) != nrow(genodata)) {
            warning("Number of samples mismatch. Using subset found in phenodata.") }
      # Order genodata rows to be the same as phenodata
      genodata <- genodata[phenodata$id,]
      cat(nrow(genodata), "samples included in analysis.\n")
      # Change which allele is
      flip.matrix<-function(x){
            zero2 <- which(x==0)
            two0 <- which(x==2)
            x[zero2] <- 2
            x[two0] <- 0 
            return(x)
      }
      nSNPs <- ncol(genodata)
      genosplit <- ceiling(nSNPs/nSplits) # number of SNPs in each subset
      snp.start <- seq(1, nSNPs, genosplit) # index of first SNP in group
      snp.stop <- pmin(snp.start+genosplit-1, nSNPs) # index of last SNP in group
      if (is.null(hosts)) {
            # On Unix this will use fork and mclapply. On Windows it # will create multiple processes on localhost.
            cl <- makeCluster(workers)
      } else {
            # The listed hosts must be accessible by the current user using # password-less ssh with R installed on all hosts, all
            # packages installed, and "rscript" is in the default PATH.
            # See docs for makeCluster() for more information.
            cl <- makeCluster(hosts, "PSOCK")
      }
      show(cl) # report number of workers and type of parallel implementation
      registerDoParallel(cl)
      foreach (part=1:nSplits) %do% {
            # Returns a standar matrix of the
            genoNum <- as(genodata[,snp.start[part]:snp.stop[part]], "numeric")
            # Flip the numeric values of genotypes to count minor allele
            if (isTRUE(flip)) genoNum <- flip.matrix(genoNum)
            # For each SNP, concatenate the genotype column to the
            # phenodata and fit a generalized linear model
            rsVec <- colnames(genoNum)
            res <- foreach(snp.name=rsVec, .combine='rbind') %dopar% {
                  a <- summary(glm(phenotype~ . - id, family=family, data=cbind(phenodata, snp=genoNum[,snp.name])))
                  a$coefficients['snp',]
            }
            # write results so far to a file
            write.table(cbind(rsVec,res), filename, append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE)
            cat(sprintf("GWAS SNPs %s-%s (%s%% finished)\n", snp.start[part], snp.stop[part], 100*part/nSplits))
      }
      stopCluster(cl)
      return(print("Done."))
}


# ---- manhattan ----
# Receives a data.frame of SNPs with Neg_logP, chr, position, and type.
# Plots Manhattan plot with significant SNPs highlighted.
GWAS_Manhattan <- function(GWAS, col.snps=c("black","gray"),
                           col.detected=c("blue"), col.imputed=c("red"), col.text="black",
                           title="GWAS Tutorial Manhattan Plot", display.text=TRUE, bonferroni.alpha=0.05,
                           bonferroni.adjustment=1000000, Lstringent.adjustment=10000) {
      bonferroni.thresh <- -log10(bonferroni.alpha / bonferroni.adjustment)
      Lstringent.thresh <- -log10(bonferroni.alpha / Lstringent.adjustment)
      xscale <- 10000000
      manhat <- GWAS[!grepl("[A-z]",GWAS$chr),]
      #sort the data by chromosome and then location
      manhat.ord <- manhat[order(as.numeric(manhat$chr),manhat$position),]
      manhat.ord <- manhat.ord[!is.na(manhat.ord$position),]
      ##Finding the maximum position for each chromosome
      max.pos <- sapply(1:21, function(i) { max(manhat.ord$position[manhat.ord$chr==i],0) })
      max.pos2 <- c(0, cumsum(max.pos))
      #Add spacing between chromosomes
      max.pos2 <- max.pos2 + c(0:21) * xscale * 10
      #defining the positions of each snp in the plot
      manhat.ord$pos <- manhat.ord$position + max.pos2[as.numeric(manhat.ord$chr)]
      # alternate coloring of chromosomes
      manhat.ord$col <- col.snps[1 + as.numeric(manhat.ord$chr) %% 2]
      # draw the chromosome label roughly in the middle of each chromosome band
      text.pos <- sapply(c(1:22), function(i) { mean(manhat.ord$pos[manhat.ord$chr==i]) })
      # Plot the data
      plot(manhat.ord$pos[manhat.ord$type=="typed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="typed"],
           pch=20, cex=.3, col= manhat.ord$col[manhat.ord$type=="typed"], xlab="Chromosome", ylab="Negative Log P-Value", axes=F, ylim=c(0,max(manhat$Neg_logP)+1))
      points(manhat.ord$pos[manhat.ord$type=="imputed"]/xscale,
             manhat.ord$Neg_logP[manhat.ord$type=="imputed"], pch=20, cex=.4, col = col.imputed)
      points(manhat.ord$pos[manhat.ord$type=="typed"]/xscale,
             manhat.ord$Neg_logP[manhat.ord$type=="typed"], pch=20, cex=.3, col = manhat.ord$col[manhat.ord$type=="typed"])
      axis(2)
      abline(h=0)
      
      SigNifSNPs <- as.character(GWAS[GWAS$Neg_logP > Lstringent.thresh & GWAS$type=="typed", "SNP"])
      #Add legend
      #legend("topright",c("Bonferroni Corrected Threshold*", "Less Stringent Threshold**"),
      #       border="black", col=c("gray60", "gray60"), pch=c(0, 0), lwd=c(1,1), lty=c(1,2), pt.cex=c(0,0), bty="n", cex=0.7)
      #Add chromosome number
      text(text.pos/xscale, -.3, seq(1,22,by=1), xpd=TRUE, cex=1)
      #Add bonferroni line
      abline(h=bonferroni.thresh, untf = FALSE, col = "gray60")
      #Add "less stringent" line
      abline(h=Lstringent.thresh, untf = FALSE, col = "gray60", lty = 2)
      #Plotting detected genes #Were any genes detected?
      if (length(SigNifSNPs)>0){
            sig.snps <- manhat.ord[,'SNP'] %in% SigNifSNPs
            points(manhat.ord$pos[sig.snps]/xscale, manhat.ord$Neg_logP[sig.snps],
                   pch=15,col=col.detected, bg=col.detected,cex=0.5)
            text(manhat.ord$pos[sig.snps]/xscale,manhat.ord$Neg_logP[sig.snps],
                 as.character(manhat.ord[sig.snps,1]), col=col.text, offset=1, adj=-.1, cex=1)
      }
}

# ---- map2gene ----
# Returns the subset of SNPs that are within extend.boundary of gene
# using the coords table of gene locations.
map2gene <- function(gene, coords, SNPs, extend.boundary = 5000) {
      coordsSub <- coords[coords$gene == gene,] #Subset coordinate file for spcified gene
      coordsSub$start <- coordsSub$start - extend.boundary # Extend gene boundaries
      coordsSub$stop <- coordsSub$stop + extend.boundary
      SNPsub <- SNPs[SNPs$position >= coordsSub$start & SNPs$position <= coordsSub$stop &
                           SNPs$chr == coordsSub$chr,] #Subset for SNPs in gene
      return(data.frame(SNPsub, gene = gene, stringsAsFactors = FALSE))
}
      
      
      