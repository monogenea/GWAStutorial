# Genome-wide association (GWA) tutorial

## Additional files

For this tutorial you will additionally need the files

- 117Malay_282lipids.txt
- 120Indian_282lipids.txt
- 122Chinese_282lipids.txt
- 105Indian_2527458snps.bed, .bim, .fam
- 108Malay_2527458snps.bed, .bim, .fam
- 110Chinese_2527458snps.bed, .bim, .fam

stored in the folders 'Lipidomic' and 'Genomics' contained in the following compressed file:
https://sphfiles.nus.edu.sg/phg/Iomics/downloads/iOmics_data.tar.gz

**UPDATE 25/06/2019: uncompressed `iOmics_data.tar.gz` now directly available as `public/`**

I noticed the URL recently changed. To avoid problems with tracking the data, I have now hosted all of them in this repo. It is no longer necessary to download from the link above.

## Instructions

1. Combine the folders 'Lipidomic' and 'Genomics' and all files from this repo in your working directory.
2. Install all packages listed on top of the scripts. `snpStats` and `SNPRelate` are deposited in BioConductor, all other packages in CRAN.

**UPDATE 25/06/2019: Linux/macOS installation of GenABEL:**
```
install.packages("GenABEL.data", repos="http://R-Forge.R-project.org")
packageurl <- "https://cran.r-project.org/src/contrib/Archive/GenABEL/GenABEL_1.8-0.tar.gz"
install.packages(packageurl, repos=NULL)
```
3. Run the scripts in their exact numbered order.

## Acknowledgements

This work was largely based on the following publications:

- *Establishing multiple omics baselines for three Southeast Asian populations in the Singapore Integrative Omics Study*, Saw et al. (2017), Nat. Comm. (data source)
- *A guide to genome-wide association analysis and post-analytic interrogation*, Reed et al. (2015), Stats. in Med. (method source)

Also, thanks to @nizzle10, @rafalcode and @bambrozio for contributing. Enjoy, all feedback is welcome!

Francisco
