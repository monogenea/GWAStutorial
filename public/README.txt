#====================================================================================================
#
#         FILE:  README.txt
#
#  DESCRIPTION:  README for resources data of Singapore Integrative Omics Study:
#                Establishing multiple omics baselines for three Southeast Asian populations
#
#=====================================================================================================

ABOUT OUR STUDY:
================
The study aims to establish baselines for each of the three ethnic groups in Singapore
- Chinese, Malays and Indians across multiple omics platform and to interrogate the 
extent of ethnic differences in each omic measurement. This iOmics can potential be used 
to investigate the degree of co-expression that exists between the different omics 
measurements. 

**Details of the quality controls for each omic platform were described in the main paper and
supplementary materials.**

FILES IN EACH FOLDERS:
======================
Info:
1. Sample QC status file across different platform 
   	- iomics_ID.csv

Genomic: 
1. Merging of genotyping file between Illumina 2.5M (2,299,708 SNPs) & Exome chips (227,750 SNPs)
	- 110Chinese_2527458snps (.bed/.bim/.fam)
	- 108Malay_2527458snps (.bed/.bim/.fam)
	- 105Indian_2527458snps (.bed/.bim/.fam)
2. Pharmacogenomic variants (4,032 SNPs)
	- 106Chinese_4032snps_pharmacogenomics (.bed/.bim/.fam)
	- 112Malay_4032snps_pharmacogenomics (.bed/.bim/.fam)
	- 115Indian_4032snps_pharmacogenomics (.bed/.bim/.fam)
3.  4-digits HLA typing alleles at eight HLA loci (198 HLA alleles)
	- 111Chinese_HLA.ped
	- 119Malay_HLA.ped
	- 120Indian_HLA.ped

Transcriptomic:
1.	21,649 probesets intensity text files
- 98Chinese_21649probesets.txt
- 75Malay_21649probesets.txt
- 96Indian_21649probesets.txt

Lipidomic:
1.	282 lipid species from four  lipid categories (Glycerophospholipids, Sphingolipids, Sterols and Glycerolipids) relative 
concentraion (in log2 pmol/ml) text files 
- 122Chinese_282lipids.txt
- 117Malay_282lipids.txt
- 120Indian_282lipids.txt

miRNA:
1.	274 miRNAs expression text files
- 117Chinese_274miRNAs.txt
- 115Malay_274miRNAs.txt
- 119Indian_274miRNAs.txt

Phenotype:
1. Questionnaire template (Lifestyle, Health survey and food frequency questionnaire)
	- IOmics_Questionnaire_and_FFQ.docx
	- IOmics_Health_Screening_form.docx
2. Data after QC ( 46 lifestyle variables, 39 clinical variables, 199 dietary variables) & Cookbook for the questionnaire and data
	- IOmics_Questionnaire_HealthScreening_and_FFQ.xlsx or IOmics_Questionnaire_HealthScreening_and_FFQ.csv
	- To note: 
		a) The uploaded data contains missingness (ie. "NA").  But in our downstream analysis, we replaced 
		them with the mean value
		b) The dietary data unit (it is a semi-quantitative food frequency questionnaire) has been converted 
		to number of servings eaten per month by multiplying number of times eaten with type of period.
	
		

DATA RELEASE POLICY:
=====================
Please cite the publication if you are using the data in any publication.
To access the raw data, please contact statyy@nus.edu.sg

Funding agencies/Acknowledgements:
===================================
This project acknowledges the support of Saw Swee Hock School of Public Health, the Yong Loo Lin School of Medicine, The National University Health System and the Life Sciences Institute. W.Y.S and Y.Y.T additionally acknowledge support from the Biomedical Research Council (grant 03/1/27/18/216), National Medical Research Council (grant 0838/2004), National Research Foundation (through the Biomedical Resarch Council, grants 05/1/21/19/425 and 11/1/21/19/678).
