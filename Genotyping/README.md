# Computational Genotyping of Mouse Embryos

## Preprocessing
The Cell Ranger pipeline version 3 (10x Genomics Inc.) was used for each scRNA-seq data set to de-multiplex the raw base call files, generate the fastq files, perform the alignment against the mouse reference genome mm10, filter the alignment and count barcodes and UMIs. Outputs from multiple sequencing runs were also combined using Cell Ranger.

## Alignment (alignment_snpsplit.sh)
For each experiment, the scRNA-seq data were aligned against an mm10 hybrid mouse genome assembly using STAR with default settings and "--outSAMattributes NH HI NM MD." The hybrid genome was prepared using SNPsplit to mask SNPs between the mouse version mm10 (GRCm38) and the CAST/EiJ strain genomes with the ambiguity base (N). Subsequently, SNPsplit was used to sort reads that cover SNPs by origin (reference genome). Unambiguous and unique alignments of WT samples were used to create a list of SNPs that were covered by reads originating from both reference genomes. Finally, reads covering these SNPs were used to determine the allele composition for each cell, i.e. fraction of CAST/EiJ specific SNPs.

## Cell to embryo assignment (cell2embryo.R)
Single cells were assigned to embryos according to the autosomal fraction of CAST SNPs, a 19-dimension vector that allowed us to estimate the number of embryos per experiment. A minimum number of 1,000 covered SNPs and SNP information for each autosome was required. k-means clustering for multiple k (kmeans function in R, k = 2-15, default parameters) was performed on cells that fulfilled this criterion and evaluated by calculating the AIC for each model. The k with the minimal AIC defined the number of detected embryos, and the kernel averages represent the SNP profile for each embryo in the pool. Cells were then assigned to the embryo based on minimum distance in their SNP profile.

## Embryo sex determination (sex2embryo.R)
Embryo sex was determined based on the expression of the following genes: Xist (ENSMUSG00000086503) to count XX contexts and Erdr1 (ENSMUSG00000096768), Ddx3y (ENSMUSG00000069045) and Eif2s3y (ENSMUSG00000069049) to reliably detect transcription from the Y chromosome. The Cell Ranger gene barcode matrices were used to obtain per cell expression counts for these 4 genes and determine the fraction of positive cells per embryo. Embryos with a high percentage of Xist expressing cells were determined to be female while embryos with higher fractions of Erdr1, Ddx3y or Eif2s3y were determined to be male.
