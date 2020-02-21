# SNPsplit
wget https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/SNPsplit_v0.3.2.tar.gz
tar -zxvf SNPsplit_v0.3.2.tar.gz

# build N-masked genome of B6 (GRCm38) and Cast based on SNP data and corresponding STAR index
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
wget ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa

perl SNPsplit_genome_preparation --vcf_file mgp.v5.merged.snps_all.dbSNP142.vcf.gz --strain CAST_EiJ --reference_genome 

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir  --genomeFastaFiles CAST_EiJ_N-masked/chr1.N-masked.fa CAST_EiJ_N-masked/chr10.N-masked.fa CAST_EiJ_N-masked/chr11.N-masked.fa CAST_EiJ_N-masked/chr12.N-masked.fa CAST_EiJ_N-masked/chr13.N-masked.fa CAST_EiJ_N-masked/chr14.N-masked.fa CAST_EiJ_N-masked/chr15.N-masked.fa CAST_EiJ_N-masked/chr16.N-masked.fa CAST_EiJ_N-masked/chr17.N-masked.fa CAST_EiJ_N-masked/chr18.N-masked.fa CAST_EiJ_N-masked/chr19.N-masked.fa CAST_EiJ_N-masked/chr2.N-masked.fa CAST_EiJ_N-masked/chr3.N-masked.fa CAST_EiJ_N-masked/chr4.N-masked.fa CAST_EiJ_N-masked/chr5.N-masked.fa CAST_EiJ_N-masked/chr6.N-masked.fa CAST_EiJ_N-masked/chr7.N-masked.fa CAST_EiJ_N-masked/chr8.N-masked.fa CAST_EiJ_N-masked/chr9.N-masked.fa CAST_EiJ_N-masked/chrX.N-masked.fa CAST_EiJ_N-masked/chrY.N-masked.fa

# convert SNP annotation file to bed
zcat all_SNPs_CAST_EiJ_GRCm38.txt.gz | perl -ane '$s=$F[2]-1; print "chr$F[1]\t$s\t$F[2]\t$F[0];$F[3];$F[4]\n"' >all_SNPs_CAST_EiJ_GRCm38.bed

