#!/bin/bash

SAMPLE=${1}      # /path/prefix
BC=${2}             # /outs/raw_gene_bc_matrices/mm10/barcodes.tsv(.gz)
R1=${3}             # "..._R1_001.fastq.gz ..._R1_001.fastq.gz ..._R1_001.fastq.gz"

R2=`echo ${R1} | sed 's/R1/R2/g'`

# get BCs from 10X pipeline
if file --mime-type ${BC} | grep -q gzip$; then
  zcat ${BC} | cut -f1 -d'-' >${SAMPLE}_10X_barcodes.tsv
else
  cut -f1 -d'-' ${BC} >${SAMPLE}_10X_barcodes.tsv
fi

# assignment of BC 2 read, reduce to 10X BCs
zcat ${R1} | sed 's/ .*//' | paste - - - - | cut -f1,2 | sed 's/^@//'  | perl -ane '$F[1]=substr($F[1], 0, 16); print "$F[0]\t$F[1]\n"' >${SAMPLE}_read.BC.tmp
fgrep -f ${SAMPLE}_10X_barcodes.tsv ${SAMPLE}_read.BC.tmp | sort >${SAMPLE}_read.BC.tsv
rm ${SAMPLE}_read.BC.tmp

# join fastq files
cat ${R2} >${SAMPLE}_R2.fastq.gz

# alignment using STAR (parameters suggested by SNPsplit manual)
STAR --runThreadN 20 --genomeDir SNPsplit_v0.3.2_hybrid_genome --readFilesIn ${SAMPLE}_R2.fastq.gz --readFilesCommand zcat --alignEndsType EndToEnd --outSAMattributes NH HI NM MD --outSAMtype BAM Unsorted --outFileNamePrefix ${SAMPLE}_

# assignment of reads to genome1 (e.g. B6) or genome2 (e.g. CAST) using SNPsplit
perl SNPsplit --snp_file all_SNPs_CAST_EiJ_GRCm38.txt.gz ${SAMPLE}_Aligned.out.bam

# clean up
rm ${SAMPLE}_Aligned.out.SNPsplit_sort.txt ${SAMPLE}_Aligned.out.SNPsplit_report.txt ${SAMPLE}_Log.final.out ${SAMPLE}_Log.out ${SAMPLE}_Log.progress.out ${SAMPLE}_SJ.out.tab ${SAMPLE}_R2.fastq.gz ${SAMPLE}_Aligned.out.unassigned.bam ${SAMPLE}_Aligned.out.bam
rm -r ${SAMPLE}__STARtmp


# filter for unambiguous alignments (no different chromosomes, no different reference assignments) and conversion to bed format
samtools merge -n - ${SAMPLE}_Aligned.out.genome1.bam ${SAMPLE}_Aligned.out.genome2.bam | samtools view -h - | sed 's/XX:Z:./XX:i:/' | perl -ane 'if($_=~m/NH:i:/){if($_=~m/NH:i:1\t/){print $_}}else{print $_}' | samtools view -bS - | bamToBed -split -tag XX -i - | perl -ane 'print "chr$F[0]\t$F[2]\t$F[3]\t$F[4]\tG$F[5]\t$F[6]\n"' | sort -k4,4 >${SAMPLE}_splitread.bed

less ${SAMPLE}_splitread.bed | bedtools groupby -g 4 -c 1,2,3,5 -o distinct,collapse,collapse,distinct | perl -ane 'if($F[1]!~m/,/ && $F[4]!~m/,/){@s=split(/,/,$F[2]); @e=split(/,/,$F[3]); for($i=0; $i<scalar(@s); $i++){print "$F[1]\t$s[$i]\t$e[$i]\t$F[0]\t$F[4]\n"}}' >${SAMPLE}_unambiguous.tmp
bedtools sort -i ${SAMPLE}_unambiguous.tmp >${SAMPLE}_unambiguous.bed


# assignemnt of reads to SNP (reduce to SNPs variants confirmed for both genomes)
bedtools intersect -sorted -wa -wb -a ${SAMPLE}_unambiguous.bed -b all_SNPs_CAST_EiJ_GRCm38.bed | perl -ane 'print "$F[8]\t$F[5]\t$F[3]\t$F[4]\n"' >${SAMPLE}_SNP.tsv


# get only white list SNPs
fgrep -w -f WT_SNP_white.list.tsv ${SAMPLE}_SNP.tsv | perl -ane 'print "$F[1].$F[3]\t$F[2]\t$F[0]\n"' | sort -k2,2 | uniq >${SAMPLE}_whitelistSNP.tsv
cat WT_SNP_white.list.tsv ${SAMPLE}_SNP.tsv | perl -ane 'print "$F[1].$F[3]\t$F[2]\t$F[0]\n"' | sort -k2,2 | uniq >${SAMPLE}_allSNP.tsv

# convert read to BC information
join -1 2 -2 1 ${SAMPLE}_whitelistSNP.tsv ${SAMPLE}_read.BC.tsv | perl -ane 'print "$F[3]\t$F[1]\t$F[2]\n"' | sort | bedtools groupby -g 1,2 -c 3 -o  count_distinct >${SAMPLE}_SNPcount.tsv
join -1 2 -2 1 ${SAMPLE}_allSNP.tsv ${SAMPLE}_read.BC.tsv | perl -ane 'print "$F[3]\t$F[1]\t$F[2]\n"' | sort | bedtools groupby -g 1,2 -c 3 -o  count_distinct >${SAMPLE}_allSNPcount.tsv

# clean up
rm ${SAMPLE}_SNP.tsv ${SAMPLE}_unambiguous.bed ${SAMPLE}_unambiguous.tmp ${SAMPLE}_Aligned.out.genome1.bam ${SAMPLE}_Aligned.out.genome2.bam
