#!/bin/bash
#SBATCH --mem=60g
#SBATCH -t 48:00:00
#SBATCH -n 16
#SBATCH -p general

#FastQC v0.11.9
#Python 3.9.7
#UMI-tools version: 1.1.2
#cutadapt-4.1
#fastqc-tools v0.8.3

module load anaconda
source activate Clipper_env

Fastq="/work/users/g/g/ggiri/PCBP1_Final/eCLIP/Fastq"
OutDir="/work/users/g/g/ggiri/PCBP1_Final/eCLIP"
RepDir="/proj/RNA_lab/Gilbert/GENOME/Repeats_New/STAR"
gDir="/work/users/g/g/ggiri/GENOME/Human/GENECODE/STAR"
file=$1


#Extracting UMI
umi_tools extract --random-seed 777 \
        --bc-pattern NNNNNNNNNN \
        -L ${Fastq}/${file}.metrics \
        --stdin ${Fastq}/${file}_CKDL240027093-1A_22N5C7LT3_L3_1.fq.gz \
        --stdout ${OutDir}/Trimmed/${file}_1_umi.fq.gz
#
###trimming
cutadapt -O 1 \
       --match-read-wildcards \
       --times 1 \
       -e 0.1 \
       --quality-cutoff 6 \
       -m 18 \
       -o ${OutDir}/Trimmed/${file}_1_umi_tr.fq.gz \
       -a AGATCGGAAGAGCAC \
       -a GATCGGAAGAGCACA \
       -a ATCGGAAGAGCACAC \
       -a TCGGAAGAGCACACG \
       -a CGGAAGAGCACACGT \
       -a GGAAGAGCACACGTC \
       -a GAAGAGCACACGTCT \
       -a AAGAGCACACGTCTG \
       -a AGAGCACACGTCTGA \
       -a GAGCACACGTCTGAA \
       -a AGCACACGTCTGAAC \
       -a GCACACGTCTGAACT \
       -a CACACGTCTGAACTC \
       -a ACACGTCTGAACTCC \
       -a CACGTCTGAACTCCA \
       -a ACGTCTGAACTCCAG \
       -a CGTCTGAACTCCAGT \
       -a GTCTGAACTCCAGTC \
       -a TCTGAACTCCAGTCA \
       -a CTGAACTCCAGTCAC \
       ${OutDir}/Trimmed/${file}_1_umi.fq.gz > ${OutDir}/Trimmed/${file}_1_metrics.txt
#
##trimmingx2
cutadapt -O 5 \
       --match-read-wildcards \
       --times 1 \
       -e 0.1 \
       --quality-cutoff 6 \
       -m 18 \
       -o ${OutDir}/Trimmed/${file}_1_umi_tr2.fq.gz \
       -a AGATCGGAAGAGCAC \
       -a GATCGGAAGAGCACA \
       -a ATCGGAAGAGCACAC \
       -a TCGGAAGAGCACACG \
       -a CGGAAGAGCACACGT \
       -a GGAAGAGCACACGTC \
       -a GAAGAGCACACGTCT \
       -a AAGAGCACACGTCTG \
       -a AGAGCACACGTCTGA \
       -a GAGCACACGTCTGAA \
       -a AGCACACGTCTGAAC \
       -a GCACACGTCTGAACT \
       -a CACACGTCTGAACTC \
       -a ACACGTCTGAACTCC \
       -a CACGTCTGAACTCCA \
       -a ACGTCTGAACTCCAG \
       -a CGTCTGAACTCCAGT \
       -a GTCTGAACTCCAGTC \
       -a TCTGAACTCCAGTCA \
       -a CTGAACTCCAGTCAC \
       ${OutDir}/Trimmed/${file}_1_umi_tr.fq.gz > ${OutDir}/Trimmed/${file}_2_metrics.txt


#FastQC
gunzip ${OutDir}/Trimmed/${file}_1_umi_tr2.fq.gz
fastqc -t 2 --extract -k 7 ${OutDir}/Trimmed/${file}_1_umi_tr2.fq -o ${OutDir}/Trimmed
fastqc -t 2 --extract -k 7 ${Fastq}/${file}_CKDL240027093-1A_22N5C7LT3_L3_1.fq.gz -o ${OutDir}/Trimmed

#Sorting the Fastq
/work/users/g/g/ggiri/Tools/fastq-tools/src/fastq-sort --id ${OutDir}/Trimmed/${file}_1_umi_tr2.fq > ${OutDir}/Trimmed/${file}_1_umi_tr2_sorted.fq
#
##mapping to repeat elements
STAR \
       --alignEndsType EndToEnd \
       --genomeDir ${RepDir} \
       --genomeLoad NoSharedMemory \
       --outBAMcompression 10 \
       --outFileNamePrefix ${OutDir}/STAR_Rep/${file} \
       --outFilterMultimapNmax 30 \
       --outFilterMultimapScoreRange 1 \
       --outFilterScoreMin 10 \
       --outFilterType BySJout \
       --outReadsUnmapped Fastx \
       --outSAMattrRGline ID:foo \
       --outSAMattributes All \
       --outSAMmode Full \
       --outSAMtype BAM Unsorted \
       --outSAMunmapped Within \
       --outStd Log \
       --readFilesIn ${OutDir}/Trimmed/${file}_1_umi_tr2_sorted.fq \
       --runMode alignReads \
       --runThreadN 32

#
##Mapping to Human Ref Genome
STAR \
       --alignEndsType EndToEnd \
       --genomeDir ${gDir} \
       --genomeLoad NoSharedMemory \
       --outBAMcompression 10 \
       --outFileNamePrefix ${OutDir}/STAR_Hum/${file} \
       --outFilterMultimapNmax 1 \
       --outFilterMultimapScoreRange 1 \
       --outFilterScoreMin 10 \
       --outFilterType BySJout \
       --outReadsUnmapped Fastx \
       --outSAMattrRGline ID:foo \
       --outSAMattributes All \
       --outSAMmode Full \
       --outSAMtype BAM Unsorted \
       --outSAMunmapped Within \
       --outStd Log \
       --readFilesIn ${OutDir}/STAR_Rep/${file}Unmapped.out.mate1 \
       --runMode alignReads \
       --runThreadN 32


samtools sort -@ 32 ${OutDir}/STAR_Hum/${file}Aligned.out.bam > ${OutDir}/STAR_Hum/${file}Aligned.sorted.bam
samtools index -@ 32 ${OutDir}/STAR_Hum/${file}Aligned.sorted.bam

#Dedup
umi_tools dedup --random-seed 777 \
       -I ${OutDir}/STAR_Hum/${file}Aligned.sorted.bam \
       --method unique \
       --output-stats ${OutDir}/STAR_Hum/${file}_Dedup_Stats \
       -S ${OutDir}/STAR_Hum/${file}Aligned.sorted.dedup.bam

samtools index -@ 32 ${OutDir}/STAR_Hum/${file}Aligned.sorted.dedup.bam

#module load clipper/2.1.2
clipper -b ${OutDir}/STAR_Hum/${file}Aligned.sorted.dedup.bam \
       -o ${OutDir}/Bed/${file}.bed \
       -s GRCh38

samtools view -cF 4 ${OutDir}/STAR_Hum/${file}Aligned.sorted.dedup.bam > ${OutDir}/STAR_Hum/${file}Aligned.sorted.dedup.txt
