#!/bin/bash
#SBATCH --mem=120g
#SBATCH -t 24:00:00
#SBATCH -n 16
#SBATCH -p general

module load anaconda
source activate Gilbert-work

gDir="Human/STAR"
RSEM="Human/RSEM/hg38"
gtf="Homo_sapiens.GRCh38.108.chr.gtf"
fastq="RIP/HEK"
outDir="RIP-seq/HEK"
file=$1

STAR \
        --runThreadN 32 \
        --sjdbGTFfile ${gtf} \
        --genomeDir ${gDir} \
        --readFilesIn ${fastq}/${file}.fastq.gz \
        --readFilesCommand zcat \
        --outFilterType BySJout \
        --genomeLoad NoSharedMemory \
        --quantMode TranscriptomeSAM \
        --outFileNamePrefix ${outDir}/STAR/${file} \
        --runRNGseed 777 \
        --outSAMtype BAM Unsorted

samtools sort -@ 32 ${outDir}/STAR/${file}Aligned.out.bam > ${outDir}/STAR/${file}Aligned.sorted.bam
samtools index -@ 32 ${outDir}/STAR/${file}Aligned.sorted.bam
rm -f ${outDir}/${file}Aligned.out.bam


#RSEM
rsem-calculate-expression --num-threads 32 \
        --alignments \
        --seed 77 \
        --no-bam-output \
        ${outDir}/STAR/${file}Aligned.toTranscriptome.out.bam \
        ${RSEM} \
        ${outDir}/RSEM/${file}
