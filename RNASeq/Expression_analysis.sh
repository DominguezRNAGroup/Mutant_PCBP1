gDir="GENOME/Human/STAR"
RSEM="GENOME/Human/RSEM/hg38"
gtf="GENOME/Human/Homo_sapiens.GRCh38.108.chr.gtf"
fastq="Fastq"
outDir="Exp"
file=$1

#Mapping
STAR \
        --runThreadN 32 \
        --alignEndsType EndToEnd \
        --outSAMstrandField intronMotif \
        --alignIntronMax 299999 \
        --sjdbGTFfile ${gtf} \
        --genomeDir ${gDir} \
        --readFilesIn ${fastq}/${file}_1.fastq.gz ${fastq}/${file}_2.fastq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix ${outDir}/STAR/${file} \
        --runRNGseed 777 \
        --outSAMtype BAM Unsorted \
        --quantMode TranscriptomeSAM

samtools sort -@ 32 ${outDir}/STAR/${file}Aligned.out.bam > ${outDir}/STAR/${file}Aligned.sorted.bam
samtools index -@ 32 ${outDir}/STAR/${file}Aligned.sorted.bam
rm -rf ${outDir}/STAR/${file}Aligned.out.bam

#Quantification
rsem-calculate-expression --paired-end \
        --num-threads 32 \
        --alignments \
        --seed 77 \
        --no-bam-output \
        ${outDir}/STAR/${file}Aligned.toTranscriptome.out.bam \
        ${RSEM} \
        ${outDir}/RSEM/${file}
