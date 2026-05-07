module load anaconda
source activate Gilbert-work

gDir="Human/STAR"
RSEM="Human/RSEM/hg38"
gtf="Human/Homo_sapiens.GRCh38.108.chr.gtf"
outDir="Polysome"
file=$1


STAR \
        --runThreadN 32 \
        --outSAMstrandField intronMotif \
        --alignIntronMax 299999 \
        --sjdbGTFfile ${gtf} \
        --genomeDir ${gDir} \
        --readFilesIn ${outDir}/Fastq/${file}_unmapped_R1_sorted.fq ${outDir}/Fastq/${file}_unmapped_R2_sorted.fq \
        --outFileNamePrefix ${outDir}/STAR/${file} \
       --runRNGseed 777 \
        --outSAMtype BAM Unsorted \
        --quantMode TranscriptomeSAM

samtools sort -@ 32 ${outDir}/STAR/${file}Aligned.out.bam > ${outDir}/STAR/${file}Aligned.sorted.bam
samtools index -@ 32 ${outDir}/STAR/${file}Aligned.sorted.bam


#RSEM
rsem-calculate-expression --paired-end \
        --num-threads 32 \
        --alignments \
        --seed 77 \
        --no-bam-output \
        ${outDir}/STAR/${file}.Aligned.toTranscriptome.out.bam \
        ${RSEM} \
        ${outDir}/RSEM/${file}
