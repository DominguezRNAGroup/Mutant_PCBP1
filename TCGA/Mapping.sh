gDir="Human/STAR"
gtf="Human/Homo_sapiens.GRCh38.108.chr.gtf"
fastq="Fastq"
outDir="STAR_Mapped"
file=$1

#single-end mapping
STAR --runThreadN 32 \
        --sjdbGTFfile ${gtf} \
        --outSAMstrandField intronMotif \
        --alignIntronMax 299999 \
        --alignEndsType EndToEnd \
        --genomeDir ${gDir} \
        --readFilesIn ${fastq}/${file}_1.fq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix ${outDir}/${file} \
        --runRNGseed 777 \
        --outSAMtype BAM Unsorted \
        --quantMode TranscriptomeSAM

#paired-end mapping
STAR --runThreadN 32 \
        --sjdbGTFfile ${gtf} \
        --alignEndsType EndToEnd \
        --outSAMstrandField intronMotif \
        --alignIntronMax 299999 \
        --genomeDir ${gDir} \
        --readFilesIn ${fastq}/${file}_1.fq.gz ${fastq}/${file}_2.fq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix ${outDir}/${file} \
        --runRNGseed 777 \
        --outSAMtype BAM Unsorted \
        --quantMode TranscriptomeSAM
