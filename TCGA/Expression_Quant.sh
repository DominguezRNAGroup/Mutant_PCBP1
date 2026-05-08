index="RSEM/hg38"
Bam="STAR_Mapped"
outDir="RSEM"
file=$1

#single end
rsem-calculate-expression --num-threads 32 \
        --alignments \
        --seed 77 \
        --no-bam-output \
        ${Bam}/${file}Aligned.toTranscriptome.out.bam \
        ${index} \
        ${outDir}/${file}

#paired end 
rsem-calculate-expression --paired-end \
        --num-threads 32 \
        --alignments \
        --seed 77 \
        --no-bam-output \
        ${Bam}/${file}Aligned.toTranscriptome.out.bam \
        ${index} \
        ${outDir}/${file}
