fastq="Fastq/"
outDir="Polysome"
file=$1
threads=32

module load bwa-mem2/2.2
module load fastqc/0.12.1
module load fastp/0.23.2


fastp -i ${outDir}/Fastq/${file}_CKDL250003877-1A_22YYLLLT3_L6_1.fq \
       -o ${outDir}/Fastq/${file}_CKDL250003877-1A_22YYLLLT3_L6_1_tr.fq \
       -I ${outDir}/Fastq/${file}_CKDL250003877-1A_22YYLLLT3_L6_2.fq \
       -O ${outDir}/Fastq/${file}_CKDL250003877-1A_22YYLLLT3_L6_2_tr.fq \
       --length_required 30 \
       --trim_poly_g

fastqc -t 2 --extract -k 7 ${outDir}/Fastq/${file}_CKDL250003877-1A_22YYLLLT3_L6_1_tr.fq -o ${outDir}/Fastq
fastqc -t 2 --extract -k 7 ${outDir}/Fastq/${file}_CKDL250003877-1A_22YYLLLT3_L6_2_tr.fq -o ${outDir}/Fastq

#mapping to SF9 to account for spike-ins
bwa-mem2 mem -t ${threads} GENOME/SF9/sf9 ${outDir}/Fastq/${file}_CKDL250003877-1A_22YYLLLT3_L6_1_tr.fq ${outDir}/Fastq/${file}_CKDL250003877-1A_22YYLLLT3_L6_2_tr.fq -o ${outDir}/Fly/${file}.sam
samtools view -@ 4 -bh ${outDir}/Fly/${file}.sam > ${outDir}/Fly/${file}.bam
rm -f ${outDir}/Fly/${file}.sam

samtools view -@ 4 -f 4 -bh ${outDir}/Fly/${file}.bam > ${outDir}/Fly/${file}_Unmapped.bam
samtools sort -@ 4 -n ${outDir}/Fly/${file}_Unmapped.bam > ${outDir}/Fly/${file}_Unmapped_sorted.bam
samtools flagstat ${outDir}/Fly/${file}.bam > ${outDir}/Fly/${file}_flagstat.txt

module load bedtools/2.30
bedtools bamtofastq -i ${outDir}/Fly/${file}_Unmapped_sorted.bam -fq ${outDir}/Fastq/${file}_unmapped_R1.fq -fq2 ${outDir}/Fastq/${file}_unmapped_R2.fq

fastq-sort --id ${outDir}/Fastq/${file}_unmapped_R1.fq > ${outDir}/Fastq/${file}_unmapped_R1_sorted.fq
fastq-sort --id ${outDir}/Fastq/${file}_unmapped_R2.fq > ${outDir}/Fastq/${file}_unmapped_R2_sorted.fq
