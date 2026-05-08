 #!/bin/bash
BedDir="CLIP/Bed"
BamDir="eCLIP/STAR_Hum"
#samples="samples_matched.txt"
scriptDir="Clipper-perl/merge_peaks/bin/perl"

source activate clipper-perl
export PERL5LIB=/nas/longleaf/home/ggiri/perl5/lib/perl5:$PERL5LIB


while FS=$"\t" read -r -a line;
do
        ip=${line[0]}
        input=${line[1]}
       	echo "${ip} ${input}"
        sbatch -J ${ip} -o ${ip}.o -e ${ip}.e --mem=120g -t 24:00:00 --wrap="perl ${scriptDir}/overlap_peakfi_with_bam.pl ${BamDir}/${ip}Aligned.sorted.dedup.bam ${BamDir}/${input}Aligned.sorted.dedup.bam ${BedDir}/${ip}.bed ${BamDir}/${ip}Aligned.sorted.dedup.txt ${BamDir}/${input}Aligned.sorted.dedup.txt ${BedDir}/${ip}.norm.bed"
done < ${samples}


while FS=$"\t" read -r -a line;
do
	ip=${line[0]}
	sbatch -J ${ip} -o ${ip}.o -e ${ip}.e --mem=120g -t 24:00:00 --wrap="perl ${scriptDir}/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.pl ${BedDir}/${ip}.norm.bed ${BedDir}/${ip}.norm.compressed.bed"
done < ${samples}


while FS=$"\t" read -r -a line;
do
       ip=${line[0]}
       sbatch -J ${ip} -o ${ip}.o -e ${ip}.e --mem=16g -t 01:00:00 --wrap="perl ${scriptDir}/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl ${BedDir}/${ip}.norm.bed.full ${BedDir}/${ip}.norm.compressed2.bed ${BedDir}/${ip}.norm.compressed.bed.full"
done < ${samples}



#make information content
while FS=$"\t" read -r -a line;
do
	ip=${line[0]}
	input=${line[1]}
	echo "perl ${scriptDir}/make_informationcontent_from_peaks.pl ${BedDir}/${BedDir}/${ip}.norm.compressed.bed.full ${BamDir}/${ip}Aligned.sorted.dedup.txt ${BamDir}/${input}Aligned.sorted.dedup.txt ${BedDir}/${ip}.norm.compressed.bed.entropy.full ${BedDir}/${ip}.norm.compressed.bed.entropy.excessreads"
	sbatch -J ${ip} -o ${ip}.o -e ${ip}.e --mem=16g -t 01:00:00 --wrap="perl ${scriptDir}/make_informationcontent_from_peaks.pl ${BedDir}/${ip}.norm.compressed.bed.full ${BamDir}/${ip}Aligned.sorted.dedup.txt ${BamDir}/${input}Aligned.sorted.dedup.txt ${BedDir}/${ip}.norm.compressed.bed.entropy.full ${BedDir}/${ip}.norm.compressed.bed.entropy.excessreads"
done < ${samples}

#################AFTER IDR
#parsing idr peaks
samples="samples_reps.txt"
while FS=$"\t" read -r -a line;
do
       rep1=${line[0]}
       rep2=${line[1]}
       sbatch -J ${rep1} -o ${rep1}.o -e ${rep1}.e --mem=16g -t 01:00:00 --wrap="perl ${scriptDir}/parse_idr_peaks.pl ${BedDir}/${rep1}_${rep2}_idr.bed ${BedDir}/${rep1}.norm.compressed.bed.entropy.full ${BedDir}/${rep2}.norm.compressed.bed.entropy.full ${BedDir}/${rep1}_${rep2}_idr_parsed.bed"
done < ${samples}


samples="samples_matched.txt"
while FS=$"\t" read -r -a line;
do
        ip=${line[0]}
        input=${line[1]}
	idr=${line[2]}
	sbatch -J ${ip} -o ${ip}.o -e ${ip}.e --mem=32g -t 01:00:00 --wrap="perl ${scriptDir}/overlap_peakfi_with_bam.pl ${BamDir}/${ip}Aligned.sorted.dedup.bam ${BamDir}/${input}Aligned.sorted.dedup.bam ${BedDir}/${idr}_idr_parsed.bed ${BamDir}/${ip}Aligned.sorted.dedup.txt ${BamDir}/${input}Aligned.sorted.dedup.txt ${BedDir}/${ip}.norm.idr.merged.bed"
done < ${samples}

samples="samples_matched.txt"
while FS=$"\t" read -r -a line;
do
	ip=${line[0]}
	input=${line[1]}
	idr=${line[2]}
	sbatch -J ${ip} -o ${ip}.o -e ${ip}.e --mem=32g -t 01:00:00 --wrap="perl ${scriptDir}/overlap_peakfi_wi$
done < ${samples}



samples="samples_reps.txt"
while FS=$"\t" read -r -a line;
do
       rep1=${line[0]}
       rep2=${line[1]}
       sbatch -J ${rep1} -o ${rep1}.o -e ${rep1}.e --mem=32g -t 01:00:00 --wrap="perl ${scriptDir}/get_reproducing_peaks.pl ${BedDir}/${rep1}.norm.idr.merged.bed.full ${BedDir}/${rep2}.norm.idr.merged.bed.full ${BedDir}/IDR/${rep1}_reprod.bed.full ${BedDir}/IDR/${rep2}_reprod.bed.full ${BedDir}/IDR/${rep1}_${rep2}_reprod.bed ${BedDir}/IDR/${rep1}_${rep2}.custombed ${BedDir}/${rep1}.norm.compressed.bed.entropy.full ${BedDir}/${rep2}.norm.compressed.bed.entropy.full ${BedDir}/${rep1}_${rep2}_idr.bed"
done < ${samples}
