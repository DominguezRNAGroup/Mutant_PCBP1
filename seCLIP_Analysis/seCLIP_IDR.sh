module load anaconda
source activate Clipper_env

samples="samples_reps.txt"
scriptDir="Clipper-perl/merge_peaks/bin/"
BedDir="eCLIP/Bed"

#make information content
while FS=$"\t" read -r -a line;
do
       ip=${line[0]}
	input=${line[1]}
	sbatch -J ${ip} -o ${ip}.o -e ${ip}.e --mem=16g -t 01:00:00 --wrap="python ${scriptDir}/full_to_bed.py --input ${BedDir}/${ip}.norm.compressed.bed.entropy.full --output ${BedDir}/${ip}.norm.compressed.bed.entropy.bed"
done < ${samples}



#Getting IDR peaks
while FS=$"\t" read -r -a line;
do
  	rep1=${line[0]}
        rep2=${line[1]}
        sbatch -J ${rep1} -o ${rep1}.o -e ${rep1}.e --mem=16g -t 01:00:00 --wrap="idr-2.0.2/bin/idr --samples ${BedDir}/${rep1}.norm.compressed.bed.entropy.bed ${BedDir}/${rep2}.norm.compressed.bed.entropy.bed --input-file-type bed --rank 5 --peak-merge-method max --plot -o ${BedDir}/${rep1}_${rep2}_idr.bed"
done < ${samples}
