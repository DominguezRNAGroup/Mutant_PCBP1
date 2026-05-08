#RKO
outdir="RKO/rMATS"
echo "run_rmats --b1 RKO_parent.txt --b2 RKO_HDR.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/Parent_HDR --tmp ${outdir}/Parent_HDR_tmp --individual-counts"
sbatch -J Parent_HDR -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e Parent_HDR.e -o Parent_HDR.o --wrap="run_rmats --b1 RKO_parent.txt --b2 RKO_HDR.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/Parent_HDR --tmp ${outdir}/Parent_HDR_tmp --individual-counts"
sbatch -J Parent_KO -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e Parent_KO.e -o Parent_KO.o --wrap="run_rmats --b1 RKO_parent.txt --b2 RKO_KO.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/Parent_KO --tmp ${outdir}/Parent_KO_tmp --individual-counts"
sbatch -J Parent_L100Q -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e Parent_L100Q.e -o Parent_L100Q.o --wrap="run_rmats --b1 RKO_parent.txt --b2 RKO_L100Q.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/Parent_L100Q --tmp ${outdir}/Parent_L100Q_tmp --individual-counts"

#HCT-New
outdir="HCT_new/rMATS"
sbatch -J Parent_HDR2 -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e Parent_HDR2.e -o Parent_HDR2.o --wrap="run_rmats --b1 HCT_new_parent.txt --b2 HCT_new_HDR2.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/Parent_HDR2 --tmp ${outdir}/Parent_HDR2_tmp --individual-counts"
sbatch -J Parent_HDR -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e Parent_HDR.e -o Parent_HDR.o --wrap="run_rmats --b1 HCT_new_parent.txt --b2 HCT_new_HDR.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/Parent_HDR --tmp ${outdir}/Parent_HDR_tmp --individual-counts"
sbatch -J Parent_KO -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e Parent_KO.e -o Parent_KO.o --wrap="run_rmats --b1 HCT_new_parent.txt --b2 HCT_new_KO.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/Parent_KO --tmp ${outdir}/Parent_KO_tmp --individual-counts"
sbatch -J Parent_F10 -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e Parent_F10.e -o Parent_F10.o --wrap="run_rmats --b1 HCT_new_parent.txt --b2 HCT_new_L100Q_F10.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/Parent_F10 --tmp ${outdir}/Parent_F10_tmp --individual-counts"
sbatch -J Parent_L100Q -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e Parent_L100Q.e -o Parent_L100Q.o --wrap="run_rmats --b1 HCT_new_parent.txt --b2 HCT_new_L100Q.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/Parent_L100Q --tmp ${outdir}/Parent_L100Q_tmp --individual-counts"

#HCT-old
outdir="HCT_old/rMATS"
sbatch -J WT_L100 -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e WT_L100.e -o WT_L100.o --wrap="run_rmats --b1 HCT_old_WT.txt --b2 HCT_old_L100Q.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/WT_L100 --tmp ${outdir}/WT_L100_tmp --individual-counts"
sbatch -J WT_L102 -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e WT_L102.e -o WT_L102.o --wrap="run_rmats --b1 HCT_old_WT.txt --b2 HCT_old_L102Q.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/WT_L102 --tmp ${outdir}/WT_L102_tmp --individual-counts"

#Porter
outdir="Porter/rMATS"
sbatch -J WT_KD -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e WT_KD.e -o WT_KD.o --wrap="run_rmats --b1 Porter_WT.txt --b2 Porter_KD.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/WT_KD --tmp ${outdir}/WT_KD_tmp --individual-counts"
sbatch -J WT_KD_WT -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e WT_KD_WT.e -o WT_KD_WT.o --wrap="run_rmats --b1 Porter_WT.txt --b2 Porter_KD_WT.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/WT_KD_WT --tmp ${outdir}/WT_KD_WT_tmp --individual-counts"
sbatch -J WT_KD_L100 -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e WT_KD_L100.e -o WT_KD_L100.o --wrap="run_rmats --b1 Porter_WT.txt --b2 Porter_KD_L100Q.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/WT_KD_L100 --tmp ${outdir}/WT_KD_L100_tmp --individual-counts"

#C2BB
outdir="C2BBe1/rMATS"
sbatch -J GFP_WT -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e GFP_WT.e -o GFP_WT.o --wrap="run_rmats --b1 C2BBe1_GFP.txt --b2 C2BBe1_WT.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/GFP_WT --tmp ${outdir}/GFP_WT_tmp --individual-counts"
sbatch -J GFP_L100Q -p general -t 4:00:00 --mem=120g -N 1 -n 16 -e GFP_WT.e -o GFP_L100Q.o --wrap="run_rmats --b1 C2BBe1_GFP.txt --b2 C2BBe1_L100Q.txt --gtf ${gtf} -t paired --readLength 150 --nthread 16 --od ${outdir}/GFP_L100Q --tmp ${outdir}/GFP_L100Q_tmp --individual-counts"
