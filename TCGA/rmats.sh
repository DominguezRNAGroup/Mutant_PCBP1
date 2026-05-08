run_rmats --b1 single_samples_mapped_ol.txt --gtf Homo_sapiens.GRCh38.108.chr.gtf -t single --readLength 76 --nthread 32 --od rMATS/Single --tmp rMATS/Single_tmp --statoff --individual-counts
run_rmats --b1 paired_samples_mapped_ol.txt --gtf Homo_sapiens.GRCh38.108.chr.gtf -t paired --readLength 48 --nthread 32 --od rMATS_2/Paired --tmp rMATS_2/Paired_tmp --statoff --individual-counts
