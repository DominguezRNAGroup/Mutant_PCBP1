#!/bin/bash

# CHANGE THIS TO MATCH YOUR INSTALL PATH
enspara_cluster_path=/home/louis/miniforge3/pkgs/enspara-0.3.1-py_0/info/recipe/enspara/apps/cluster.py

# stride rate for downsampling.
# ds=10
# number of clusters to do.
# k=100
radius=0.30
# tag=rad-$radius
# number of medoid sweeps to do; note that each medoid pass 
# will take approx as long as the original kcenters.
iters=5  
# file containing mdtraj selection string for atoms to compute pairwise RMSD for
base_selstr=bb-NO.txt
# note that you could replace the 'atoms' flag line below with a straight-up string;
# this would make sense if you were going for something simple, like 'name CA'. So like
# '--atoms "name CA"'
distance=rmsd
export OMP_NUM_THREADS=4
export NUMEXPR_MAX_THREADS=4
echo OMP_NUM_THREADS $OMP_NUM_THREADS
clus_sess_dir=rmsd-clustering-$(date -I)
echo starting to fill $clus_sess_dir
# your cluster might have a different system for getting its MPI env set up. 
# Do whatever dance you need to to make it work here.
# module load openmpi/5.0.3
for protp in vas-sims/*; do 
    protn=$(basename $protp)
    if [[ $protn == UP1 ]]; then
        clustag=short-rad-$radius
        selstr_file=short-sel-up1.txt
    elif [[ $protn == PCBP1 ]]; then
        clustag=short-rad-$radius
        selstr_file=short-sel-pcbp1.txt
    else
        clustag=rad-$radius
        selstr_file=$base_selstr
    fi
    for varp in $protp/*; do
        varn=$(basename $varp)
        radius=0.30
        tag=$protn-$varn-$clustag
        # number of medoid sweeps to do; note that each medoid pass 
        # will take approx as long as the original kcenters.
        iters=5    # list of trajectories in a newline delimited text file
        traj_list=traj-paths-$protn-$varn.txt
        # the topology that matches te trajectory files
        top=$varp/filtered/filtered.pdb
        outdir=$clus_sess_dir/$tag
        if [[ ! -d outdir ]]; then
            mkdir -p $outdir
        fi
        if [[ -f stop-cluster ]]; then
            echo 'Found cluster stop file; exiting.'
        fi
        nts=$(wc -l $traj_list | awk '{print ($1 < 64)? $1 : 64}')
        echo working on $protn $varn
        # Run the clustering app with mpiexec
        mpiexec -np $nts python $enspara_cluster_path \
        --trajectories $(tr '\n' ' ' < $traj_list) \
        --topology $top \
        --cluster-distance $distance\
        --algorithm khybrid \
        --cluster-iterations $iters \
        --cluster-radius $radius \
        --atoms "$(cat $selstr_file)" \
        --distances $outdir/distances.h5 \
        --center-features $outdir/centers.pickle \
        --assignments $outdir/assigns.h5 2>&1 > rmsd-clustering/$tag.out
    done
done
