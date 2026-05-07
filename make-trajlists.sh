#!/bin/bash

for protp in vas-sims/*; do
    protn=$(basename $protp)
    for varp in $protp/*; do
        varn=$(basename $varp)
        echo making $protn $varn
        fdfind -I output.filtered.xtc $varp | sort -V > traj-paths-$protn-$varn.txt
        cat traj-paths-$protn-$varn.txt
    done
done