#!/bin/bash

base=$1
mkdir -p ${base}_surrogates

declare -i increment=15000
declare -i sur_pair_count=$increment

while test -z "$done"; do

    echo "generating surrogates for $sur_pair_count surrogate CCH's"

    if ! mpirun --mca pml ob1 --map-by core -hostfile hostfile testcch basename $base sur_pair_count $sur_pair_count surrogates; then
        echo "surrogate generation failed"
        exit 1
    fi

    echo "generating $sur_pair_count surrogate CCH's"

    if ! mpirun --mca pml ob1 --map-by core -hostfile hostfile testcch basename $base sur_pair_count $sur_pair_count ; then
        echo "surrogate CCH generation failed"
        exit 1
    fi
    
    threshold=`mpirun --mca pml ob1 -np 1 testcch basename $base sur_pair_count $sur_pair_count threshold`
    
    if [ -n "$threshold" -a "$threshold" -le $sur_pair_count ] ; then
        done=1
    else
        sur_pair_count+=$increment
    fi
    
done

echo $threshold > ${base}_surrogates/spcnt
