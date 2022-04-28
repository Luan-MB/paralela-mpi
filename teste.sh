#!/bin/bash

for threads in 1 2
do
    echo ${threads}
    for size in 10000 20000 40000 80000 160000
    do

        python3 input_generator.py ${threads} > input.in

        echo $size
        for i in {1..20..1}
        do
            echo $i
            mpirun -np ${threads} ./parallel ${size}.in >> dados/test_${size}_${threads}procs.txt
        done
    done
done
