#!/bin/sh

for i in {0..728}
do
    cd job_$i
    mpiexec -n 4 --allow-run-as-root ~/bin/xstelloptv2 input.HSX_test
    rm boozmn* input.HSX_test.00000 jxbout* mercier* neolog.HSX_test.00000 parvmecinfo.txt stellopt.HSX_test threed1.HSX_test timings.txt var_labels wout* xvec.dat mgrid_hsx_complete.nc
    cd ../
done
