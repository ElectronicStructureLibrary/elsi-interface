#! /bin/bash

outdir=$3
nproc=$2

curdir=$PWD

cd $1
##ls comm_stat* | grep -v tgz | $curdir/parse_avg 0 $nproc 1  > $outdir/bcastU_total.dat 2>/dev/null
##ls comm_stat* | grep -v tgz | $curdir/parse_avg 0 $nproc 3  > $outdir/bcastL_total.dat 2>/dev/null
##ls comm_stat* | grep -v tgz | $curdir/parse_avg 0 $nproc 4  > $outdir/reduceL_total.dat 2>/dev/null
#ls comm_stat* | grep -v tgz | $curdir/parse_avg 0 $nproc  > $outdir/total.dat 2>/dev/null
##ls comm_stat* | grep -v tgz | $curdir/parse_avg 1 $nproc  > $outdir/avg.dat 2>/dev/null
##ls comm_stat* | grep -v tgz | ../parse_avg 0 $nproc 9  > send_L_CD_total.dat 2>/dev/null
#ls comm_stat* | grep -v tgz | $curdir/parse_avg_send 0 $nproc  > $outdir/sender_total.dat 2>/dev/null
#ls comm_stat* | grep -v tgz | $curdir/parse_avg_send 0 $nproc 1  > $outdir/sender_bcastU_total.dat 2>/dev/null
#ls comm_stat* | grep -v tgz | $curdir/parse_avg_send 0 $nproc 3  > $outdir/sender_bcastL_total.dat 2>/dev/null
#ls comm_stat* | grep -v tgz | $curdir/parse_avg_send 0 $nproc 4  > $outdir/sender_reduceL_total.dat 2>/dev/null
ls comm_stat* | grep -v tgz | $curdir/parse_avg_send 0 $nproc 7  > $outdir/sender_reduceD_total.dat 2>/dev/null
#ls comm_stat* | grep -v tgz | ../parse_avg_send_msg 0 $nproc 1  > $outdir/sender_bcastU_total_msg.dat 2>/dev/null
cd ..
