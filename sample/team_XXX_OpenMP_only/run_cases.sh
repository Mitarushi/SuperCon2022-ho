#!/bin/bash
#PJM --rsc-list "elapse=0:20:00"
#PJM --out "log/jobscript.sh.%j.out"
#PJM -j
#PJM -L node=1

export IS_MPI_PROGRAM=0
export OMP_NUM_THREADS=48
export PLE_MPI_STD_EMPTYFILE=off
./tests/jobscript.sh

# MPI を使用しない場合のジョブスクリプトの例。
# #PJM -L node=1
# 
# export IS_MPI_PROGRAM=0
# export OMP_NUM_THREADS=Y
# の Y をある整数値に指定する（上は Y = 48 の場合）。それ以外の変更は加えないこと。
# ・node の値は、ノード数を表す。
# 
# ・IS_MPI_PROGRAM の値は、MPIを使用する場合 1、使用しない場合 0 を指定する。
# ・OMP_NUM_THREADS の値 (Y) は、OpenMPのスレッド数を表す。
