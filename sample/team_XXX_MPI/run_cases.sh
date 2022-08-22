#!/bin/bash
#PJM --rsc-list "elapse=0:20:00"
#PJM --out "log/jobscript.sh.%j.out"
#PJM -j
#PJM -L node=1
#PJM --mpi proc=4
export IS_MPI_PROGRAM=1
export OMP_NUM_THREADS=12
export PLE_MPI_STD_EMPTYFILE=off
./tests/jobscript.sh

# MPI を使用する（MPIとOpenMPを使用する場合も含む）の場合のジョブスクリプトの例。
# #PJM -L node=1
# #PJM --mpi proc=X
# export IS_MPI_PROGRAM=1
# export OMP_NUM_THREADS=Y
# の X, Y をある整数値に指定する（上は X = 4, Y = 12 の場合）。それ以外の変更は加えないこと。
# ・node の値は、ノード数を表す。
# ・proc の値 (X) は、MPIのプロセス数。
# ・IS_MPI_PROGRAM の値は、MPIを使用する場合 1、使用しない場合 0 を指定する。
# ・OMP_NUM_THREADS の値 (Y) は、OpenMPのスレッド数を表す。
