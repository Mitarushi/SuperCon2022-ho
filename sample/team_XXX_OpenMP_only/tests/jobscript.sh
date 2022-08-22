#!/bin/bash
if [ ! "${IS_MPI_PROGRAM}" == "0" ] && [ ! "${IS_MPI_PROGRAM}" == "1" ]; then
    echo "Set \"IS_MPI_PROGRAM=0\" or \"IS_MPI_PROGRAM=1\""
    exit 1
fi

date
echo "PJM_NODE:"
echo ${PJM_NODE}
echo "PJM_MPI_PROC:"
echo ${PJM_MPI_PROC}
echo "OMP_NUM_THREADS:"
echo ${OMP_NUM_THREADS}
echo "IS_MPI_PROGRAM:"
echo ${IS_MPI_PROGRAM}

CASES="random_01.in random_02.in random_03.in random_04.in random_05.in random_06.in random_07.in random_08.in random_09.in random_10.in"
# CASES="00_sample_01.in 00_sample_02.in random_01.in"

# compile the programs
cd tests
    FCC -Ofast -fopenmp -Nclang -std=c++17 generator.cpp -o generator
    ./generator
    FCC -Ofast -fopenmp -Nclang -std=c++17 judge.cpp -o judge
cd ..
if [ ${IS_MPI_PROGRAM} -eq 1 ]; then 
    mpiFCC -Ofast -fopenmp -Nclang -std=c++17 main.cpp -o log/main.${PJM_JOBID}
else
    FCC -Ofast -fopenmp -Nclang -std=c++17 main.cpp -o log/main.${PJM_JOBID}
fi

for CASE in $CASES; do
    OUT="log/${CASE}.${PJM_JOBID}.out"
    ERR="log/${CASE}.${PJM_JOBID}.err"
    JUDGE="log/${CASE}.${PJM_JOBID}.judge"
# execute the program
    if [ ${IS_MPI_PROGRAM} -eq 1 ]; then 
        timeout 70.0s mpiexec --stdin tests/$CASE --stdout $OUT --stderr $ERR ./log/main.${PJM_JOBID}
    else
        timeout 70.0s ./log/main.${PJM_JOBID} < tests/${CASE} > $OUT 2> $ERR
    fi
# judge
    ./tests/judge $OUT < tests/${CASE} > $JUDGE
done

# output the result
RESULT="log/result.${PJM_JOBID}.out"
echo "TEAM:" >> $RESULT
echo $(basename $(pwd)) >> $RESULT
echo "PJM_NODE:" >> $RESULT
echo ${PJM_NODE} >> $RESULT
echo "PJM_MPI_PROC:" >> $RESULT
echo ${PJM_MPI_PROC} >> $RESULT
echo "OMP_NUM_THREADS:" >> $RESULT
echo ${OMP_NUM_THREADS} >> $RESULT
echo "IS_MPI_PROGRAM:"  >> $RESULT
echo ${IS_MPI_PROGRAM} >> $RESULT

for CASE in $CASES; do
    JUDGE="log/${CASE}.${PJM_JOBID}.judge"
    # output the score and the runnning time (seconds)
    echo "${CASE}_SCORE:" >> $RESULT
    echo $(sed -n 1p $JUDGE) >> $RESULT
    echo "${CASE}_TIME:" >> $RESULT
    echo $(sed -n 2p $JUDGE) >> $RESULT
done
