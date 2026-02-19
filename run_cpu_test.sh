#!/bin/bash
#SBATCH -p barbun
#SBATCH -J qe_test
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --threads=1
#SBATCH --time=3-00:00:00


echo "SLURM_NODELIST $SLURM_NODELIST"
echo "NUMBER OF CORES $SLURM_NTASKS"

export PATH="/arf/home/aakdogan/opt/new/bin:$PATH"


TIME_LIMIT=$((180 * 60))
#source /arf/sw/comp/oneapi/2023.0/setvars.sh

check_exit_status() {
    local exit_code=$1
    local program_name=$2  # Name of the program being tested

    if [ $exit_code -eq 124 ]; then
        echo "EXIT $program_name: Timeout (Exceeded 3 hours) $exit_code."
    elif [ $exit_code -eq 139 ]; then
        echo "EXIT $program_name: Segmentation Fault (SIGSEGV) $exit_code."
    elif [ $exit_code -eq 137 ]; then
        echo "EXIT $program_name: Out of Memory (30GB) (Killed by OS) $exit_code."
    elif [ $exit_code -ne 0 ]; then
        echo "EXIT $program_name: Exited with error code $exit_code."
    else
        echo "EXIT $program_name: Executed successfully $exit_code."
    fi
}

num_to_bool() {
    local num=$1

    if [ $num -eq 0 ]; then
        echo "true"
    else
        echo "false"
    fi
}


i=7
cd /arf/home/aakdogan/CMPE-Cholesky-Thesis
search_dir=matrices
for entry in "$search_dir"/*
do
    rm result/result.mtx 2> /dev/null;
    rm result/order.mtx 2> /dev/null;
    echo "$entry";
    read -r MEM EXIT_CODE <<< "$(
    {
        __time -v bash -c "timeout $TIME_LIMIT make i=$i && ./opt_sequential $entry >> result_my.log" 2>&1
        echo "__EXIT_CODE__:$?"
    } 2>&1 | awk -F: '
        /Maximum resident set size/ {gsub(/[^0-9]/,"",$2); mem=$2}
        /__EXIT_CODE__/ {code=$2}
        END {print mem, code}'
    )"
    echo "Peak Memory: ${MEM} KB" >> result_my.log
    echo "" >> result_my.log;
    RESULT=$(check_exit_status $EXIT_CODE "opt_sequential")
    echo "$RESULT" >> result_my.log;
    echo "-----------" >> result_my.log;
    echo "" >> result_my.log;
    
    export LD_LIBRARY_PATH="/arf/home/aakdogan/opt/new/lib:/arf/home/aakdogan/opt/new/lib64"
    read -r MEM EXIT_CODE <<< "$(
    {
        __time -v timeout "$TIME_LIMIT" /arf/home/aakdogan/opt/new/bin/octave --eval "test_case('$entry', $i, false, $(num_to_bool $EXIT_CODE))" >> result_oct10.log
        echo "__EXIT_CODE__:$?"
    } 2>&1 | awk -F: '
        /Maximum resident set size/ {gsub(/[^0-9]/,"",$2); mem=$2}
        /__EXIT_CODE__/ {code=$2}
        END {print mem, code}'
    )"
    echo "Peak Memory: ${MEM} KB" >> result_oct10.log
    echo "" >> result_oct10.log;
    RESULT=$(check_exit_status $EXIT_CODE "octave")
    echo "$RESULT" >> result_oct10.log;
    echo "-----------" >> result_oct10.log;
    echo "" >> result_oct10.log;
done
#make i=$i && ./opt_sequential 

exit

