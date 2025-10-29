#!/bin/bash
#SBATCH -p barbun-cuda
#SBATCH -J qe_test
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --gres=gpu:1
#SBATCH --threads=1
#SBATCH --time=1-00:00:00
#SBATCH -o print_gpu.out    # Ciktinin yazilacagi dosya adi


echo "SLURM_NODELIST $SLURM_NODELIST"
echo "NUMBER OF CORES $SLURM_NTASKS"

module load lib/cuda/12.4

TIME_LIMIT=$((180 * 60))
#source /arf/sw/comp/oneapi/2023.0/setvars.sh
export PATH="/arf/home/aakdogan/opt/new/bin:$PATH"

check_exit_status() {
    local exit_code=$1
    local program_name=$2  # Name of the program being tested

    if [ $exit_code -eq 124 ]; then
        echo "EXIT $program_name: Timeout (Exceeded 1.5 hours) $exit_code."
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

i=11

search_dir=matrices
for entry in "$search_dir"/*
do
    rm result/result.mtx 2> /dev/null;
    rm result/order.mtx 2> /dev/null;
    echo "$entry"
    read -r MEM EXIT_CODE <<< "$(
    {
        __time -v bash -c "timeout $TIME_LIMIT make opt_sequential_cuda_truba i=$i && ./opt_sequential $entry >> result_gpu.log" 2>&1
        echo "__EXIT_CODE__:$?"
    } 2>&1 | awk -F: '
        /Maximum resident set size/ {gsub(/[^0-9]/,"",$2); mem=$2}
        /__EXIT_CODE__/ {code=$2}
        END {print mem, code}'
    )"
    echo "Peak Memory: ${MEM} KB" >> result_gpu.log
    echo "" >> result_gpu.log;
    RESULT=$(check_exit_status $EXIT_CODE "opt_sequential")
    echo "$RESULT" >> result_gpu.log;
    echo "-----------" >> result_gpu.log;
    echo "" >> result_gpu.log;
done
#make i=$i && ./opt_sequential 


exit

