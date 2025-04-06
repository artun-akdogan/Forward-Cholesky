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


TIME_LIMIT=$((90 * 60))
MEMORY_LIMIT=$((10 * 1024 * 1024))

ulimit -v $MEMORY_LIMIT

for i in $(seq 10 11);
do
    search_dir=matrices
    for entry in "$search_dir"/*
    do
        echo "$entry"
        timeout $TIME_LIMIT make opt_sequential_cuda i=$i && ./opt_sequential $entry >> result_my.log;
        EXIT_CODE=$?
        echo "" >> result_my.log;
        RESULT=$(check_exit_status $EXIT_CODE "opt_sequential")
        echo "$RESULT" >> result_my.log;
        echo "-----------" >> result_my.log;
        echo "" >> result_my.log;

        if [ $EXIT_CODE -eq 0 ]; then
            timeout $TIME_LIMIT octave --eval "test_case('$entry', $i, false, $(num_to_bool $EXIT_CODE))"  >> result_oct.log;
            EXIT_CODE=$?
            echo "" >> result_oct.log;
            RESULT=$(check_exit_status $EXIT_CODE "opt_sequential")
            echo "$RESULT" >> result_oct.log;
            echo "-----------" >> result_oct.log;
            echo "" >> result_oct.log;
        fi
    done
    #make i=$i && ./opt_sequential 
done

exit 0

for i in $(seq 0 3);
do
    search_dir=matrices
    for entry in "$search_dir"/*
    do
        echo "$entry"
        timeout $TIME_LIMIT make i=$i && ./opt_sequential $entry >> result_my.log;
        EXIT_CODE=$?
        echo "" >> result_my.log;
        RESULT=$(check_exit_status $EXIT_CODE "opt_sequential")
        echo "$RESULT" >> result_my.log;
        echo "-----------" >> result_my.log;
        echo "" >> result_my.log;
    done
    #make i=$i && ./opt_sequential 
done
