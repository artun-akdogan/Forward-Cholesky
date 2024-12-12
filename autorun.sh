for i in $(seq 0 7);
do
    search_dir=matrices
    for entry in "$search_dir"/*
    do
        echo "$entry"
        (make i=$i && ./opt_sequential $entry >> result2.log) && (echo "------" >> result2.log) && (octave --eval "test_case('$entry', $i, false)"  >> result2.log);
        echo "-----------" >> result2.log;
        echo "" >> result2.log;
    done
    #make i=$i && ./opt_sequential 
done
