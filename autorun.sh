for i in $(seq 6 7);
do
    search_dir=matrices
    for entry in "$search_dir"/*
    do
        rm /tmp/result.mtx 2> /dev/null;
        rm /tmp/order.mtx 2> /dev/null;
        echo "$entry";
        timeout 3600 make i=$i && ./opt_sequential $entry >> result_my.log;
        echo "-----------" >> result_my.log;
        echo "" >> result_my.log;

        echo "$entry";
        echo "Matlab";
        (timeout 3600 matlab -batch "test_case('$entry', $i, false)"  >> result_mat.log) && (
        echo "$entry";
        echo "Octave";
        timeout 3600 octave --eval "test_case('$entry', $i, false)"  >> result_oct.log;
        echo "-----------" >> result_oct.log;
        echo "" >> result_oct.log;);
        echo "-----------" >> result_mat.log;
        echo "" >> result_mat.log;
    done
    #make i=$i && ./opt_sequential 
done
exit 0

for i in $(seq 2 7);
do
    search_dir=matrices
    for entry in "$search_dir"/*
    do
        echo "$entry"
        timeout 3600 make i=$i && ./opt_sequential $entry >> result_my.log;
        echo "-----------" >> result_my.log;
        echo "" >> result_my.log;
    done
    #make i=$i && ./opt_sequential 
done

for i in $(seq 10 11);
do
    search_dir=matrices
    for entry in "$search_dir"/*
    do
        echo "$entry"
        timeout 3600 make opt_sequential_cuda i=$i && ./opt_sequential $entry >> result_my.log;
        echo "-----------" >> result_my.log;
        echo "" >> result_my.log;
    done
    #make i=$i && ./opt_sequential 
done
