for i in $(seq 0 7);
do
    make i=$i && ./opt_sequential
done
