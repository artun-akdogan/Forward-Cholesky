for i in $(seq 0 7);
do
    make CFLAGS=-DBUILD=$i && ./opt_sequential
done
