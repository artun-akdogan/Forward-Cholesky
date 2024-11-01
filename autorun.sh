for i in $(seq 0 8);
do
    make CFLAGS=-DBUILD=$i && ./opt_sequential
done
