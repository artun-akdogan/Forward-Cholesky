#-----------------------------------------------------------------------------
# compile the COLAMD demo
#-----------------------------------------------------------------------------

default: opt_sequential

include SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk

I = -ISuiteSparse/COLAMD/Include -ISuiteSparse/SuiteSparse_config

C = $(CXX) $(CF) $(I)

library:
	( cd SuiteSparse/COLAMD/Lib ; $(MAKE) )

#------------------------------------------------------------------------------
# Create the demo program, run it, and compare the output
#------------------------------------------------------------------------------

dist:

opt_sequential: opt_sequential.cpp library
	$(C) $(CFLAGS) -DBUILD=$(i) -o opt_sequential opt_sequential.cpp SuiteSparse/COLAMD/Lib/libcolamd.a -lm -fopenmp

opt_sequential_upper_cuda.o: opt_sequential_upper_cuda.cu
	nvcc -c opt_sequential_upper_cuda.cu -o opt_sequential_upper_cuda.o

opt_sequential_cuda: opt_sequential_upper_cuda.o opt_sequential.cpp library
	$(C) $(CFLAGS) -DBUILD=$(i) -o opt_sequential.o -c opt_sequential.cpp SuiteSparse/COLAMD/Lib/libcolamd.a -lm
	$(C) opt_sequential_upper_cuda.o opt_sequential.o -o opt_sequential -L/usr/local/cuda/lib64 -lcudart


#------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#------------------------------------------------------------------------------

clean:
	- $(RM) $(CLEAN)

purge: distclean

distclean: clean
	- $(RM) colamd_example colamd_l_example
	- $(RM) my_colamd_example.out my_colamd_l_example.out
	- $(RM) -r *.dSYM
