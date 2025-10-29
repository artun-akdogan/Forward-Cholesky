#-----------------------------------------------------------------------------
# compile the COLAMD demo
#-----------------------------------------------------------------------------

default: opt_sequential

#include SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk

I = -ISuiteSparse/COLAMD/Include -ISuiteSparse/SuiteSparse_config

C = $(CXX) $(CF) $(I)

library:
	( cd SuiteSparse/SuiteSparse_config ; $(MAKE) library )
	( cd SuiteSparse/COLAMD ; $(MAKE) library )

#------------------------------------------------------------------------------
# Create the demo program, run it, and compare the output
#------------------------------------------------------------------------------

dist:

opt_sequential: opt_sequential.cpp library
	$(C) $(CFLAGS) -O3 -DBUILD=$(i) -o opt_sequential opt_sequential.cpp SuiteSparse/COLAMD/build/libcolamd.a SuiteSparse/SuiteSparse_config/build/libsuitesparseconfig.a -lm -fopenmp

opt_sequential_upper_cuda.o: opt_sequential_upper_cuda.cu
	nvcc -c opt_sequential_upper_cuda.cu -o opt_sequential_upper_cuda.o --expt-relaxed-constexpr --extended-lambda

opt_sequential_cuda: opt_sequential_upper_cuda.o opt_sequential.cpp library
	$(C) $(CFLAGS) -O3 -DBUILD=$(i) -o opt_sequential.o -c opt_sequential.cpp -lm
	$(C) -O3 opt_sequential_upper_cuda.o opt_sequential.o  SuiteSparse/COLAMD/build/libcolamd.a  SuiteSparse/SuiteSparse_config/build/libsuitesparseconfig.a -o opt_sequential -L/usr/local/cuda/lib64 -lcudart -fopenmp

opt_sequential_cuda_truba: opt_sequential_upper_cuda.o opt_sequential.cpp library
	$(C) $(CFLAGS) -O3 -DBUILD=$(i) -o opt_sequential.o -c opt_sequential.cpp -lm
	$(C) -O3 opt_sequential_upper_cuda.o opt_sequential.o  SuiteSparse/COLAMD/build/libcolamd.a  SuiteSparse/SuiteSparse_config/build/libsuitesparseconfig.a -o opt_sequential -L/arf/sw/lib/cuda/12.4/lib64 -lcudart -fopenmp


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
