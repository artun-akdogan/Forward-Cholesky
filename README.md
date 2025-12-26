# Forward Cholesky Method For Sparse Matrices
This repository hosts the source code of the Forward method for symmetric positive-definite sparse matrices.


### Quick Run
Please ensure your PC has openBLAS, openMP and CUDA (if GPU version will be run) installed.
* sudo apt install libopenblas-dev libomp-dev

To run reordered parallel Forward algorithm, please run the following:
* make i=7 && ./opt_sequential matrices/bcsstk03.mtx
For unordered:
* make i=6 && ./opt_sequential matrices/bcsstk03.mtx


To run the GPU version of the Forward algorithm with reorder, please run the following (nvcc must exist, /usr/local/cuda/lib64 path is required):
* make opt_sequential_cuda i=11 && ./opt_sequential matrices/bcsstk03.mtx 
For unordered:
* make opt_sequential_cuda i=10 && ./opt_sequential matrices/bcsstk03.mtx 

### Testing
To obtain the relative reconstruction residude, please run the following (Second argument should be the number you run the build with):
* octave --eval "test_case('matrices/bcsstk03.mtx', 7, false, true)"

