#include <iostream>
#include <chrono>
#include <vector>
#include <set>
#include <algorithm>
#include <cstring>
#include <cuda_runtime.h>

#include "opt_sequential_upper_cuda.h"

#define CHECK_CUDA(call)                                                \
    cudaError_t err = call;                                             \
    if (err != cudaSuccess)                                             \
    {                                                                   \
        std::cerr << "CUDA error in " << __FILE__ << " at " << __LINE__ \
                  << ": " << cudaGetErrorString(err) << std::endl;      \
        exit(EXIT_FAILURE);                                             \
    }

__global__ void columnCalculate(const mattype row, const indtype *d_m_rows, const mattype *d_m_cols,
                                const dattype *d_m_values, const indtype *d_r_rows,
                                const mattype *d_r_cols, dattype *d_r_values)
{
    int iter = blockIdx.x * blockDim.x + threadIdx.x;
    dattype res_diag = sqrt((float)(d_m_values[d_m_rows[row]] - d_r_values[d_r_rows[row]]));
    if (iter == 0)
    {
        d_r_values[d_r_rows[row]] = res_diag;
        return;
    }

    indtype __max_iter_calc = d_r_rows[row + 1] - d_r_rows[row];
    dattype temp_mat = 0;
    indtype tgt_ind = col_find(d_m_cols, d_r_cols[d_r_rows[row] + iter], d_m_rows[row], d_m_rows[row + 1]);
    if (tgt_ind != COLMAX)
    {
        temp_mat = d_m_values[tgt_ind];
    }
    d_r_values[d_r_rows[row] + iter] = (temp_mat - d_r_values[d_r_rows[row] + iter]) / res_diag;
};

__global__ void forwardIteration(const mattype row, const indtype *d_m_rows, const mattype *d_m_cols,
                                 const dattype *d_m_values, const indtype *d_r_rows,
                                 const mattype *d_r_cols, dattype *d_r_values)
{
    int iter = blockIdx.x * blockDim.x + threadIdx.x;

    indtype max_iter = d_r_rows[row + 1] - d_r_rows[row] - 1;
    indtype __max_iter = max_iter * (max_iter + 1) >> 1;
    indtype max_iter_const = 4 * max_iter * max_iter + 4 * max_iter + 1;

    indtype i = (indtype)((max_iter << 1) + 1 - std::sqrt(max_iter_const - (iter << 3))) >> 1;
    indtype j = i + iter - ((i * (2 * max_iter - i + 1)) >> 1);
    // std::cout << iter << " " <<i << " " << j << " " << max_iter << std::endl;
    indtype fi_ind = d_r_rows[row] + 1 + i;
    indtype se_ind = d_r_rows[row] + 1 + j;
    indtype tgt_ind = col_find(d_r_cols, d_r_cols[se_ind], d_r_rows[d_r_cols[fi_ind]], d_r_rows[d_r_cols[fi_ind] + 1]);
    if (tgt_ind == COLMAX)
    {
        std::cout << "Error: " << row << std::endl;
        exit(0);
    }
    d_r_values[tgt_ind] += d_r_values[fi_ind] * d_r_values[se_ind];
};

void upper_cholesky_calculate(const mattype num_rows,
                              const indtype *m_rows,
                              const mattype *m_cols,
                              const dattype *m_values,
                              const indtype *r_rows,
                              const mattype *r_cols,
                              dattype *r_values)
{
    const int threadsPerBlock = 256;
    indtype *d_m_rows, *d_r_rows;
    mattype *d_m_cols, *d_r_cols;
    dattype *d_m_values, *d_r_values;

    CHECK_CUDA(cudaMalloc(&d_m_rows, num_rows * sizeof(indtype)));
    CHECK_CUDA(cudaMemcpy(d_m_rows, m_rows, num_rows * sizeof(indtype), cudaMemcpyHostToDevice));

    CHECK_CUDA(cudaMalloc(&d_m_cols, m_rows[num_rows] * sizeof(mattype)));
    CHECK_CUDA(cudaMemcpy(d_m_cols, m_cols, m_rows[num_rows] * sizeof(mattype), cudaMemcpyHostToDevice));

    CHECK_CUDA(cudaMalloc(&d_m_values, num_rows * sizeof(mattype)));
    CHECK_CUDA(cudaMemcpy(d_m_values, m_values, m_rows[num_rows] * sizeof(mattype), cudaMemcpyHostToDevice));

    CHECK_CUDA(cudaMalloc(&d_r_rows, num_rows * sizeof(indtype)));
    CHECK_CUDA(cudaMemcpy(d_r_rows, r_rows, num_rows * sizeof(indtype), cudaMemcpyHostToDevice));

    CHECK_CUDA(cudaMalloc(&d_r_cols, r_rows[num_rows] * sizeof(mattype)));
    CHECK_CUDA(cudaMemcpy(d_r_cols, r_cols, r_rows[num_rows] * sizeof(mattype), cudaMemcpyHostToDevice));

    CHECK_CUDA(cudaMalloc(&d_r_values, num_rows * sizeof(mattype)));

    for (mattype row = 0; row < num_rows; row++)
    {
        if (!(row % 1000))
        {
            std::cout << row << " " << num_rows << std::endl;
        }
        const int blocksPerGrid = (r_rows[row + 1] - r_rows[row] + threadsPerBlock - 1) / threadsPerBlock;
        columnCalculate<<<blocksPerGrid, threadsPerBlock>>>(row, d_m_rows, d_m_cols, d_m_values, d_r_rows, d_r_cols, d_r_values);
        forwardIteration<<<blocksPerGrid, threadsPerBlock>>>(row, d_m_rows, d_m_cols, d_m_values, d_r_rows, d_r_cols, d_r_values);
    }

    CHECK_CUDA(cudaPeekAtLastError());
    CHECK_CUDA(cudaMemcpy(r_values, d_r_values, r_rows[num_rows] * sizeof(mattype), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaFree(d_m_rows));
    CHECK_CUDA(cudaFree(d_m_cols));
    CHECK_CUDA(cudaFree(d_m_values));
    CHECK_CUDA(cudaFree(d_r_rows));
    CHECK_CUDA(cudaFree(d_r_cols));
    CHECK_CUDA(cudaFree(d_r_values));
}
