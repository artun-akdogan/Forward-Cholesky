#include <thrust/device_vector.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cooperative_groups.h>

namespace cg = cooperative_groups;

#define CUDA

#include "opt_sequential_upper_cuda.h"

/*
__device__ double sqrt_newton(double x)
{
    if (x < 0.0)
        return NAN;
    if (x == 0.0)
        return 0.0;

    // Get an approximate 1/sqrt(x) using single-precision rsqrtf
    float approx = rsqrtf((float)x); // still fast, even if x is double
    double y = (double)approx;

    // Perform Newton-Raphson iterations to refine
    // First iteration
    y = y * (1.5 - 0.5 * x * y * y);
    // Second iteration (optional, for very high precision)
    y = y * (1.5 - 0.5 * x * y * y);

    // Now compute sqrt(x)
    double result = x * y;
    return result;
}
*/

__global__ void upper_cholesky_calculate_algorithm(
    const mattype num_rows, const indtype *d_m_rows, const mattype *d_m_cols,
    const dattype *d_m_values, const indtype *d_r_rows, const mattype *d_r_cols,
    dattype *d_r_values, const mattype *d_topologicMatrix, const mattype *d_topologicRows,
    const mattype d_topologicDepth)
{
    cg::grid_group grid = cg::this_grid();

    int _idx = threadIdx.x + blockIdx.x * blockDim.x;

    for (mattype cur_iter = 0; cur_iter < d_topologicDepth; cur_iter++)
    {
        // Step 1: Process first element only

        indtype total_work = d_topologicRows[cur_iter + 1] - d_topologicRows[cur_iter];
        for (int idx = _idx; idx < total_work; idx += grid.size())
        {
            mattype row = d_topologicMatrix[idx + d_topologicRows[cur_iter]];
            d_r_values[d_r_rows[row]] = sqrt(d_m_values[d_m_rows[row]] - d_r_values[d_r_rows[row]]);
        }
        grid.sync(); // Sync all threads before moving to step 2

        total_work = 0;
        for (mattype row_ind = d_topologicRows[cur_iter]; row_ind < d_topologicRows[cur_iter + 1]; row_ind++)
        {
            mattype row = d_topologicMatrix[row_ind];
            total_work += d_r_rows[row + 1] - d_r_rows[row] - 1;
        }

        for (int gid = _idx; gid < total_work; gid += grid.size())
        {
            indtype acc = 0;
            for (mattype row_ind = d_topologicRows[cur_iter]; row_ind < d_topologicRows[cur_iter + 1]; row_ind++)
            {
                mattype row = d_topologicMatrix[row_ind];
                indtype max_iter = d_r_rows[row + 1] - d_r_rows[row] - 1;
                if (gid < acc + max_iter)
                {
                    int idx = gid - acc;

                    dattype res_diag = d_r_values[d_r_rows[row]];
                    dattype temp_mat = 0;
                    indtype tgt_ind = col_find(d_m_cols, d_r_cols[idx + d_r_rows[row] + 1], d_m_rows[row] + 1, d_m_rows[row + 1]);
                    if (tgt_ind != COLMAX)
                    {
                        temp_mat = d_m_values[tgt_ind];
                    }

                    d_r_values[idx + d_r_rows[row] + 1] = (temp_mat - d_r_values[idx + d_r_rows[row] + 1]) / res_diag;
                    break;
                }
                acc += max_iter;
            }
        }
        grid.sync(); // Sync before moving to step 3

        for (mattype row_ind = d_topologicRows[cur_iter]; row_ind < d_topologicRows[cur_iter + 1]; row_ind++)
        {
            mattype row = d_topologicMatrix[row_ind];
            indtype max_iter = d_r_rows[row + 1] - d_r_rows[row] - 1;
            // indtype fwd_max_iter = max_iter * (max_iter + 1) >> 1;
            indtype fwd_max_iter = max_iter * max_iter;

            for (int idx = _idx; idx < fwd_max_iter; idx += grid.size())
            {
                /*
                indtype max_iter_const = 4 * max_iter * max_iter + 4 * max_iter + 1;
                indtype i = (indtype)((max_iter << 1) + 1 - sqrtf(max_iter_const - (idx << 3))) >> 1;
                indtype j = i + idx - ((i * (2 * max_iter - i + 1)) >> 1);
                */
                indtype i = idx / max_iter;
                indtype j = idx % max_iter;
                if (i > j)
                    continue;
                // std::cout << iter << " " <<i << " " << j << " " << max_iter << std::endl;
                indtype fi_ind = d_r_rows[row] + 1 + i;
                indtype se_ind = d_r_rows[row] + 1 + j;
                indtype tgt_ind = d_r_rows[d_r_cols[fi_ind]];
                tgt_ind = col_find(d_r_cols, d_r_cols[se_ind], tgt_ind, d_r_rows[d_r_cols[fi_ind] + 1]);
                if (tgt_ind == COLMAX)
                {
                    printf("Error: %d\n", row);
                }
                d_r_values[tgt_ind] += d_r_values[fi_ind] * d_r_values[se_ind];
            }
            grid.sync(); // Final sync before the next iteration
        }
    }
}

void upper_cholesky_calculate(mattype num_rows,
                              const indtype *m_rows,
                              const mattype *m_cols,
                              const dattype *m_values,
                              const indtype *r_rows,
                              const mattype *r_cols,
                              dattype *r_values,
                              const mattype *topologicMatrix,
                              const mattype *topologicRows,
                              mattype topologicDepth)
{

    int supportsCoopLaunch = 0;
    cudaDeviceGetAttribute(&supportsCoopLaunch, cudaDevAttrCooperativeLaunch, 0);
    if (!supportsCoopLaunch)
    {
        std::cout << "Error: Your GPU does not support cooperative kernels!" << std::endl;
        abort();
    }
    /*
int maxBlocks;
cudaDeviceGetAttribute(&maxBlocks, cudaDevAttrMaxBlocksPerMultiprocessor, 0);
std::cout << "Max blocks per SM: " << maxBlocks << std::endl;
int maxThreadsPerSM, numSMs, maxThreads;
cudaDeviceGetAttribute(&maxThreadsPerSM, cudaDevAttrMaxThreadsPerMultiProcessor, 0);
cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, 0);
*/
    int maxThreadsPerBlock, maxThreadsPerSM, numSMs;
    cudaDeviceGetAttribute(&maxThreadsPerBlock, cudaDevAttrMaxThreadsPerBlock, 0);
    cudaDeviceGetAttribute(&maxThreadsPerSM, cudaDevAttrMaxThreadsPerMultiProcessor, 0);
    cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, 0);

    size_t free_mem, total_mem;
    cudaMemGetInfo(&free_mem, &total_mem);

    size_t required = sizeof(indtype) * (num_rows + 1) * 2 +
                      sizeof(mattype) * (m_rows[num_rows] + r_rows[num_rows] + topologicRows[topologicDepth] + topologicDepth + 1) +
                      sizeof(dattype) * (m_rows[num_rows], r_rows[num_rows]);
    if (required >= free_mem)
    {
        std::cout << "Not enough GPU memory" << std::endl;
        abort();
    }

    indtype max_row = 0;
    for (int i = 0; i < num_rows; i++)
    {
        if (max_row < (r_rows[i + 1] - r_rows[i]))
        {
            max_row = r_rows[i + 1] - r_rows[i];
        }
    }
    // indtype __max_iter = max_row * (max_row - 1) >> 1;
    indtype __max_iter = max_row * max_row;

    thrust::device_vector<indtype> d_m_rows(m_rows, m_rows + num_rows + 1);
    thrust::device_vector<mattype> d_m_cols(m_cols, m_cols + m_rows[num_rows]);
    thrust::device_vector<dattype> d_m_values(m_values, m_values + m_rows[num_rows]);
    thrust::device_vector<indtype> d_r_rows(r_rows, r_rows + num_rows + 1);
    thrust::device_vector<mattype> d_r_cols(r_cols, r_cols + r_rows[num_rows]);
    thrust::device_vector<dattype> d_r_values(r_rows[num_rows], 0);
    thrust::device_vector<mattype> d_topologicMatrix(topologicMatrix, topologicMatrix + topologicRows[topologicDepth]);
    thrust::device_vector<mattype> d_topologicRows(topologicRows, topologicRows + topologicDepth + 1);

    indtype *raw_d_m_rows = thrust::raw_pointer_cast(d_m_rows.data());
    mattype *raw_d_m_cols = thrust::raw_pointer_cast(d_m_cols.data());
    dattype *raw_d_m_values = thrust::raw_pointer_cast(d_m_values.data());
    indtype *raw_d_r_rows = thrust::raw_pointer_cast(d_r_rows.data());
    mattype *raw_d_r_cols = thrust::raw_pointer_cast(d_r_cols.data());
    dattype *raw_d_r_values = thrust::raw_pointer_cast(d_r_values.data());
    mattype *raw_d_topologicMatrix = thrust::raw_pointer_cast(d_topologicMatrix.data());
    mattype *raw_d_topologicRows = thrust::raw_pointer_cast(d_topologicRows.data());

    void *kernelArgs[] = {&num_rows, &raw_d_m_rows, &raw_d_m_cols, &raw_d_m_values, &raw_d_r_rows, &raw_d_r_cols, &raw_d_r_values, &raw_d_topologicMatrix, &raw_d_topologicRows, &topologicDepth};

    int threads = maxThreadsPerBlock;

    // Compute optimal number of blocks
    int maxTotalThreads = maxThreadsPerSM * numSMs;
    std::cout << "Running with " << maxTotalThreads << " threads" << std::endl;
    int blocks = min((__max_iter + threads - 1) / threads, maxTotalThreads / threads);

    cudaError_t err = cudaLaunchCooperativeKernel((void *)upper_cholesky_calculate_algorithm, dim3(blocks), dim3(threads), kernelArgs);
    if (err != cudaSuccess)
    {
        std::cout << "CUDA Error: " << cudaGetErrorString(err) << " (" << err << ")" << std::endl;
        abort();
    }

    cudaDeviceSynchronize();

    cudaError_t last = cudaGetLastError();
    if (last != cudaSuccess)
    {
        std::cout << "Runtime kernel error: "
                  << cudaGetErrorString(last)
                  << std::endl;
        abort();
    }

    thrust::copy(d_r_values.begin(), d_r_values.end(), r_values);
}
