#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/copy.h>
#include <thrust/functional.h>
#include <thrust/sequence.h>
#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>

#include <thrust/execution_policy.h>

#define CUDA

#include "opt_sequential_upper_cuda.h"

struct columnCalculate

{
    const mattype row;
    const indtype *d_m_rows;
    const mattype *d_m_cols;
    const dattype *d_m_values;
    const indtype *d_r_rows;
    const mattype *d_r_cols;
    dattype *d_r_values;

    __host__ __device__
    columnCalculate(const mattype row, const indtype *d_m_rows, const mattype *d_m_cols,
                    const dattype *d_m_values, const indtype *d_r_rows, const mattype *d_r_cols,
                    dattype *d_r_values) : row(row), d_m_rows(d_m_rows), d_m_cols(d_m_cols), d_m_values(d_m_values),
                                           d_r_rows(d_r_rows), d_r_cols(d_r_cols), d_r_values(d_r_values) {}
    __device__ void operator()(const indtype &res_ind)
    {
        if (res_ind == d_r_values[d_r_rows[row]])
        {
            d_r_values[d_r_rows[row]] = std::sqrt(d_m_values[d_m_rows[row]] - d_r_values[d_r_rows[row]]);
            return;
        }
        __syncthreads();
        dattype res_diag = d_r_values[d_r_rows[row]];
        /*dattype res_diag = std::sqrt(d_m_values[d_m_rows[row]] - d_r_values[d_r_rows[row]]);

        if(res_ind==d_r_rows[row]){
            d_r_values[d_r_rows[row]] = res_diag;
            return;
        }*/

        dattype temp_mat = 0;
        indtype tgt_ind = col_find(d_m_cols, d_r_cols[res_ind], d_m_rows[row] + 1, d_m_rows[row + 1]);
        if (tgt_ind != COLMAX)
        {
            temp_mat = d_m_values[tgt_ind];
        }

        d_r_values[res_ind] = (temp_mat - d_r_values[res_ind]) / res_diag;
    }
};

struct forwardIteration
{

    const mattype row;
    const indtype *d_m_rows;
    const mattype *d_m_cols;
    const dattype *d_m_values;
    const indtype *d_r_rows;
    const mattype *d_r_cols;
    dattype *d_r_values;

    __host__ __device__
    forwardIteration(const mattype row, const indtype *d_m_rows, const mattype *d_m_cols,
                     const dattype *d_m_values, const indtype *d_r_rows, const mattype *d_r_cols,
                     dattype *d_r_values) : row(row), d_m_rows(d_m_rows), d_m_cols(d_m_cols), d_m_values(d_m_values),
                                            d_r_rows(d_r_rows), d_r_cols(d_r_cols), d_r_values(d_r_values) {}

    __device__ void operator()(const indtype &iter)
    {
        indtype max_iter = d_r_rows[row + 1] - d_r_rows[row] - 1;
        indtype max_iter_const = 4 * max_iter * max_iter + 4 * max_iter + 1;
        indtype i = (indtype)((max_iter << 1) + 1 - std::sqrt(max_iter_const - (iter << 3))) >> 1;
        indtype j = i + iter - ((i * (2 * max_iter - i + 1)) >> 1);
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
};

/*
void upper_cholesky_calculate(const mattype num_rows,
                              const indtype *m_rows,
                              const mattype *m_cols,
                              const dattype *m_values,
                              const indtype *r_rows,
                              const mattype *r_cols,
                              dattype *r_values)
{

    thrust::device_vector<indtype> __d_m_rows(m_rows, m_rows + num_rows + 1);
    thrust::device_vector<mattype> __d_m_cols(m_cols, m_cols + m_rows[num_rows]);
    thrust::device_vector<dattype> __d_m_values(m_values, m_values + m_rows[num_rows]);
    thrust::device_vector<indtype> __d_r_rows(r_rows, r_rows + num_rows + 1);
    thrust::device_vector<mattype> __d_r_cols(r_cols, r_cols + r_rows[num_rows]);
    thrust::device_vector<dattype> __d_r_values(r_rows[num_rows], 0);

    thrust::device_vector<mattype> __d_r_cache(r_rows[num_rows], -1);

    indtype *d_m_rows = thrust::raw_pointer_cast(__d_m_rows.data());
    mattype *d_m_cols = thrust::raw_pointer_cast(__d_m_cols.data());
    dattype *d_m_values = thrust::raw_pointer_cast(__d_m_values.data());
    indtype *d_r_rows = thrust::raw_pointer_cast(__d_r_rows.data());
    mattype *d_r_cols = thrust::raw_pointer_cast(__d_r_cols.data());
    dattype *d_r_values = thrust::raw_pointer_cast(__d_r_values.data());

    /*
    mattype *d_r_cache = thrust::raw_pointer_cast(__d_r_cache.data());

    thrust::counting_iterator<indtype> rows_begin_t(0);
    thrust::counting_iterator<indtype> rows_end_t(num_rows);

    thrust::for_each(
    rows_begin_t, rows_end_t,
    [=] __device__(const indtype &row)
    {
        indtype res_ind = d_r_rows[row];
        indtype mat_ind = d_m_rows[row];
        for (; res_ind < d_r_rows[row + 1]; res_ind++){
            if (mat_ind < d_m_rows[row + 1] && d_m_cols[mat_ind] == d_r_cols[res_ind])
            {
                d_r_cache[res_ind] = mat_ind;
                mat_ind++;
            }
        }
    });*//*

    for (mattype row = 0; row < num_rows; row++)
    {
        //if (!(row % 1000))
        //{
        //    std::cout << row << " " << num_rows << std::endl;
        //}
        // timer.start(1);
        //  std::cout << "->ok " << row << std::endl;
        // d_r_values[d_r_rows[row]] = std::sqrt(d_m_values[d_m_rows[row]] - d_r_values[d_r_rows[row]]);
        // timer.stop(1);
        // timer.start(2);
        // std::cout << "->ok " << d_r_rows[row] + 1 << " " << d_r_rows[row + 1] << std::endl;
        // std::cout << "->ok " << std::endl;
        indtype max_iter = r_rows[row + 1] - r_rows[row] - 1;
        indtype __max_iter = max_iter * (max_iter + 1) >> 1;
        thrust::counting_iterator<indtype> rows_begin_r(0);
        thrust::counting_iterator<indtype> rows_end_r(__max_iter);
        thrust::counting_iterator<indtype> rows_begin_c(r_rows[row]);
        thrust::counting_iterator<indtype> rows_end_c(r_rows[row + 1]);

        thrust::for_each(
            rows_begin_c, rows_end_c,
            [=] __device__(const indtype &res_ind)
            {
                if (res_ind == d_r_values[d_r_rows[row]])
                {
                    d_r_values[d_r_rows[row]] = std::sqrt(d_m_values[d_m_rows[row]] - d_r_values[d_r_rows[row]]);
                    return;
                }
                __syncthreads();
                dattype res_diag = d_r_values[d_r_rows[row]];
                //dattype res_diag = std::sqrt(d_m_values[d_m_rows[row]] - d_r_values[d_r_rows[row]]);

                //if(res_ind==d_r_rows[row]){
                //    d_r_values[d_r_rows[row]] = res_diag;
                //    return;
                //}

                dattype temp_mat = 0;

                indtype tgt_ind = col_find(d_m_cols, d_r_cols[res_ind], d_m_rows[row] + 1, d_m_rows[row + 1]);
                if (tgt_ind != COLMAX)
                {
                    temp_mat = d_m_values[tgt_ind];
                }
                //if (d_r_cache[res_ind] != -1)
                //{
                //    temp_mat = d_m_values[d_r_cache[res_ind]];
                //}

                d_r_values[res_ind] = (temp_mat - d_r_values[res_ind]) / res_diag;
            });
        // timer.stop(2);
        // timer.start(3);
        // std::cout << "->ok " << res_diag << std::endl;
        // std::cout << "--->ok " << std::endl;
        // std::cout << "->ok " << d_r_rows[row] + 1 << " " << d_r_rows[row + 1] << std::endl;
        thrust::for_each(
            rows_begin_r, rows_end_r,
            [=] __device__(const indtype &iter)
            {
                indtype max_iter = d_r_rows[row + 1] - d_r_rows[row] - 1;
                indtype max_iter_const = 4 * max_iter * max_iter + 4 * max_iter + 1;
                indtype i = (indtype)((max_iter << 1) + 1 - std::sqrt(max_iter_const - (iter << 3))) >> 1;
                indtype j = i + iter - ((i * (2 * max_iter - i + 1)) >> 1);
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
            });
        // timer.stop(3);
    }
    thrust::copy(__d_r_values.begin(), __d_r_values.end(), r_values);
}
*/

void upper_cholesky_calculate(const mattype num_rows,
                              const indtype *m_rows,
                              const mattype *m_cols,
                              const dattype *m_values,
                              const indtype *r_rows,
                              const mattype *r_cols,
                              dattype *r_values)
{

    thrust::device_vector<indtype> d_m_rows(m_rows, m_rows + num_rows + 1);
    thrust::device_vector<mattype> d_m_cols(m_cols, m_cols + m_rows[num_rows]);
    thrust::device_vector<dattype> d_m_values(m_values, m_values + m_rows[num_rows]);
    thrust::device_vector<indtype> d_r_rows(r_rows, r_rows + num_rows + 1);
    thrust::device_vector<mattype> d_r_cols(r_cols, r_cols + r_rows[num_rows]);
    thrust::device_vector<dattype> d_r_values(r_rows[num_rows], 0);

    for (mattype row = 0; row < num_rows; row++)
    {
        //if (!(row % 1000))
        //{
        //    std::cout << row << " " << num_rows << std::endl;
        //}
        //timer.start(1);
        // std::cout << "->ok " << row << std::endl;
        //d_r_values[d_r_rows[row]] = std::sqrt(d_m_values[d_m_rows[row]] - d_r_values[d_r_rows[row]]);
        //timer.stop(1);
        //timer.start(2);
        thrust::counting_iterator<indtype> rows_begin_c(r_rows[row]+1);
        thrust::counting_iterator<indtype> rows_end_c(r_rows[row + 1]);
        indtype max_iter = r_rows[row + 1] - r_rows[row] - 1;
        indtype __max_iter = max_iter * (max_iter + 1) >> 1;
        // std::cout << "->ok " << res_diag << std::endl;
        thrust::counting_iterator<indtype> rows_begin(0);
        thrust::counting_iterator<indtype> rows_end(__max_iter);
        // std::cout << "->ok " << d_r_rows[row] + 1 << " " << d_r_rows[row + 1] << std::endl;
        thrust::for_each(rows_begin_c, rows_end_c, columnCalculate(row, thrust::raw_pointer_cast(d_m_rows.data()), thrust::raw_pointer_cast(d_m_cols.data()), thrust::raw_pointer_cast(d_m_values.data()), thrust::raw_pointer_cast(d_r_rows.data()), thrust::raw_pointer_cast(d_r_cols.data()), thrust::raw_pointer_cast(d_r_values.data())));

        //timer.stop(2);
        //timer.start(3);
        // std::cout << "->ok " << d_r_rows[row] + 1 << " " << d_r_rows[row + 1] << std::endl;
        thrust::for_each(rows_begin, rows_end, forwardIteration(row, thrust::raw_pointer_cast(d_m_rows.data()), thrust::raw_pointer_cast(d_m_cols.data()), thrust::raw_pointer_cast(d_m_values.data()), thrust::raw_pointer_cast(d_r_rows.data()), thrust::raw_pointer_cast(d_r_cols.data()), thrust::raw_pointer_cast(d_r_values.data())));

        //timer.stop(3);
    }
    thrust::copy(d_r_values.begin(), d_r_values.end(), r_values);
}
