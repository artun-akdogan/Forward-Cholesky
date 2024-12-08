#ifndef OPT_SEQUENTIAL_LOWER
#define OPT_SEQUENTIAL_LOWER

#include <iostream>
#include <chrono>
#include <omp.h>

#include "opt_sequential_common.h"

void lower_cholesky_calculate(const mattype num_rows,
                              const indtype *m_rows,
                              const mattype *m_cols,
                              const dattype *m_values,
                              const indtype *r_rows,
                              const mattype *r_cols,
                              dattype *r_values,
                              const std::vector<std::vector<mattype>> &topological,
                              const mattype *topologicOrder)
{
    for (mattype cur_iter = topological.size() - 1; cur_iter >= 0; cur_iter--)
    {
        if (!(cur_iter % 1000))
        {
            std::cout << cur_iter << " " << topological.size() << std::endl;
        } /*
         for (mattype row_ind = 0; row_ind < topological[cur_iter].size(); row_ind++)
         {
             std::cout << topological[cur_iter][row_ind] << " ";
         }
         std::cout << std::endl;*/

#pragma omp parallel for
        for (mattype row_ind = 0; row_ind < topological[cur_iter].size(); row_ind++)
        {
            mattype row = topological[cur_iter][row_ind];
            const indtype *__m_rows = m_rows;
            const mattype *__m_cols = m_cols;
            const dattype *__m_values = m_values;
            const indtype *__r_rows = r_rows;
            const mattype *__r_cols = r_cols;
            dattype *__r_values = r_values;
            const mattype *__topologicOrder = topologicOrder;
            // std::cout << row << " " << num_rows << std::endl;
            for (indtype fi = __r_rows[row]; fi < __r_rows[row + 1] - 1; fi++)
            {
                dattype tot = 0, mat_val = 0, res_val = 0;

                // std::cout << "->ok " << fi << " " << r_cols[fi] << std::endl;
                //  for (idx_t l = m_rows[i]; l < m_rows[i + 1]; l++)

                if (__topologicOrder[row] < __topologicOrder[__r_cols[fi]])
                {
                    for (indtype a = __r_rows[r_cols[fi]], b = __r_rows[row]; a < __r_rows[__r_cols[fi] + 1] && b < fi;)
                    {
                        // std::cout << "par: " << r_values[fi].column << " " << row << " " << r_values[a].column << " " << r_values[b].column << " " << tot << std::endl;
                        if (__r_cols[a] < __r_cols[b])
                        {
                            a++;
                        }
                        else if (__r_cols[a] > __r_cols[b])
                        {
                            b++;
                        }
                        else
                        {
                            tot += __r_values[a] * __r_values[b];
                            a++;
                            b++;
                        }
                    }
                }
                indtype mat_it = col_find(__m_cols, __r_cols[fi], __m_rows[row], __m_rows[row + 1]);
                if (mat_it != COLMAX)
                {
                    mat_val = __m_values[mat_it];
                }
                // std::cout << "ok" << std::endl;
                indtype res_it = col_find(__r_cols, __r_cols[fi], __r_rows[r_cols[fi]], __r_rows[r_cols[fi] + 1]);
                if (res_it == COLMAX)
                {
                    std::cout << "RError: " << row << std::endl;
                    exit(0);
                }
                res_val = __r_values[res_it];
                if ((mat_val - tot) / res_val != 0)
                {
                    __r_values[fi] = (mat_val - tot) / res_val;
                }
                // std::cout << "!ok" << std::endl;
            }
            // std::cout << "!!ok" << std::endl;
            dattype tot = 0, mat_val = 0;
            for (indtype a = __r_rows[row]; a < __r_rows[row + 1] - 1; a++)
            {
                tot += __r_values[a] * __r_values[a];
            }
            indtype mat_it = col_find(__m_cols, row, __m_rows[row], __m_rows[row + 1]);

            // std::cout << "!!!ok" << std::endl;
            if (mat_it != COLMAX)
            {
                mat_val = __m_values[mat_it];
            }
            // std::cout << "!!!ok" << r_rows[row + 1] - 1 << std::endl;
            if (mat_val - tot != 0)
            {
                __r_values[__r_rows[row + 1] - 1] = std::sqrt(mat_val - tot);
            }
            // // std::cout << mat_val << " " << tot  << std::endl;
            // std::cout << "!!!ok" << std::endl;
        }
    }
}

#endif
