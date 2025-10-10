#ifndef OPT_SEQUENTIAL_UPPER
#define OPT_SEQUENTIAL_UPPER

#include <iostream>
#include <chrono>
#include <omp.h>

#include "opt_sequential_common.h"

void upper_cholesky_calculate(const mattype num_rows,
                              const indtype *m_rows,
                              const mattype *m_cols,
                              const dattype *m_values,
                              const indtype *r_rows,
                              const mattype *r_cols,
                              dattype *r_values,
                              const mattype *topologicMatrix,
                              const mattype *topologicRows,
                              const mattype topologicDepth)
{ /*
     for (mattype row = 0; row < num_rows; row++){
         for(indtype l = m_rows[row]; l<m_rows[row+1]; l++){
             std::cout << row << " " << m_cols[l] << " " << m_values[l] << std::endl;
         }
     }
     for (mattype row = 0; row < num_rows; row++){
         for(indtype l = r_rows[row]; l<r_rows[row+1]; l++){
             std::cout << row << " " << r_cols[l] << " " << r_values[l] << std::endl;
         }
     }*/
    // mattype perc_bit = (1<< (25-__builtin_clz(num_rows))) -1;
    for (mattype cur_iter = 0; cur_iter < topologicDepth; cur_iter++)
    {   /*
           if((row & perc_bit) == perc_bit){
               std::cout << "Progress: " << (row + 1) * 100 / num_rows << "%\r" << std::flush;
           }*/
        /*if (!(row % 1000))
        {
            std::cout << row << " " << num_rows << std::endl;
        }*/
        // std::cout << "->ok " << row << std::endl;

#pragma omp parallel for firstprivate(m_values, m_rows, r_rows, m_cols, r_cols) //num_threads(4)
        for (mattype row_ind = topologicRows[cur_iter]; row_ind < topologicRows[cur_iter+1]; row_ind++)
        {
            mattype row = topologicMatrix[row_ind];
            indtype res_ind = r_rows[row];
            indtype mat_ind = m_rows[row];
            r_values[res_ind] = std::sqrt(m_values[mat_ind] - r_values[res_ind]);
            dattype res_diag = r_values[res_ind]; // r_values[r_rows[row]]
            // timer.start(1);
            for (res_ind++, mat_ind++; res_ind < r_rows[row + 1]; res_ind++)
            {
                dattype temp_mat = 0;
                if (mat_ind < m_rows[row + 1] && m_cols[mat_ind] == r_cols[res_ind])
                {
                    temp_mat = m_values[mat_ind];
                    mat_ind++;
                }
                r_values[res_ind] = (temp_mat - r_values[res_ind]) / res_diag;
            }
        }
        // timer.stop(1);
        //  std::cout << "->ok " << res_diag << std::endl;
        // timer.start(2);

        for (mattype row_ind = topologicRows[cur_iter]; row_ind < topologicRows[cur_iter+1]; row_ind++)
        {
            mattype row = topologicMatrix[row_ind];
#pragma omp parallel for firstprivate(r_rows, r_cols, row) //num_threads(4)
            for (indtype fi_ind = r_rows[row] + 1; fi_ind < r_rows[row + 1]; fi_ind++)
            {
                indtype tgt_ind = r_rows[r_cols[fi_ind]];
                // const indtype end_tgt_ind = r_rows[r_cols[fi_ind] + 1];
                for (indtype se_ind = fi_ind; se_ind < r_rows[row + 1]; se_ind++)
                {
                    // std::cout << "--->ok " << r_cols[fi_ind] << " " << r_cols[se_ind] << std::endl;
                    // tgt_ind = col_find(r_cols, r_cols[se_ind], tgt_ind, r_rows[r_cols[fi_ind] + 1]);
                    // tgt_ind = col_find_iter(r_cols, r_cols[se_ind], tgt_ind, r_rows[r_cols[fi_ind] + 1]);
                    // tgt_ind = col_find_custom(r_cols, r_cols[se_ind], tgt_ind, r_rows[r_cols[fi_ind] + 1]);
                    // tgt_ind = col_find_custom_iter(r_cols, r_cols[se_ind], tgt_ind, r_rows[r_cols[fi_ind] + 1]);

                    tgt_ind = col_find_iter(r_cols, r_cols[se_ind], tgt_ind, r_rows[r_cols[fi_ind] + 1]);
                    if (tgt_ind == COLMAX)
                    {
                        std::cout << "Error: " << row << std::endl;
                        exit(0);
                    }
                    r_values[tgt_ind] += r_values[fi_ind] * r_values[se_ind];
                }
            }
        }
        // timer.stop(2);
    }
}

#endif
