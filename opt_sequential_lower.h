/* Previous experiment, unused! */

#ifndef OPT_SEQUENTIAL_LOWER
#define OPT_SEQUENTIAL_LOWER

#include <iostream>
#include <chrono>

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
    for (mattype row = 0; row < num_rows; row++)
    {
        // std::cout << row << " " << num_rows << std::endl;
        /*if (!(row % 1000))
        {
            std::cout << row << " " << num_rows << std::endl;
        }*/
        for (indtype fi = r_rows[row]; fi < r_rows[row + 1] - 1; fi++)
        {
            dattype tot = 0, mat_val = 0, res_val = 0;

            // std::cout << "->ok " << fi << " " << r_cols[fi] << std::endl;
            //  for (idx_t l = m_rows[i]; l < m_rows[i + 1]; l++)
            if (topologicOrder[row] < topologicOrder[r_cols[fi]])
            {
                for (indtype a = r_rows[r_cols[fi]], b = r_rows[row]; a < r_rows[r_cols[fi] + 1] && b < fi;)
                {
                    // std::cout << "par: " << r_values[fi].column << " " << row << " " << r_values[a].column << " " << r_values[b].column << " " << tot << std::endl;
                    if (r_cols[a] < r_cols[b])
                    {
                        a++;
                    }
                    else if (r_cols[a] > r_cols[b])
                    {
                        b++;
                    }
                    else
                    {
                        tot += r_values[a] * r_values[b];
                        a++;
                        b++;
                    }
                }
            }
            indtype mat_it = col_find(m_cols, r_cols[fi], m_rows[row], m_rows[row + 1]);
            if (mat_it != COLMAX)
            {
                mat_val = m_values[mat_it];
            }
            // std::cout << "ok" << std::endl;
            indtype res_it = col_find(r_cols, r_cols[fi], r_rows[r_cols[fi]], r_rows[r_cols[fi] + 1]);
            if (res_it == COLMAX)
            {
                std::cout << "RError: " << row << std::endl;
                exit(0);
            }
            res_val = r_values[res_it];
            if ((mat_val - tot) / res_val != 0)
            {
                r_values[fi] = (mat_val - tot) / res_val;
            }
            // std::cout << "!ok" << std::endl;
        }
        // std::cout << "!!ok" << std::endl;
        dattype tot = 0, mat_val = 0;
        for (indtype a = r_rows[row]; a < r_rows[row + 1] - 1; a++)
        {
            tot += r_values[a] * r_values[a];
        }
        indtype mat_it = col_find(m_cols, row, m_rows[row], m_rows[row + 1]);

        // std::cout << "!!!ok" << std::endl;
        if (mat_it != COLMAX)
        {
            mat_val = m_values[mat_it];
        }
        // std::cout << "!!!ok" << r_rows[row + 1] - 1 << std::endl;
        if (mat_val - tot != 0)
        {
            r_values[r_rows[row + 1] - 1] = std::sqrt(mat_val - tot);
        }
        // // std::cout << mat_val << " " << tot  << std::endl;
        // std::cout << "!!!ok" << std::endl;
    }
}

#endif
