#ifndef OPT_SEQUENTIAL_UPPER
#define OPT_SEQUENTIAL_UPPER

#include <iostream>
#include <chrono>
#include <vector>
#include <set>
#include <algorithm>
#include <cstring>

#include "opt_sequential_common.h"

void upper_mat_init(const mattype num_rows,
                    const indtype num_lines,
                    indtype *&m_rows,
                    mattype *&m_cols,
                    dattype *&m_values,
                    const std::vector<std::vector<sparse_raw>> &matrix)
{
    std::vector<std::vector<sparse_raw>> rev_matrix(num_rows);
    for (mattype l = 0; l < num_rows; l++)
    {
        for (mattype i = 0; i < matrix[l].size(); i++)
        {
            rev_matrix[matrix[l][i].column].push_back({matrix[l][i].data, matrix[l][i].column, matrix[l][i].row});
        }
    }
    for (mattype l = 0; l < num_rows; l++)
    {
        std::sort(rev_matrix[l].begin(), rev_matrix[l].end(), [](const sparse_raw &lhs, const sparse_raw &rhs)
                  { return lhs.column < rhs.column; });
    }
    m_rows = new indtype[num_rows + 1];
    m_cols = new mattype[num_lines];
    m_values = new dattype[num_lines];

    m_rows[0] = 0;
    for (mattype l = 0; l < num_rows; l++)
    {
        for (mattype i = 0; i < rev_matrix[l].size(); i++)
        {
            m_cols[m_rows[l] + i] = rev_matrix[l][i].column;
            m_values[m_rows[l] + i] = rev_matrix[l][i].data;
        }
        m_rows[l + 1] = rev_matrix[l].size() + m_rows[l];
    }
}

void upper_res_init(const mattype num_rows,
                    indtype *&r_rows,
                    mattype *&r_cols,
                    dattype *&r_values,
                    const std::vector<std::set<mattype>> &rev_graph)
{

    auto end3 = std::chrono::high_resolution_clock::now();
    std::cout << "Reversing Graph" << std::endl;
    indtype tot_nonz = 0;
    std::vector<std::vector<mattype>> graph(num_rows);

    for (mattype i = 0; i < num_rows; i++)
    {
        auto fit = rev_graph[i].begin();

        for (; fit != rev_graph[i].end(); fit++)
        {
            graph[i].push_back(*fit);
            tot_nonz++;
        }
        tot_nonz++;
    }
    std::cout << "Total nonzero: " << tot_nonz << std::endl;

    auto end4 = std::chrono::high_resolution_clock::now();

    std::cout << "Normalized graph in " << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - end3).count() << " ms" << std::endl;

    std::cout << "Initializing result matrix" << std::endl;

    r_rows = new indtype[num_rows + 1];
    r_cols = new mattype[tot_nonz];
    r_values = new dattype[tot_nonz];

    std::memset(r_values, 0, tot_nonz * sizeof(dattype));

    r_rows[0] = 0;
    for (mattype l = 0; l < num_rows; l++)
    {
        r_cols[r_rows[l]] = l;
        for (mattype i = 0; i < graph[l].size(); i++)
        {
            r_cols[r_rows[l] + i + 1] = graph[l][i];
            // std::cout << l << " " << r_cols[r_rows[l] + i] <<  std::endl;
        }
        r_rows[l + 1] = r_rows[l] + graph[l].size() + 1;
    }
    std::cout << "R Mat " << r_rows[num_rows] << " " << tot_nonz << std::endl;
}

void upper_cholesky_calculate(const mattype num_rows,
                              const indtype *m_rows,
                              const mattype *m_cols,
                              const dattype *m_values,
                              const indtype *r_rows,
                              const mattype *r_cols,
                              dattype *r_values)
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
    for (mattype row = 0; row < num_rows; row++)
    {
        if (!(row % 1000))
        {
            std::cout << row << " " << num_rows << std::endl;
        }
        // std::cout << "->ok " << row << std::endl;
        indtype res_ind = r_rows[row];
        indtype mat_ind = m_rows[row];
        r_values[res_ind] = std::sqrt(m_values[mat_ind] - r_values[res_ind]);
        dattype res_diag = r_values[res_ind]; // r_values[r_rows[row]]
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
        // std::cout << "->ok " << res_diag << std::endl;
        for (indtype fi_ind = r_rows[row] + 1; fi_ind < r_rows[row + 1]; fi_ind++)
        {
            indtype tgt_ind = r_rows[r_cols[fi_ind]];
            for (indtype se_ind = fi_ind; se_ind < r_rows[row + 1]; se_ind++)
            {
                // std::cout << "--->ok " << r_cols[fi_ind] << " " << r_cols[se_ind] << std::endl;
                // tgt_ind = col_find(r_cols, r_cols[se_ind], tgt_ind, r_rows[r_cols[fi_ind] + 1]);
                tgt_ind = col_find_iter(r_cols, r_cols[se_ind], tgt_ind, r_rows[r_cols[fi_ind] + 1]);
                // tgt_ind = col_find_custom(r_cols, r_cols[se_ind], tgt_ind, r_rows[r_cols[fi_ind] + 1]);
                // tgt_ind = col_find_custom_iter(r_cols, r_cols[se_ind], tgt_ind, r_rows[r_cols[fi_ind] + 1]);
                if (tgt_ind == COLMAX)
                {
                    std::cout << "Error: " << row << std::endl;
                    exit(0);
                }
                r_values[tgt_ind] += r_values[fi_ind] * r_values[se_ind];
            }
        }
    }
}

#endif
