#ifndef OPT_SEQUENTIAL_LOWER
#define OPT_SEQUENTIAL_LOWER

#include <iostream>
#include <chrono>
#include <vector>
#include <set>
#include <cstring>

#include "opt_sequential_common.h"

void lower_mat_init(const mattype num_rows,
                    const indtype num_lines,
                    indtype *&m_rows,
                    mattype *&m_cols,
                    dattype *&m_values,
                    const std::vector<std::vector<sparse_raw>> &matrix)
{
    m_rows = new indtype[num_rows + 1];
    m_cols = new mattype[num_lines];
    m_values = new dattype[num_lines];

    m_rows[0] = 0;
    for (mattype l = 0; l < num_rows; l++)
    {
        for (mattype i = 0; i < matrix[l].size(); i++)
        {
            m_cols[m_rows[l] + i] = matrix[l][i].column;
            m_values[m_rows[l] + i] = matrix[l][i].data;
        }
        m_rows[l + 1] = matrix[l].size() + m_rows[l];
    }
}

void lower_res_init(const mattype num_rows,
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
            graph[*fit].push_back(i);
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
        for (mattype i = 0; i < graph[l].size(); i++)
        {
            r_cols[r_rows[l] + i] = graph[l][i];
            //std::cout << l << " " << graph[l][i] << std::endl;
        }
        r_cols[r_rows[l] + graph[l].size()] = l;
        //std::cout << l << " " << l << std::endl;
        r_rows[l + 1] = r_rows[l] + graph[l].size() + 1;
    }
    std::cout << "R Mat " << r_rows[num_rows] << " " << tot_nonz << std::endl;
}

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
        if (!(row % 1000))
        {
            std::cout << row << " " << num_rows << std::endl;
        }
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
