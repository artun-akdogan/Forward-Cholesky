#ifndef OPT_SEQUENTIAL_LOWER_COMMON
#define OPT_SEQUENTIAL_LOWER_COMMON

#include <iostream>
#include <chrono>
#include <vector>
#include <set>
#include <algorithm>
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
            // std::cout << l << " " << graph[l][i] << std::endl;
        }
        r_cols[r_rows[l] + graph[l].size()] = l;
        // std::cout << l << " " << l << std::endl;
        r_rows[l + 1] = r_rows[l] + graph[l].size() + 1;
    }
    std::cout << "R Mat " << r_rows[num_rows] << " " << tot_nonz << std::endl;
}

#endif