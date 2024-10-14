#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <set>
#include <algorithm>

#include "opt_sequential_common.h"
#include "opt_sequential_lower.h"
// #include "opt_sequential_upper.h"

struct sparse_raw
{
    dattype data;
    mattype row;
    mattype column;
};

int main()
{
    auto beg = std::chrono::high_resolution_clock::now();

    std::ifstream file("ct20stif.mtx");
    // std::ifstream file("bcsstk03.mtx");
    mattype num_row, num_col;
    indtype num_lines;

    // Ignore comments headers
    while (file.peek() == '%')
        file.ignore(2048, '\n');

    // Read number of rows and columns
    file >> num_row >> num_col >> num_lines;

    std::vector<std::vector<sparse_raw>> matrix(num_row);
    std::vector<std::vector<sparse_raw>> result(num_row);
    std::vector<std::set<mattype>> rev_graph(num_row);
    std::vector<std::vector<mattype>> graph(num_row);

    for (indtype l = 0; l < num_lines; l++)
    {
        dattype data;
        mattype row, col;
        file >> row >> col >> data;
        matrix[row - 1].push_back({data, row - 1, col - 1});
    }

    file.close();

    for (mattype l = 0; l < num_row; l++)
    {
        std::sort(matrix[l].begin(), matrix[l].end(), [](const sparse_raw &lhs, const sparse_raw &rhs)
                  { return lhs.column < rhs.column; });
    }

    indtype *m_rows = new indtype[num_row + 1];
    mattype *m_cols = new mattype[num_lines];
    dattype *m_values = new dattype[num_lines];

    m_rows[0] = 0;
    for (mattype l = 0; l < num_row; l++)
    {
        for (mattype i = 0; i < matrix[l].size(); i++)
        {
            m_cols[m_rows[l] + i] = matrix[l][i].column;
            m_values[m_rows[l] + i] = matrix[l][i].data;
        }
        m_rows[l + 1] = matrix[l].size() + m_rows[l];
    }

    std::cout << "Total nonzero " << num_lines << std::endl;
    auto end1 = std::chrono::high_resolution_clock::now();

    std::cout << "File read in " << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - beg).count() << " ms" << std::endl;

    auto end11 = std::chrono::high_resolution_clock::now();
    std::cout << "Reorder in " << std::chrono::duration_cast<std::chrono::milliseconds>(end11 - end1).count() << " ms" << std::endl;

    for (mattype i = 0; i < num_row; i++)
    {
        for (indtype l = m_rows[i]; l < m_rows[i + 1]; l++)
        {
            if (m_cols[l] < i)
            {
                rev_graph[m_cols[l]].insert(i);
                // graph[matrix[i][l].column].insert(i);
            }
            else if (m_cols[l] > i)
            {
                std::cout << "wrong" << std::endl;
                return 0;
            }
            else
            {
                break;
            }
        }
    }
    auto fit = rev_graph[0].begin();

    auto end2 = std::chrono::high_resolution_clock::now();

    std::cout << "Reverse graph created in " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end1).count() << " ms" << std::endl;

    std::cout << "Finding fills" << std::endl;

    indtype inserted = 0;

    for (mattype i = 0; i < num_row; i++)
    {
        auto fit = rev_graph[i].begin();

        if (fit == rev_graph[i].end())
        {
            continue;
        }

        auto sit = fit;
        for (sit++; sit != rev_graph[i].end(); sit++)
        {
            if (*fit == i || *sit == i)
            {
                std::cout << "wrong" << std::endl;
                return 0;
            }
            rev_graph[*fit].insert(*sit);
        }
    }

    auto end3 = std::chrono::high_resolution_clock::now();

    std::cout << "Found fills in " << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - end2).count() << " ms" << std::endl;

    std::cout << "Creating elimination tree (transitive reduction)" << std::endl;

    std::vector<std::vector<mattype>> tree(num_row);
    std::vector<std::vector<mattype>> topological;

    for (mattype i = 0; i < num_row; i++)
    {
        tree[*(rev_graph[i].begin())].push_back(i);
    }
    topologicalSort(tree, topological);
    tree.clear();
    tree.shrink_to_fit();

    auto end31 = std::chrono::high_resolution_clock::now();

    std::cout << "Tree with topological depth " << topological.size() << " for rows " << num_row
              << " created in " << std::chrono::duration_cast<std::chrono::milliseconds>(end31 - end3).count() << "ms" << std::endl;

    std::cout << "Reversing Graph" << std::endl;
    indtype tot_nonz = 0;

    for (mattype i = 0; i < num_row; i++)
    {
        auto fit = rev_graph[i].begin();

        for (; fit != rev_graph[i].end(); fit++)
        {
            graph[*fit].push_back(i);
            tot_nonz++;
        }
        tot_nonz++;
    }
    rev_graph.clear();
    rev_graph.shrink_to_fit();
    std::cout << "Total nonzero: " << tot_nonz << std::endl;

    auto end4 = std::chrono::high_resolution_clock::now();

    std::cout << "Normalized graph in " << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - end3).count() << " ms" << std::endl;

    std::cout << "Initializing result matrix" << std::endl;

    indtype *r_rows = new indtype[num_row + 1];
    mattype *r_cols = new mattype[tot_nonz];
    dattype *r_values = new dattype[tot_nonz];

    r_rows[0] = 0;
    for (mattype l = 0; l < num_row; l++)
    {
        for (mattype i = 0; i < graph[l].size(); i++)
        {
            r_cols[r_rows[l] + i] = graph[l][i];
        }
        r_cols[r_rows[l] + graph[l].size()] = l;
        r_rows[l + 1] = r_rows[l] + graph[l].size() + 1;
    }

    std::cout << "R Mat " << r_rows[num_row] << " " << tot_nonz << std::endl;

    std::cout << "Graph created in " << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - end1).count() << " ms" << std::endl;

    auto end41 = std::chrono::high_resolution_clock::now();

    lower_cholesky_calculate(num_row, m_rows, m_cols, m_values, r_rows, r_cols, r_values);

    auto end5 = std::chrono::high_resolution_clock::now();

    std::cout << "Operation completed for " << tot_nonz << " nonzeros in " << std::chrono::duration_cast<std::chrono::milliseconds>(end5 - end41).count() << " ms" << std::endl;

    std::ofstream ofile("result.mtx");
    ofile << num_row << " " << num_col << " " << tot_nonz << "\n";

    for (mattype i = 0; i < num_row - 1; i++)
    {
        for (indtype l = r_rows[i]; l < r_rows[i + 1]; l++)
        {
            ofile << i + 1 << " " << r_cols[l] + 1 << " " << r_values[l] << "\n";
        }
    }

    ofile.close();

    auto end6 = std::chrono::high_resolution_clock::now();

    std::cout << "File written in " << std::chrono::duration_cast<std::chrono::milliseconds>(end6 - end5).count() << " ms" << std::endl;

    std::cout << "Program completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(end6 - beg).count() << " ms" << std::endl;

    return 0;
}
