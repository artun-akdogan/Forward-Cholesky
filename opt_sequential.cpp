#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <set>
#include <algorithm>

#include "opt_sequential_common.h"

// #define UPPER
#define PARALLEL
#define REORDER

#ifdef REORDER
#include "SuiteSparse/COLAMD/Include/colamd.h"
#endif


#ifdef UPPER

#ifdef PARALLEL
#include "opt_sequential_upper_parallel.h"
#else
#include "opt_sequential_upper.h"
#endif

#else

#ifdef PARALLEL
#include "opt_sequential_lower_parallel.h"
#else
#include "opt_sequential_lower.h"
#endif

#endif


void lower_mat_init_structure(const mattype num_rows,
                    const indtype num_lines,
                    indtype *&m_rows,
                    mattype *&m_cols,
                    const std::vector<std::vector<sparse_raw>> &matrix)
{
    m_rows = new indtype[num_rows + 1];
    m_cols = new mattype[num_lines];

    m_rows[0] = 0;
    for (mattype l = 0; l < num_rows; l++)
    {
        for (mattype i = 0; i < matrix[l].size(); i++)
        {
            m_cols[m_rows[l] + i] = matrix[l][i].column;
        }
        m_rows[l + 1] = matrix[l].size() + m_rows[l];
    }
}

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
    std::vector<std::set<mattype>> rev_graph(num_row);

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

    auto end1 = std::chrono::high_resolution_clock::now();

    std::cout << "File read in " << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - beg).count() << " ms" << std::endl;

    indtype *m_rows, *r_rows;
    mattype *m_cols, *r_cols;
    dattype *m_values, *r_values;


#ifdef REORDER
    int perm [num_row+1] ;		/* note the size is N+1 */
    int stats [COLAMD_STATS] ;	/* for colamd and symamd output statistics */
    lower_mat_init_structure(num_row, num_lines, m_rows, m_cols, matrix);
    int ok = symamd (num_row, m_cols, m_rows, perm, (double *) NULL, stats, &calloc, &free) ;
    symamd_report (stats) ;

    if (!ok)
    {
	printf ("symamd error!\n") ;
	exit (1) ;
    }

    delete [] m_rows;
    delete [] m_cols;
    auto end11 = std::chrono::high_resolution_clock::now();
    std::cout << "Reorder in " << std::chrono::duration_cast<std::chrono::milliseconds>(end11 - end1).count() << " ms" << std::endl;
#endif

    for (mattype l = 0; l < num_row; l++)
    {
        for (mattype i = 0; i < matrix[l].size(); i++)
        {
            if (matrix[l][i].column < l)
            {
                rev_graph[matrix[l][i].column].insert(l);
                // graph[matrix[i][l].column].insert(i);
            }
            else if (matrix[l][i].column > l)
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
    mattype *topologicOrder = new mattype[num_row];
    topologicalSort(tree, topological, topologicOrder);
    tree.clear();
    tree.shrink_to_fit();

    auto end31 = std::chrono::high_resolution_clock::now();

    std::cout << "Tree with topological depth " << topological.size() << " for rows " << num_row
              << " created in " << std::chrono::duration_cast<std::chrono::milliseconds>(end31 - end3).count() << "ms" << std::endl;

#ifdef UPPER
    upper_mat_init(num_row, num_lines, m_rows, m_cols, m_values, matrix);
#else
    lower_mat_init(num_row, num_lines, m_rows, m_cols, m_values, matrix);
#endif

    matrix.clear();
    matrix.shrink_to_fit();

#ifdef UPPER
    upper_res_init(num_row, r_rows, r_cols, r_values, rev_graph);
#else
    lower_res_init(num_row, r_rows, r_cols, r_values, rev_graph);
#endif

    rev_graph.clear();
    rev_graph.shrink_to_fit();

    auto end41 = std::chrono::high_resolution_clock::now();

    std::cout << "Graph created in " << std::chrono::duration_cast<std::chrono::milliseconds>(end41 - end1).count() << " ms" << std::endl;

#ifdef UPPER
    upper_cholesky_calculate(num_row, m_rows, m_cols, m_values, r_rows, r_cols, r_values);
#else
    lower_cholesky_calculate(num_row, m_rows, m_cols, m_values, r_rows, r_cols, r_values, topological, topologicOrder);
#endif

    std::cout << "Total iteration: " << total << std::endl;

    delete[] m_rows;
    delete[] m_cols;
    delete[] m_values;
    delete[] topologicOrder;

    auto end5 = std::chrono::high_resolution_clock::now();

    std::cout << "Operation completed for " << r_rows[num_row] << " nonzeros in " << std::chrono::duration_cast<std::chrono::milliseconds>(end5 - end41).count() << " ms" << std::endl;

    std::ofstream ofile("result.mtx");
    ofile << num_row << " " << num_col << " " << r_rows[num_row] << "\n";

    for (mattype i = 0; i < num_row - 1; i++)
    {
        for (indtype l = r_rows[i]; l < r_rows[i + 1]; l++)
        {
            ofile << i + 1 << " " << r_cols[l] + 1 << " " << r_values[l] << "\n";
        }
    }

    ofile.close();

    delete[] r_rows;
    delete[] r_cols;
    delete[] r_values;

    auto end6 = std::chrono::high_resolution_clock::now();

    std::cout << "File written in " << std::chrono::duration_cast<std::chrono::milliseconds>(end6 - end5).count() << " ms" << std::endl;

    std::cout << "Program completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(end6 - beg).count() << " ms" << std::endl;

    return 0;
}
