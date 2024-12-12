#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <set>
#include <algorithm>
#include <filesystem>

#include "opt_sequential_common.h"

#if BUILD & 1
#define REORDER
#endif

#if BUILD & 2
#define UPPER
#endif

#if BUILD & 4
#define PARALLEL
#endif

#if BUILD & 8
#define CUDA
#endif

#ifdef REORDER
#include "SuiteSparse/COLAMD/Include/colamd.h"
#endif

#ifdef UPPER
#include "opt_sequential_upper_common.h"

#ifdef CUDA
#include "opt_sequential_upper_cuda.h"
#elif defined(PARALLEL)
#include "opt_sequential_upper_parallel.h"
#else
#include "opt_sequential_upper.h"
#endif

#else
#include "opt_sequential_lower_common.h"

#ifdef CUDA
static_assert(false, "Not implemented")
#elif defined(PARALLEL)
#include "opt_sequential_lower_parallel.h"
#else
#include "opt_sequential_lower.h"
#endif

#endif

void upper_mat_init_structure(const mattype num_rows,
                              const indtype num_lines,
                              indtype *&m_rows,
                              mattype *&m_cols,
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

    m_rows[0] = 0;
    for (mattype l = 0; l < num_rows; l++)
    {
        for (mattype i = 0; i < rev_matrix[l].size(); i++)
        {
            m_cols[m_rows[l] + i] = rev_matrix[l][i].column;
        }
        m_rows[l + 1] = rev_matrix[l].size() + m_rows[l];
    }
}

void operation_main(const char *matrix_name, bool save)
{
    auto beg = std::chrono::high_resolution_clock::now();
    //timer.start(0);

    std::ifstream file(matrix_name);
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
    int perm[num_row + 5];   /* note the size is N+1 */
    int stats[COLAMD_STATS]; /* for colamd and symamd output statistics */
    upper_mat_init_structure(num_row, num_lines, m_rows, m_cols, matrix);
    int ok = symamd(num_row, m_cols, m_rows, perm, (double *)NULL, stats, &calloc, &free);
    // symamd_report(stats);

    if (!ok)
    {
        printf("symamd error!\n");
        exit(1);
    }

    int *iperm = new int[num_row];
    for (int i = 0; i < num_row; i++)
    {
        iperm[perm[i]] = i;
    }

    if (save)
    {
        std::ofstream ofile("order.mtx");

        for (mattype i = 0; i < num_row; i++)
        {
            ofile << perm[i]+1 << "\n";
        }
        ofile.close();
    }

    delete[] m_rows;
    delete[] m_cols;

    std::vector<std::vector<sparse_raw>> reordered_matrix(num_row);
    for (mattype l = 0; l < num_row; l++)
    {
        for (mattype i = 0; i < matrix[l].size(); i++)
        { /*
             if(matrix[l][i].data<0.00000001 && matrix[l][i].data>-0.00000001){
                 std::cout << matrix[l][i].data;
                 exit(-1);
             }*/
            if (iperm[matrix[l][i].column] > iperm[matrix[l][i].row])
            {
                // std::cout << "Here: " <<iperm[matrix[l][i].row] << " " << iperm[matrix[l][i].column] << std::endl;
                reordered_matrix[iperm[matrix[l][i].column]].push_back({matrix[l][i].data, iperm[matrix[l][i].column], iperm[matrix[l][i].row]});
            }
            else
            {
                reordered_matrix[iperm[matrix[l][i].row]].push_back({matrix[l][i].data, iperm[matrix[l][i].row], iperm[matrix[l][i].column]});
            }
        }
    }

    for (mattype l = 0; l < num_row; l++)
    {
        std::sort(reordered_matrix[l].begin(), reordered_matrix[l].end(), [](const sparse_raw &lhs, const sparse_raw &rhs)
                  { return lhs.column < rhs.column; });
    }
    /*
        for (mattype l = 0; l < num_row; l++)
        {
            for (mattype i = 1; i < matrix[l].size(); i++){
                if (matrix[l][i-1].column>= matrix[l][i].column){
                    if (matrix[l][i-1].data!=matrix[l][i].data){
                        std::cout << "MError " << matrix[l][i-1].data << " " <<matrix[l][i].data << " " << matrix[l][i-1].row << " " <<matrix[l][i-1].column << " "<< matrix[l][i].column << std::endl;
                        exit(-1);
                    }else{
                        std::cout << "MError equal data" << std::endl;
                        exit(-1);
                    }
                }
            }
        }
        for (mattype l = 0; l < num_row; l++)
        {
            for (mattype i = 1; i < reordered_matrix[l].size(); i++){
                if (reordered_matrix[l][i-1].column>= reordered_matrix[l][i].column){
                    if (reordered_matrix[l][i-1].data!=reordered_matrix[l][i].data){
                        std::cout << "Error " << reordered_matrix[l][i-1].data << " " <<reordered_matrix[l][i].data << " " << reordered_matrix[l][i-1].row << " " <<reordered_matrix[l][i-1].column << " "<< reordered_matrix[l][i].column << std::endl;
                        exit(-1);
                    }else{
                        std::cout << "Error equal data" << std::endl;
                        exit(-1);
                    }
                }
            }
        }*/
    matrix = reordered_matrix;
    /*
        for (mattype l = 0; l < num_row; l++)
        {
            for (mattype k=0; k<matrix[l].size(); k++){
                std::cout << matrix[l][k].row +1 << " "<< matrix[l][k].column +1 << " "<< matrix[l][k].data << std::endl;
            }
        }*/

    auto end11 = std::chrono::high_resolution_clock::now();
    std::cout << "Reorder in " << std::chrono::duration_cast<std::chrono::milliseconds>(end11 - end1).count() << " ms" << std::endl;
#else
    auto end11 = std::chrono::high_resolution_clock::now();
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
                return;
            }
            else
            {
                break;
            }
        }
    }

    auto end2 = std::chrono::high_resolution_clock::now();

    std::cout << "Reverse graph created in " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end1).count() << " ms" << std::endl;

    std::cout << "Finding fills" << std::endl;

    indtype inserted = 0;
    int __approx__nnz = 0;

    for (mattype i = 0; i < num_row; i++)
    {
        __approx__nnz += rev_graph[i].size();
    }

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
                return;
            }
            if (rev_graph[*fit].insert(*sit).second)
            {
                __approx__nnz++;
            }
        }
        if (__approx__nnz > MAX_NNZ)
        {
            std::ofstream logfile("result.log", std::ios::app);
            logfile << matrix_name << " " << num_lines << " cannot computer more than " << MAX_NNZ << "\n";
            logfile.close();
            std::cout << "More than " << MAX_NNZ << " non-zeros. Skipping..." << std::endl;
            exit(1);
        }
    }

    auto end3 = std::chrono::high_resolution_clock::now();

    std::cout << "Found fills in " << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - end2).count() << " ms" << std::endl;

#ifndef UPPER
    std::cout << "Creating elimination tree (transitive reduction)" << std::endl;

    std::vector<std::vector<mattype>> tree(num_row);
    std::vector<std::vector<mattype>> topological;

    for (mattype i = 0; i < num_row; i++)
    {
        if (rev_graph[i].begin() != rev_graph[i].end())
        {
            tree[*(rev_graph[i].begin())].push_back(i);
        }
    }
    mattype *topologicOrder = new mattype[num_row];
    topologicalSort(tree, topological, topologicOrder);
    tree.clear();
    tree.shrink_to_fit();

    auto end31 = std::chrono::high_resolution_clock::now();

    std::cout << "Tree with topological depth " << topological.size() << " for rows " << num_row
              << " created in " << std::chrono::duration_cast<std::chrono::milliseconds>(end31 - end3).count() << "ms" << std::endl;
#endif

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
    delete[] topologicOrder;
#endif

    //std::cout << "Total iteration: " << total << std::endl;
    //total = 0;

    delete[] m_rows;
    delete[] m_cols;
    delete[] m_values;

    auto end5 = std::chrono::high_resolution_clock::now();

    std::cout << "Operation completed for " << r_rows[num_row] << " nonzeros in " << std::chrono::duration_cast<std::chrono::milliseconds>(end5 - end41).count() << " ms" << std::endl;

    std::cout << "Total time except disk io: " << std::chrono::duration_cast<std::chrono::milliseconds>(end5 - end11).count() << " ms" << std::endl;

    std::ofstream logfile("result.log", std::ios::app);
    logfile << matrix_name << " " << num_lines << " " << r_rows[num_row] << " " << std::chrono::duration_cast<std::chrono::milliseconds>(end5 - end1).count() << "\n";
    logfile.close();

    if (save)
    {
        std::ofstream ofile("result.mtx");
        ofile << num_row << " " << num_col << " " << r_rows[num_row] << "\n";

        for (mattype i = 0; i <= num_row - 1; i++)
        {
            for (indtype l = r_rows[i]; l < r_rows[i + 1]; l++)
            {

#ifdef UPPER
                ofile << r_cols[l] + 1 << " " << i + 1 << (r_values[l] < 0 ? " " : "  ") << r_values[l] << "\n";
#else
                ofile << i + 1 << " " << r_cols[l] + 1 << " " << r_values[l] << "\n";
#endif
            }
        }

        ofile.close();
    }

    delete[] r_rows;
    delete[] r_cols;
    delete[] r_values;

    auto end6 = std::chrono::high_resolution_clock::now();

    if (save)
    {
        std::cout << "File written in " << std::chrono::duration_cast<std::chrono::milliseconds>(end6 - end5).count() << " ms" << std::endl;
    }
    std::cout << "Program completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(end6 - beg).count() << " ms" << std::endl;

    //timer.stop(0);
    //timer.print_summary();
}

int main(int argc, char *argv[])
{
    //std::ofstream logfile("result.log", std::ios::app);
    std::cout <<
#ifdef REORDER
        "Reorder"
#else
        "Original"
#endif
            << " " <<
#ifdef CUDA
        "Cuda"
#elif defined(PARALLEL)
        "Parallel"
#else
        "Sequential"
#endif
            << " " <<
#ifdef UPPER
        "Upper"
#else
        "Lower"
#endif
            << std::endl;

    if (argc == 2)
    {
        std::cout << "\n\nStart for " << argv[1] << std::endl;
        operation_main(argv[1], false);
        return 0;
    }

    for (const auto &entry : std::filesystem::directory_iterator("matrices"))
    {
        std::cout << "\n\nStart for " << entry.path().c_str() << std::endl;
        operation_main(entry.path().c_str(), false);
    }

    return 0;
}