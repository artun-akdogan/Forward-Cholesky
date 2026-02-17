#include <cuda_runtime.h>
#include <cusparse.h>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <set>
#include <algorithm>
#include <filesystem>
#include <limits>
#include <string>

#include "opt_sequential_common.h"
#include "opt_sequential_lower_common.h"


#if BUILD & 1
#define REORDER
#endif

#define CUDA_CHECK(x)       \
    if ((x) != cudaSuccess) \
        throw std::runtime_error(std::string("CUDA error ") + std::to_string(x));

#define CUSPARSE_CHECK(x)               \
    if ((x) != CUSPARSE_STATUS_SUCCESS) \
        throw std::runtime_error(std::string("cuSPARSE error ") + std::to_string(x));

void compute_numeric_L(
    int num_row,
    const int *m_rows,
    const int *m_cols,
    const double *m_values,
    const int *r_rows,
    const int *r_cols,
    double *r_values
)
{
    //std::cout << "ok1" << std::endl;
    int nnzL = r_rows[num_row];

    //std::cout << "ok1" << std::endl;
    // initialize r_values
    for (int i = 0; i < num_row; i++)
    {
        int a_ptr = m_rows[i];
        int a_end = m_rows[i + 1];

        for (int l_ptr = r_rows[i]; l_ptr < r_rows[i + 1]; l_ptr++)
        {
            int j = r_cols[l_ptr];

            // advance A pointer until >= j
            while (a_ptr < a_end && m_cols[a_ptr] < j)
                a_ptr++;

            if (a_ptr < a_end && m_cols[a_ptr] == j)
                r_values[l_ptr] = m_values[a_ptr];
            else
                r_values[l_ptr] = 0.0; // fill
        }
    }

    //std::cout << "ok3" << std::endl;

    int *d_r_rows, *d_r_cols;
    double *d_r_values;
    //std::cout << "ok" << std::endl;

    CUDA_CHECK(cudaMalloc(&d_r_rows, (num_row + 1) * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_r_cols, nnzL * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_r_values, nnzL * sizeof(double)));

    CUDA_CHECK(cudaMemcpy(d_r_rows, r_rows,
                          (num_row + 1) * sizeof(int),
                          cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMemcpy(d_r_cols, r_cols,
                          nnzL * sizeof(int),
                          cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMemcpy(d_r_values, r_values,
                          nnzL * sizeof(double),
                          cudaMemcpyHostToDevice));
    //std::cout << "ok4" << std::endl;

    // ---------------------------------------
    // Step 4: cuSPARSE setup
    // ---------------------------------------
    cusparseHandle_t handle;
    CUSPARSE_CHECK(cusparseCreate(&handle));
    //std::cout << "ok" << std::endl;

    cusparseMatDescr_t descr;
    CUSPARSE_CHECK(cusparseCreateMatDescr(&descr));
    CUSPARSE_CHECK(cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO));
    CUSPARSE_CHECK(cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL));
    CUSPARSE_CHECK(cusparseSetMatFillMode(descr, CUSPARSE_FILL_MODE_LOWER));
    CUSPARSE_CHECK(cusparseSetMatDiagType(descr, CUSPARSE_DIAG_TYPE_NON_UNIT));
    //std::cout << "ok" << std::endl;

    csric02Info_t info;
    CUSPARSE_CHECK(cusparseCreateCsric02Info(&info));
    //std::cout << "ok" << std::endl;

    int bufferSize;
    CUSPARSE_CHECK(
        cusparseDcsric02_bufferSize(
            handle,
            num_row,
            nnzL,
            descr,
            d_r_values,
            d_r_rows,
            d_r_cols,
            info,
            &bufferSize));

    void *buffer;
    CUDA_CHECK(cudaMalloc(&buffer, bufferSize));
    //std::cout << "ok5" << std::endl;

    CUSPARSE_CHECK(
        cusparseDcsric02_analysis(
            handle,
            num_row,
            nnzL,
            descr,
            d_r_values,
            d_r_rows,
            d_r_cols,
            info,
            CUSPARSE_SOLVE_POLICY_NO_LEVEL,
            buffer));

    //std::cout << "ok6" << std::endl;

    // Numeric factorization
    CUSPARSE_CHECK(
        cusparseDcsric02(
            handle,
            num_row,
            nnzL,
            descr,
            d_r_values,
            d_r_rows,
            d_r_cols,
            info,
            CUSPARSE_SOLVE_POLICY_NO_LEVEL,
            buffer));

    //std::cout << "ok7" << std::endl;


    CUDA_CHECK(cudaMemcpy(r_values,
                          d_r_values,
                          nnzL * sizeof(double),
                          cudaMemcpyDeviceToHost));

    //std::cout << "ok8" << std::endl;
    cusparseDestroyCsric02Info(info);
    cusparseDestroyMatDescr(descr);
    cusparseDestroy(handle);

    cudaFree(buffer);
    cudaFree(d_r_values);
    cudaFree(d_r_cols);
    cudaFree(d_r_rows);
}

Timer timer(5);

const char * rsplit(const char * start, const char delim){
    const char *newstart = start;
    while((*start)!='\0'){
        if((*start)==delim){
            newstart = start;
        }
        start++;
    }
    return newstart;
}

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
    // timer.start(0);

    std::ifstream file(matrix_name);
    // std::ifstream file("bcsstk03.mtx");
    mattype num_row, num_col;
    indtype init_num_lines, num_lines=0;

    // Ignore comments headers
    while (file.peek() == '%')
        file.ignore(2048, '\n');

    // Read number of rows and columns
    file >> num_row >> num_col >> init_num_lines;

    std::vector<std::vector<sparse_raw>> matrix(num_row);
    std::vector<std::set<mattype>> rev_graph(num_row);

    for (indtype l = 0; l < init_num_lines; l++)
    {
        dattype data;
        mattype row, col;
        file >> row >> col >> data;
        if (std::abs(data) >= std::numeric_limits<dattype>::epsilon()){
            matrix[row - 1].push_back({data, row - 1, col - 1});
            num_lines++;
        }
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
    
    std::ifstream ifile("result/order.mtx");
    for (mattype i = 0; i < num_row; i++)
    {
        ifile >> perm[i];
        perm[i]--;
    }
    ifile.close();

    int *iperm = new int[num_row];
    for (int i = 0; i < num_row; i++)
    {
        iperm[perm[i]] = i;
    }

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
    matrix = reordered_matrix;
    auto end11 = std::chrono::high_resolution_clock::now();
    std::cout << "Reorder in " << std::chrono::duration_cast<std::chrono::milliseconds>(end11 - end1).count() << " ms" << std::endl;
#else
    auto end11 = std::chrono::high_resolution_clock::now();
#endif
    //return;

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
        } /*
         if (__approx__nnz > MAX_NNZ)
         {
             std::ofstream logfile("result.log", std::ios::app);
             logfile << matrix_name << " " << num_lines << " cannot computer more than " << MAX_NNZ << "\n";
             logfile.close();
             std::cout << "More than " << MAX_NNZ << " non-zeros. Skipping..." << std::endl;
             exit(1);
         }*/
    }

    auto end3 = std::chrono::high_resolution_clock::now();

    std::cout << "Found fills in " << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - end2).count() << " ms" << std::endl;

    lower_mat_init(num_row, num_lines, m_rows, m_cols, m_values, matrix);

    matrix.clear();
    matrix.shrink_to_fit();

    lower_res_init(num_row, r_rows, r_cols, r_values, rev_graph);

    rev_graph.clear();
    rev_graph.shrink_to_fit();

    auto end41 = std::chrono::high_resolution_clock::now();

    std::cout << "Graph created in " << std::chrono::duration_cast<std::chrono::milliseconds>(end41 - end1).count() << " ms" << std::endl;

    compute_numeric_L(num_row, m_rows, m_cols, m_values, r_rows, r_cols, r_values);

    // std::cout << "Total iteration: " << total << std::endl;
    // total = 0;

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
        std::ofstream ofile("result/result.mtx");
        ofile << "%%MatrixMarket matrix coordinate real symmetric\n";
        ofile << num_row << " " << num_col << " " << r_rows[num_row] << "\n";

        for (mattype i = 0; i <= num_row - 1; i++)
        {
            for (indtype l = r_rows[i]; l < r_rows[i + 1]; l++)
            {
                ofile  << i + 1<< " " << r_cols[l] + 1  <<  (r_values[l] < 0 ? " " : "  ") << r_values[l] << "\n";
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

    // timer.stop(0);
    // timer.print_summary();
}

int main(int argc, char *argv[])
{
    // std::ofstream logfile("result.log", std::ios::app);
    std::cout <<
#ifdef REORDER
        "Reorder"
#else
        "Original"
#endif
              << " CuSparse Test" << std::endl;

    if (argc == 2)
    {
        std::cout << "\n\nStart for " << argv[1] << std::endl;
        operation_main(argv[1], true);
        return 0;
    }

    for (const auto &entry : std::filesystem::directory_iterator("matrices"))
    {
        std::cout << "\n\nStart for " << entry.path().c_str() << std::endl;
        operation_main(entry.path().c_str(), false);
    }

    return 0;
}