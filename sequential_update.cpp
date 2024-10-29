#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <cmath>
#include <algorithm>
#include <stack>
#include <chrono>
// #include <omp.h>
 #include <metis.h>

typedef unsigned long long ullong;

struct sparse
{
    double data;
    idx_t row;
    idx_t column;
};
/*
template<class InputIt1, class InputIt2>
bool intersect(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2)
{
    while (first1 != last1)
        if (std::binary_search(first2, last2, *first1++))
            return true;
    return false;
}
*/

const int tab32[32] = {
    0, 9, 1, 10, 13, 21, 2, 29,
    11, 14, 16, 18, 22, 25, 3, 30,
    8, 12, 20, 28, 15, 17, 24, 7,
    19, 27, 23, 6, 26, 5, 4, 31};

int log2_32(idx_t value)
{
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    return tab32[(idx_t)(value * 0x07C4ACDD) >> 27];
}

// Function to perform DFS and topological sorting
void topologicalSortUtil(int v, int d, const std::vector<std::vector<idx_t>> &adj,
                         std::vector<bool> &visited,
                         std::vector<std::vector<idx_t>> &Stack)
{
    // Mark the current node as visited
    visited[v] = true;

    // Recur for all adjacent vertices
    for (auto x : adj[v])
    {
        if (!visited[x])
            topologicalSortUtil(x, d + 1, adj, visited, Stack);
    }

    // Push current vertex to stack which stores the result
    while (Stack.size() <= d)
    {
        Stack.push_back(std::vector<idx_t>());
    }
    Stack[d].push_back(v);
}

// Function to perform Topological Sort
void topologicalSort(const std::vector<std::vector<idx_t>> &tree, std::vector<std::vector<idx_t>> &topological)
{
    std::vector<bool> visited(tree.size(), false);

    // Call the recursive helper function to store
    // Topological Sort starting from all vertices one by
    // one
    for (int i = tree.size() - 1; i >= 0; i--)
    {
        if (!visited[i])
            topologicalSortUtil(i, 0, tree, visited, topological);
    }
    /*
        for (int i = 0; i < topological.size(); i++){
            std::sort(topological[i].begin(), topological[i].end());
        }*/
}

idx_t calc_row(idx_t row,
               const std::vector<std::vector<idx_t>> &rev_graph,
               std::vector<idx_t> &m_rows,
               std::vector<idx_t> &m_cols,
               std::vector<double> &m_values,
               std::vector<idx_t> &r_rows,
               std::vector<idx_t> &r_cols,
               std::vector<double> &r_values)
{ /*
     int temp;
     if (row == 110)
     { // Segmentation
         std::cin >> temp;
     }*/
    // std::cout << row <<"--------------" << std::endl;
    idx_t tot_lines = 0;
    for (idx_t fi = r_rows[row]; fi < r_rows[row + 1] - 1; fi++)
    {
        double tot = 0, mat_val = 0, res_val = 0;
        // for (idx_t l = m_rows[i]; l < m_rows[i + 1]; l++)
        for (idx_t a = r_rows[r_cols[fi]], b = r_rows[row]; a < r_rows[r_cols[fi] + 1] && b < fi;)
        {
            // std::cout << "par: " << r_cols[fi] << " " << row << " " << r_cols[a] << " " << r_cols[b] << " " << tot << std::endl;
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

        // std::cout << "ok" << std::endl;

        auto mat_it = std::lower_bound(m_cols.begin() + m_rows[row], m_cols.begin() + m_rows[row + 1], r_cols[fi]) - m_cols.begin();
        if (m_cols[mat_it] == r_cols[fi])
        {
            mat_val = m_values[mat_it];
        }
        // std::cout << "ok" << std::endl;
        auto res_it = std::lower_bound(r_cols.begin() + r_rows[r_cols[fi]], r_cols.begin() + r_rows[r_cols[fi] + 1], r_cols[fi]) - r_cols.begin();
        if (r_cols[res_it] == r_cols[fi])
        {
            res_val = r_values[res_it];
        }
        if ((mat_val - tot) / res_val != 0)
        {
            r_values[fi] = (mat_val - tot) / res_val;
            tot_lines++;
        }
        // std::cout << mat_val << " " << tot << " " << res_val << std::endl;

        // std::cout << "!ok" << std::endl;
    }

    // std::cout << "!!!ok" << std::endl;

    double tot = 0, mat_val = 0;
    for (idx_t a = r_rows[row]; a < r_rows[row + 1] - 1; a++)
    {
        tot += r_values[a] * r_values[a];
    }
    auto mat_it = std::lower_bound(m_cols.begin() + m_rows[row], m_cols.begin() + m_rows[row + 1], row) - m_cols.begin();

    // std::cout << "!!!ok" << std::endl;
    if (m_cols[mat_it] == row)
    {
        mat_val = m_values[mat_it];
    }
    // std::cout << "!!!ok" << r_rows[row + 1] - 1 << " " << r_values.size() << std::endl;
    if (mat_val - tot != 0)
    {
        r_values[r_rows[row + 1] - 1] = std::sqrt(mat_val - tot);
        tot_lines++;
    }
    // std::cout << mat_val << " " << tot  << std::endl;
    // std::cout << "!!!ok" << std::endl;
    return tot_lines;
}

int main()
{

    auto beg = std::chrono::high_resolution_clock::now();

    std::ifstream file("ct20stif.mtx");
    int num_row, num_col, num_lines;

    // Ignore comments headers
    while (file.peek() == '%')
        file.ignore(2048, '\n');

    // Read number of rows and columns
    file >> num_row >> num_col >> num_lines;

    std::vector<std::vector<sparse>> matrix(num_row);
    std::vector<std::vector<sparse>> result(num_row);
    std::vector<std::set<idx_t>> graph(num_row);
    std::vector<std::vector<idx_t>> tree(num_row);
    std::vector<std::vector<idx_t>> topological;
    std::vector<std::vector<idx_t>> rev_graph(num_row);

    std::vector<idx_t> m_rows;
    std::vector<idx_t> m_cols;
    std::vector<double> m_values;

    std::vector<idx_t> r_rows;
    std::vector<idx_t> r_cols;
    std::vector<double> r_values;

    for (idx_t l = 0; l < num_lines; l++)
    {
        double data;
        idx_t row, col;
        file >> row >> col >> data;
        matrix[row - 1].push_back({data, row - 1, col - 1});
    }

    file.close();

    for (idx_t l = 0; l < num_row; l++)
    {
        std::sort(matrix[l].begin(), matrix[l].end(), [](const sparse &lhs, const sparse &rhs)
                  { return lhs.column < rhs.column; });
    }

    m_rows.push_back(0);
    for (idx_t l = 0; l < matrix.size(); l++)
    {
        for (idx_t i = 0; i < matrix[l].size(); i++)
        {
            m_cols.push_back(matrix[l][i].column);
            m_values.push_back(matrix[l][i].data);
        }
        m_rows.push_back(m_cols.size());
    }

    std::cout << "Total nonzero " << m_cols.size() << std::endl;
    auto end1 = std::chrono::high_resolution_clock::now();

    std::cout << "File read in " << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - beg).count() << " ms" << std::endl;

    auto end11 = std::chrono::high_resolution_clock::now();
    std::cout << "Reorder in " << std::chrono::duration_cast<std::chrono::milliseconds>(end11 - end1).count() << " ms" << std::endl;

    for (idx_t i = 0; i < m_rows.size() - 1; i++)
    {
        for (idx_t l = m_rows[i]; l < m_rows[i + 1]; l++)
        {
            if (m_cols[l] < i)
            {
                graph[m_cols[l]].insert(i);
                // rev_graph[matrix[i][l].column].insert(i);
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
    auto fit = graph[0].begin();

    auto end2 = std::chrono::high_resolution_clock::now();

    std::cout << "Base graph created in " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end1).count() << " ms" << std::endl;

    std::cout << "Finding fills" << std::endl;

    idx_t inserted = 0;

    for (idx_t i = 0; i < num_row; i++)
    {
        auto fit = graph[i].begin();

        if (fit == graph[i].end())
        {
            continue;
        }

        auto sit = fit;
        for (sit++; sit != graph[i].end(); sit++)
        {
            if (*fit == i || *sit == i)
            {
                std::cout << "wrong" << std::endl;
                return 0;
            }
            graph[*fit].insert(*sit);
        }
    }

    auto end3 = std::chrono::high_resolution_clock::now();

    std::cout << "Found fills in " << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - end2).count() << " ms" << std::endl;

    std::cout << "Creating elimination tree (transitive reduction)" << std::endl;

    for (idx_t i = 0; i < num_row; i++)
    {
        tree[*(graph[i].begin())].push_back(i);
    }
    topologicalSort(tree, topological);
    tree.clear();
    tree.shrink_to_fit();

    auto end31 = std::chrono::high_resolution_clock::now();

    std::cout << "Tree with topological depth " << topological.size() << " for rows " << num_row
              << " created in " << std::chrono::duration_cast<std::chrono::milliseconds>(end31 - end3).count() << "ms" << std::endl;

    std::cout << "Reversing Graph" << std::endl;
    idx_t tot = 0;

    for (idx_t i = 0; i < num_row; i++)
    {
        auto fit = graph[i].begin();

        for (; fit != graph[i].end(); fit++)
        {
            rev_graph[*fit].push_back(i);
            tot++;
        }
        tot++;
    }
    std::cout << "Total nonzero: " << tot << std::endl;

    auto end4 = std::chrono::high_resolution_clock::now();

    std::cout << "Reversed graph in " << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - end3).count() << " ms" << std::endl;

    std::cout << "Initializing result matrix" << std::endl;

    r_rows.push_back(0);
    for (idx_t l = 0; l < rev_graph.size(); l++)
    {
        for (idx_t i = 0; i < rev_graph[l].size(); i++)
        {
            r_cols.push_back(rev_graph[l][i]);
            // m_values.push_back(matrix[l][i].data);
        }
        r_cols.push_back(l);
        r_rows.push_back(r_cols.size());
    }
    r_values.resize(r_cols.size(), 0);

    std::cout << "R Mat " << r_rows[0] << r_rows[1] << std::endl;

    std::cout << "Graph created in " << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - end1).count() << " ms" << std::endl;

    graph.clear();
    graph.shrink_to_fit();
    auto end41 = std::chrono::high_resolution_clock::now();

    idx_t tot_lines = 0;

    for (idx_t i = 0; i < num_row; i++)
    {
        if (!(i % 1000))
        {
            std::cout << i << " " << num_row << std::endl;
        }
        tot_lines += calc_row(i, rev_graph, m_rows, m_cols, m_values, r_rows, r_cols, r_values);
    }

    auto end5 = std::chrono::high_resolution_clock::now();

    std::cout << "Operation completed for " << tot_lines << " nonzeros in " << std::chrono::duration_cast<std::chrono::milliseconds>(end5 - end41).count() << " ms" << std::endl;

    std::ofstream ofile("result.mtx");
    ofile << num_row << " " << num_col << " " << r_cols.size() << "\n";

    for (idx_t i = 0; i < r_rows.size() - 1; i++)
    {
        for (idx_t l = r_rows[i]; l < r_rows[i + 1]; l++)
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