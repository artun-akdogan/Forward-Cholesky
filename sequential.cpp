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

typedef unsigned int uint;
typedef unsigned long long ullong;

struct sparse
{
    double data;
    uint row;
    uint column;
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

int log2_32(uint32_t value)
{
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    return tab32[(uint32_t)(value * 0x07C4ACDD) >> 27];
}

// Function to perform DFS and topological sorting
void topologicalSortUtil(int v, int d, const std::vector<uint> &adj,
                         std::vector<bool> &visited,
                         std::vector<std::vector<uint>> &Stack)
{
    // Mark the current node as visited
    visited[v] = true;

    // Recur for all adjacent vertices

    if (!visited[adj[v]])
        topologicalSortUtil(adj[v], d + 1, adj, visited, Stack);

    // Push current vertex to stack which stores the result
    while (Stack.size() <= d)
    {
        Stack.push_back(std::vector<uint>());
    }
    Stack[d].push_back(v);
}

// Function to perform Topological Sort
void topologicalSort(const std::vector<uint> &tree, std::vector<std::vector<uint>> &topological)
{
    std::vector<bool> visited(tree.size(), false);

    // Call the recursive helper function to store
    // Topological Sort starting from all vertices one by
    // one
    for (int i = tree.size() - 1; i >=0; i--)
    {
        if (!visited[i])
            topologicalSortUtil(i, 0, tree, visited, topological);
    }
}

int main()
{

    auto beg = std::chrono::high_resolution_clock::now();

    std::ifstream file("bcsstk03.mtx");
    int num_row, num_col, num_lines, tot_lines = 0;

    // Ignore comments headers
    while (file.peek() == '%')
        file.ignore(2048, '\n');

    // Read number of rows and columns
    file >> num_row >> num_col >> num_lines;

    std::vector<std::vector<sparse>> matrix(num_row);
    std::vector<std::vector<sparse>> result(num_row);
    std::vector<std::set<uint>> graph(num_row);
    std::vector<uint> tree(num_row);
    std::vector<std::vector<uint>> topological;
    std::vector<std::vector<uint>> rev_graph(num_row);

    for (uint l = 0; l < num_lines; l++)
    {
        double data;
        uint row, col;
        file >> row >> col >> data;
        matrix[row - 1].push_back({data, row - 1, col - 1});
    }

    file.close();

    for (uint l = 0; l < num_row; l++)
    {
        std::sort(matrix[l].begin(), matrix[l].end(), [](const sparse &lhs, const sparse &rhs)
                  { return lhs.column < rhs.column; });
    }

    auto end1 = std::chrono::high_resolution_clock::now();

    std::cout << "File read in " << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - beg).count() << " ms" << std::endl;

    for (uint i = 0; i < num_row; i++)
    {
        for (uint l = 0; l < matrix[i].size(); l++)
        {
            if (matrix[i][l].column < i)
            {
                graph[matrix[i][l].column].insert(i);
                // rev_graph[matrix[i][l].column].insert(i);
            }
            else if (matrix[i][l].column > i)
            {
                std::cout << "wrong" << std::endl;
                return 0;
            }
        }
    }
    auto fit = graph[0].begin();

    auto end2 = std::chrono::high_resolution_clock::now();

    std::cout << "Base graph created in " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end1).count() << " ms" << std::endl;

    std::cout << "Finding fills" << std::endl;

    uint inserted = 0;

    for (uint i = 0; i < num_row; i++)
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

    for (uint i = 0; i < num_row; i++)
    {
        tree[*(graph[i].begin())] = i;
    }
    topologicalSort(tree, topological);
    tree.clear();
    tree.shrink_to_fit();

    auto end31 = std::chrono::high_resolution_clock::now();

    std::cout << "Tree with topological depth " << topological.size() << " for rows " << num_row
              << " created in " << std::chrono::duration_cast<std::chrono::milliseconds>(end31 - end3).count() << "ms" << std::endl;

    std::cout << "Reversing Graph" << std::endl;
    uint tot = 0;

    for (uint i = 0; i < num_row; i++)
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

    std::cout << "Graph created in " << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - end1).count() << " ms" << std::endl;

    graph.clear();
    graph.shrink_to_fit();

    for (uint i = 0; i < num_row; i++)
    {
        result[i].reserve(rev_graph[i].size() + 20);
        for (auto fit = rev_graph[i].begin(); fit != rev_graph[i].end(); fit++)
        {
            double tot = 0, mat_val = 0, res_val = 0;
            for (uint a = 0, b = 0; a < result[i].size() && b < result[*fit].size() && result[i][a].column < *fit && result[*fit][b].column < *fit;)
            {
                if (result[i][a].column < result[*fit][b].column)
                {
                    a++;
                }
                else if (result[i][a].column > result[*fit][b].column)
                {
                    b++;
                }
                else
                {
                    tot += result[i][a].data * result[*fit][b].data;
                    a++;
                    b++;
                }
            }

            auto mat_it = std::lower_bound(matrix[i].begin(), matrix[i].end(), (uint)*fit,
                                           [](const sparse &lhs, const uint rhs)
                                           { return lhs.column < rhs; });
            if (mat_it->column == *fit)
            {
                mat_val = mat_it->data;
            }
            auto res_it = std::lower_bound(result[*fit].begin(), result[*fit].end(), (uint)*fit,
                                           [](const sparse &lhs, const uint rhs)
                                           { return lhs.column < rhs; });
            if (res_it->column == *fit)
            {
                res_val = res_it->data;
            }
            if ((mat_val - tot) / res_val != 0)
            {
                result[i].push_back({(mat_val - tot) / res_val, i, *fit});
                tot_lines++;
            }
        }

        if (!(i % 1000))
        {
            std::cout << i << std::endl;
        }

        double tot = 0, mat_val = 0;
        for (uint a = 0; a < result[i].size(); a++)
        {
            tot += result[i][a].data * result[i][a].data;
        }
        auto mat_it = std::lower_bound(matrix[i].begin(), matrix[i].end(), i,
                                       [](const sparse &lhs, const uint rhs)
                                       { return lhs.column < rhs; });
        if (mat_it->column == i)
        {
            mat_val = mat_it->data;
        }
        if (mat_val - tot != 0)
        {
            result[i].push_back({std::sqrt(mat_val - tot), i, i});
            tot_lines++;
        }
    }

    auto end5 = std::chrono::high_resolution_clock::now();

    std::cout << "Operation completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(end5 - end4).count() << " ms" << std::endl;

    std::ofstream ofile("result.mtx");
    ofile << num_row << " " << num_col << " " << tot_lines << "\n";

    for (uint i = 0; i < num_row; i++)
    {
        for (uint l = 0; l < result[i].size(); l++)
        {
            ofile << i + 1 << " " << result[i][l].column + 1 << " " << result[i][l].data << "\n";
        }
    }

    ofile.close();

    auto end6 = std::chrono::high_resolution_clock::now();

    std::cout << "File written in " << std::chrono::duration_cast<std::chrono::milliseconds>(end6 - end5).count() << " ms" << std::endl;

    std::cout << "Program completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(end6 - beg).count() << " ms" << std::endl;

    return 0;
}