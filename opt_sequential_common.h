#ifndef OPT_SEQUENTIAL_COMMON
#define OPT_SEQUENTIAL_COMMON

#include <climits>
#include <cmath>

#define COLMAX INT_MAX

typedef int mattype;
typedef int indtype;
typedef double dattype;


struct sparse_raw
{
    dattype data;
    mattype row;
    mattype column;
};


inline indtype col_find(const mattype *arr, mattype target, indtype start, indtype end)
{
    while (end > start)
    {
        indtype mid = (start + end) / 2;

        if (arr[mid] > target)
        {
            end = mid;
        }
        else if (arr[mid] < target)
        {
            start = mid + 1;
        }
        else
        {
            return mid;
        }
    }
    return COLMAX;
}

inline indtype col_find_it(const mattype *arr, mattype target, indtype start, indtype end)
{
    indtype iter = 1;
    while (end > start)
    {
        indtype mid = start + iter - 1;

        if (mid >= end)
        {
            return col_find(arr, target, start, end);
        }
        if (arr[mid] > target)
        {
            return col_find(arr, target, start, mid);
        }
        if (arr[mid] == target)
        {
            return mid;
        }
        start = mid + 1;
        iter <<= 1;
    }
    return COLMAX;
}

#include <vector>

// Function to perform DFS and topological sorting
void topologicalSortUtil(int v, int d, const std::vector<std::vector<mattype>> &adj,
                         std::vector<bool> &visited,
                         std::vector<std::vector<mattype>> &Stack)
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
        Stack.push_back(std::vector<mattype>());
    }
    Stack[d].push_back(v);
}

// Function to perform Topological Sort
void topologicalSort(const std::vector<std::vector<mattype>> &tree, std::vector<std::vector<mattype>> &topological)
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

#endif
