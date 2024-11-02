#ifndef OPT_SEQUENTIAL_COMMON
#define OPT_SEQUENTIAL_COMMON

#include <climits>
#include <cmath>

#define COLMAX INT_MAX

#define MAX_NNZ 400000000

typedef int mattype;
typedef int indtype;
typedef double dattype;

struct sparse_raw
{
    dattype data;
    mattype row;
    mattype column;
};

long long total = 0;

inline indtype col_find(const mattype *arr, const mattype target, indtype start, indtype end)
{
    while (end > start)
    {
        total++;
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

inline indtype col_find_iter(const mattype *arr, const mattype target, indtype start, indtype end)
{
    while (end > start)
    {
        total++;
        if (arr[start] == target)
        {
            return start;
        }
        start++;
    }
    return COLMAX;
}

inline indtype col_find_custom(const mattype *arr, const mattype target, indtype start, indtype end)
{
    indtype iter = 1;
    while (end > start)
    {
        total++;
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

inline indtype col_find_custom_iter(const mattype *arr, const mattype target, indtype start, indtype end)
{
    indtype iter = 1;
    while (end > start)
    {
        total++;
        indtype mid = start + iter - 1;

        if (mid >= end)
        {
            iter = 1;
        }
        else if (arr[mid] > target)
        {
            end = mid;
            iter = 1;
        }
        else if (arr[mid] == target)
        {
            return mid;
        }
        else
        {
            start = mid + 1;
            iter <<= 1;
        }
    }
    return COLMAX;
}

#include <vector>
#include <stack>


// Function to perform Topological Sort
void topologicalSort(const std::vector<std::vector<mattype>> &tree,
                     std::vector<std::vector<mattype>> &topological,
                     mattype *topologicOrder)
{
    std::vector<bool> visited(tree.size(), false);
    std::stack<std::pair<int, int>> visit_stack;

    for (int i = 0; i < tree.size(); i++)
    {
        if (!visited[i]){
            visit_stack.push({i, 0});
        }
    }
    while (!visit_stack.empty())
    {
        std::pair<int, int>__recur_val = visit_stack.top();
        int v = __recur_val.first, d = __recur_val.second;
        visit_stack.pop();
        // Mark the current node as visited
        if(visited[v]){
            continue;
        }
        visited[v] = true;

        // Recur for all adjacent vertices
        for (auto x : tree[v])
        {
            if (!visited[x])
            {
                visit_stack.push({x, d + 1});
            }
        }

        // Push current vertex to stack which stores the result
        while (topological.size() <= d)
        {
            topological.push_back(std::vector<mattype>());
        }
        topological[d].push_back(v);
        topologicOrder[v] = d;
    }
}

#endif
