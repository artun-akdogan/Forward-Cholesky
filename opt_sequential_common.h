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

// long long total = 0;

class Timer
{
private:
    using Clock = std::chrono::high_resolution_clock;
    std::vector<std::chrono::duration<double>> durations; // Stores durations for each block
    std::vector<Clock::time_point> start_times;           // Stores start times for each block
    int num_blocks;                                       // Total number of blocks

public:
    Timer(int num_blocks) : num_blocks(num_blocks)
    {
        durations.resize(num_blocks, std::chrono::duration<double>::zero());
        start_times.resize(num_blocks);
    }

    // Start timing for a specific block (index-based)
    inline void start(int block_index)
    {
        start_times[block_index] = Clock::now();
    }

    // Stop timing for a specific block and accumulate its duration
    inline void stop(int block_index)
    {
        durations[block_index] += Clock::now() - start_times[block_index];
    }

    // Print all timings as percentages of total (excluding index 0)
    void print_summary() const
    {

        std::cout << "Timing Summary:\n";
        for (int i = 1; i < num_blocks; ++i)
        { // Skip index 0
            double percentage = (durations[0].count() > 0.0) ? (durations[i].count() / durations[0].count() * 100.0) : 0.0;
            std::cout << "  Block " << i << ": " << durations[i].count() << " seconds ("
                      << percentage << "%)\n";
        }

        // Include "Total" as a separate summary
        std::cout << "  Total: " << durations[0].count() << " seconds (100.00%)\n";
    }
} timer(5);

inline indtype col_find(const mattype *arr, const mattype target, indtype start, indtype end)
{
    while (end > start)
    {
        // total++;
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
        // total++;
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
        // total++;
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
        // total++;
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
        if (!visited[i])
        {
            visit_stack.push({i, 0});
        }
    }
    while (!visit_stack.empty())
    {
        std::pair<int, int> __recur_val = visit_stack.top();
        int v = __recur_val.first, d = __recur_val.second;
        visit_stack.pop();
        // Mark the current node as visited
        if (visited[v])
        {
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
