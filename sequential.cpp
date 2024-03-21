#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <cmath>
#include <algorithm>

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
int main()
{
    std::ifstream file("ct20stif.mtx");
    int num_row, num_col, num_lines, tot_lines=0;

    // Ignore comments headers
    while (file.peek() == '%')
        file.ignore(2048, '\n');

    // Read number of rows and columns
    file >> num_row >> num_col >> num_lines;

    std::vector<std::vector<sparse>> matrix(num_row);
    std::vector<std::vector<sparse>> result(num_row);
    std::vector<std::set<uint>> graph(num_row);
    std::vector<std::set<uint>> rev_graph(num_row);

    for (uint l = 0; l < num_lines; l++)
    {
        double data;
        uint row, col;
        file >> row >> col >> data;
        matrix[row - 1].push_back({data, row - 1, col - 1});
    }

    file.close();

    std::cout << "File read" << std::endl;


    for (uint i = 0; i < num_row; i++)
    {
        for (uint l = 0; l < matrix[i].size(); l++)
        {
            if(matrix[i][l].column < i){
                graph[matrix[i][l].column].insert(i);
                //rev_graph[matrix[i][l].column].insert(i);
            }else if(matrix[i][l].column > i){
                std::cout << "wrong" << std::endl;
                return 0;
            }
        }
    }
    auto fit = graph[0].begin();

    for (; fit != graph[0].end(); fit++)
    {
        std::cout << *fit << " ";
    }
    std::cout << rev_graph[0].size() << " " << graph[0].size() << std::endl;
    std::cout << "Finding fills" << std::endl;

    uint inserted = 0;

    for (uint i = 0; i < num_row; i++)
    {
        std::cout <<inserted << " " << i << " " << num_row << " " << graph[i].size() << std::endl;
        auto fit = graph[i].begin();

        for (; fit != graph[i].end(); fit++)
        {
            auto sit = fit;
            for(sit++; sit != graph[i].end(); sit++){
                if(*fit==i || *sit==i){
                std::cout << "wrong" << std::endl;
                return 0;
                }
                graph[*fit].insert(*sit);
            }
        }
    }
    std::cout << rev_graph[0].size() << " " << graph[0].size() << std::endl;
    std::cout << "Reversing Graph" << std::endl;

    for (uint i = 0; i < num_row; i++)
    {
        auto fit = graph[i].begin();

        for (; fit != graph[i].end(); fit++)
        {
            rev_graph[*fit].insert(i);
        }
    }
    std::cout << rev_graph[0].size() << " " << graph[0].size() << std::endl;

    std::cout << "Graph created" << std::endl;

    for (uint i = 0; i < num_row; i++)
    {
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
        if (std::sqrt(mat_val - tot) != 0)
        {
            result[i].push_back({std::sqrt(mat_val - tot), i, i});
            tot_lines++;
        }
    }

    std::cout << "Operation completed" << std::endl;

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

    std::cout << "File written" << std::endl;

    return 0;
}