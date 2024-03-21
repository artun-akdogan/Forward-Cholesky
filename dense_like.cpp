#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

typedef unsigned int uint;

struct sparse
{
    double data;
    uint row;
    uint column;
};

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
        for (uint l = 0; l < i; l++)
        {
            double tot = 0, mat_val = 0, res_val = 0;
            for (uint a = 0, b = 0; a < result[i].size() && b < result[l].size() && result[i][a].column < l && result[l][b].column < l;)
            {

                if (result[i][a].column < result[l][b].column)
                {
                    a++;
                }
                else if (result[i][a].column > result[l][b].column)
                {
                    b++;
                }
                else
                {
                    tot += result[i][a].data * result[l][b].data;
                    a++;
                    b++;
                }
            }
            auto mat_it = std::lower_bound(matrix[i].begin(), matrix[i].end(), l,
                                           [](const sparse &lhs, const uint rhs)
                                           { return lhs.column < rhs; });
            if (mat_it->column == l)
            {
                mat_val = mat_it->data;
            }
            auto res_it = std::lower_bound(result[l].begin(), result[l].end(), l,
                                           [](const sparse &lhs, const uint rhs)
                                           { return lhs.column < rhs; });
            if (res_it->column == l)
            {
                res_val = res_it->data;
            }
            if ((mat_val - tot) / res_val != 0)
            {
                result[i].push_back({(mat_val - tot) / res_val, i, l});
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