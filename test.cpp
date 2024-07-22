#include <iostream>
#include <metis.h>

int main() {
    // Example symmetric matrix in CSR format (lower triangle)
    idx_t n = 4; // Number of vertices
    idx_t nnz = 6; // Number of non-zeros in lower triangle
    idx_t *xadj = new idx_t[n + 1]{0, 2, 4, 6, 6}; // Pointers to the start of each vertex's adjacency list
    idx_t *adjncy = new idx_t[nnz]{1, 2, 2, 3, 3, 3}; // Adjacency list of each vertex

    // METIS options
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    // METIS output
    idx_t *perm = new idx_t[n]; // Permutation array
    idx_t *iperm = new idx_t[n]; // Inverse permutation array

    // Perform fill-reducing matrix ordering
    METIS_NodeND(&n, xadj, adjncy, nullptr, options, perm, iperm);

    // Print the permutation
    std::cout << "Permutation: ";
    for (idx_t i = 0; i < n; ++i) {
        std::cout << perm[i] << " ";
    }
    std::cout << std::endl;

    // Clean up
    delete[] xadj;
    delete[] adjncy;
    delete[] perm;
    delete[] iperm;

    return 0;
}