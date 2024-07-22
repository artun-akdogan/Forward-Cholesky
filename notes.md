# Cholesky Formulas

$$ 
{\displaystyle L_{j,j}={\sqrt {A_{j,j}-\sum_{k=1}^{j-1}L^{2}_{j,k}}},} 
$$

$$ 
{\displaystyle L_{i,j}={\frac {1}{L_{j,j}}}\left(A_{i,j}-\sum_{k=1}^{j-1}L_{i,k}L_{j,k}\right)\quad {\text{for }}i>j.} 
$$

# METIS_NodeND
Let A be the original matrix and A′ be the permuted matrix. The arrays perm and iperm are defined as follows. Row (column) i of A′ is the perm(i) row (column) of A, and row (column) i of A is the iperm(i) row (column) of A′. the numbering of this vector starts from either 0 or 1, depending on the value of options(METIS_OPTION_NUMBERING).