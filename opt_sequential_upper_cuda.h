#ifndef OPT_SEQUENTIAL_UPPER
#define OPT_SEQUENTIAL_UPPER

#include "opt_sequential_common.h"

void upper_cholesky_calculate(const mattype num_rows,
                              const indtype *m_rows,
                              const mattype *m_cols,
                              const dattype *m_values,
                              const indtype *r_rows,
                              const mattype *r_cols,
                              dattype *r_values);

#endif
