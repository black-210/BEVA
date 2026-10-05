#ifndef BEVA_QR_H
#define BEVA_QR_H
#include "matrix.h"
beva_status beva_qr_mgs(const beva_matrix*, beva_matrix*, beva_matrix*);
beva_status beva_qr_solve(const beva_matrix*, const beva_vector*, beva_vector*);
#endif
