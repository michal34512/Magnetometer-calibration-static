#ifndef MATRIX_H
#define MATRIX_H

#include "magcalib.h"
#include <stdbool.h>

typedef struct {
    mag_real_t m[3][3];
} matrix3_t;

typedef struct {
    mag_real_t m[6][6];
} matrix6_t;

typedef struct {
    mag_real_t m[6][4];
} matrix6x4_t;

typedef struct {
    mag_real_t m[4][4];
} matrix4_t;

typedef struct {
    mag_real_t m[10][10];
} matrix10_t;

bool mat_inv_4x4(matrix4_t *mat, matrix4_t *mat_out);

#endif  // MATRIX_H