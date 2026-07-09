#include <stdint.h>
#include <stdio.h>

#include "magcalib_defs.h"
#include "matrix.h"

bool mat_inv_4x4(matrix4_t *mat, matrix4_t *mat_out) {
    mag_real_t inv[16];
    mag_real_t m[16];

    for (uint8_t r = 0; r < 4; r++) {
        for (uint8_t c = 0; c < 4; c++) {
            m[r * 4 + c] = mat->m[r][c];
        }
    }

    inv[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14]
             + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];

    inv[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14]
             - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];

    inv[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13]
             + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];

    inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13]
              - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];

    inv[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14]
             - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];

    inv[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14]
             + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];

    inv[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13]
             - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];

    inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13]
              + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];

    inv[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] + m[13] * m[2] * m[7]
             - m[13] * m[3] * m[6];

    inv[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14]
             - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];

    inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13]
              + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];

    inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13]
              - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] - m[9] * m[2] * m[7]
             + m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] + m[8] * m[2] * m[7]
             - m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[9] - m[8] * m[1] * m[7]
              + m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[9] + m[8] * m[1] * m[6]
              - m[8] * m[2] * m[5];

    mag_real_t det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
    if (det == 0.0) {
        return false;
    }

    det = 1.0 / det;
    for (uint8_t r = 0; r < 4; r++) {
        for (uint8_t c = 0; c < 4; c++) {
            mat_out->m[r][c] = (mag_real_t)(inv[r * 4 + c] * det);
        }
    }

    return true;
}

void mat_print_3(matrix3_t *mat) {
    printf("---- 3x3 Matrix ---\n");
    for (uint8_t r = 0; r < 3; r++) {
        for (uint8_t c = 0; c < 3; c++) {
            printf("%.3f ", mat->m[r][c]);
        }
        printf("\n");
    }
    printf("-------------------\n");
}

void mat_print_4(matrix4_t *mat) {
    printf("---- 4x4 Matrix ---\n");
    for (uint8_t r = 0; r < 4; r++) {
        for (uint8_t c = 0; c < 4; c++) {
            printf("%.3f ", mat->m[r][c]);
        }
        printf("\n");
    }
    printf("-------------------\n");
}

void mat_print_6x4(matrix6x4_t *mat) {
    printf("---- 6x4 Matrix ---\n");
    for (uint8_t r = 0; r < 6; r++) {
        for (uint8_t c = 0; c < 4; c++) {
            printf("%.3f ", mat->m[r][c]);
        }
        printf("\n");
    }
    printf("-------------------\n");
}

void mat_print_6(matrix6_t *mat) {
    printf("---- 6x6 Matrix ---\n");
    for (uint8_t r = 0; r < 6; r++) {
        for (uint8_t c = 0; c < 6; c++) {
            printf("%.3f ", mat->m[r][c]);
        }
        printf("\n");
    }
    printf("-------------------\n");
}