#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "eig.h"
#include "magcalib.h"
#include "matrix.h"

static inline void mat_mul_6(matrix6_t *mat, matrix6_t *q, matrix6_t *r) {
    for (uint8_t i = 0; i < 6; i++) {
        for (uint8_t j = 0; j < 6; j++) {
            mat->m[j][i] = 0.0;
            for (uint8_t k = 0; k < 6; k++) {
                mat->m[j][i] += q->m[j][k] * r->m[k][i];
            }
        }
    }
}

static inline bool eig_qr_decomposition(const matrix6_t *mat, matrix6_t *q, matrix6_t *r) {
    memset(q, 0, sizeof(matrix6_t));
    memset(r, 0, sizeof(matrix6_t));

    mag_real_t v[6];
    for (uint8_t i = 0; i < 6; i++) {
        for (uint8_t k = 0; k < 6; k++) {
            v[k] = mat->m[k][i];
        }

        for (uint8_t j = 0; j < i; j++) {
            mag_real_t currentDotProduct = 0.0f;
            for (uint8_t k = 0; k < 6; k++) {
                currentDotProduct += v[k] * q->m[k][j];
            }
            for (uint8_t k = 0; k < 6; k++) {
                v[k] -= q->m[k][j] * currentDotProduct;
            }
            r->m[j][i] = currentDotProduct;
        }
        mag_real_t norm = 0.0f;
        for (uint8_t k = 0; k < 6; k++) {
            norm += v[k] * v[k];
        }
        norm = sqrt(norm);
        if (fabs(norm) < EIGEN_SINGULARITY_TOLERANCE) {
            return false;
        }
        for (uint8_t k = 0; k < 6; k++) {
            q->m[k][i] = v[k] / norm;
        }
        r->m[i][i] = norm;
    }

    return true;
}

static inline bool mat_is_upper_triang_6(matrix6_t *mat, mag_real_t tolerance) {
    for (uint8_t i = 1; i < 6; i++) {
        for (uint8_t j = 0; j < i; j++) {
            if (fabs(mat->m[i][j]) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

static eigen_result_t eig_solve_eigenvalues_6(matrix6_t *matA, eigen6_t *eigen) {
    matrix6_t X;
    memcpy(&X, matA, sizeof(matrix6_t));
    uint16_t i = 0;
    do {
        matrix6_t Q, R;
        if (!eig_qr_decomposition(&X, &Q, &R)) {
            return EIGEN_SINGULARITY;
        }
        mat_mul_6(&X, &R, &Q);
        i++;
    } while (!mat_is_upper_triang_6(&X, EIGEN_TOLERANCE) && (i < EIGEN_MAX_ITER));

    for (uint8_t j = 0; j < 6; j++) {
        eigen->eigenvalues[j] = X.m[j][j];
        if (eigen->largest_eigenvalue_index == -1
            || eigen->eigenvalues[j] > eigen->eigenvalues[eigen->largest_eigenvalue_index]) {
            eigen->largest_eigenvalue_index = j;
        }
    }

    return i == EIGEN_MAX_ITER ? EIGEN_MAX_ITERATION : EIGEN_SUCCESS;
}

static inline void eig_linsolve_upper_triang_6(matrix6_t *A, vector6_t *b, vector6_t *x) {
    assert(A != NULL);
    assert(b != NULL);
    assert(x != NULL);

    mag_real_t back_substitute;
    for (int8_t i = 5; i >= 0; i--) {
        back_substitute = 0;
        for (uint8_t j = i + 1; j < 6; j++) {
            back_substitute += x->v[j] * A->m[i][j];
        }
        x->v[i] = (b->v[i] - back_substitute) / A->m[i][i];
    }
}

static bool eig_linsolve_qr_6(matrix6_t *A, vector6_t *b, vector6_t *x) {
    matrix6_t Q, R;
    vector6_t QTb;
    if (!eig_qr_decomposition(A, &Q, &R)) {
        return false;
    }
    for (uint8_t r = 0; r < 6; r++) {
        QTb.v[r] = 0.0f;
        for (uint8_t k = 0; k < 6; k++) {
            QTb.v[r] += Q.m[k][r] * b->v[k];
        }
    }
    eig_linsolve_upper_triang_6(&R, &QTb, x);
    return true;
}

static inline bool vec_equal_6(vector6_t *vecA, vector6_t *vecB, double tol) {
    for (int i = 0; i < 6; i++) {
        if (fabs(vecA->v[i] - vecB->v[i]) > tol)
            return false;
    }
    return true;
}

static eigen_result_t eig_solve_eigenvectors_6(matrix6_t *matA, eigen6_t *eigen) {
    for (uint8_t i = 0; i < 6; i++) {
        double lambda = eigen->eigenvalues[i] + ((double)rand() / (double)RAND_MAX) * 0.000001;

        // A - lambda*I
        matrix6_t A_minus_lambda_I;
        for (uint8_t r = 0; r < 6; r++) {
            for (uint8_t c = 0; c < 6; c++) {
                A_minus_lambda_I.m[r][c] = matA->m[r][c];
                if (r == c) {
                    A_minus_lambda_I.m[r][c] -= lambda;
                }
            }
        }

        int j = 0;
        vector6_t previous;
        for (size_t j = 0; j < 6; j++) {
            eigen->eigenvectors[i].v[j] = 1.0;
        }
        do {
            previous = eigen->eigenvectors[i];
            if (!eig_linsolve_qr_6(&A_minus_lambda_I, &previous, &eigen->eigenvectors[i])) {
                return EIGEN_SINGULARITY;
            }

            if (eigen->eigenvectors[i].v[0] < 0) {
                for (size_t j = 0; j < 6; j++) {
                    eigen->eigenvectors[i].v[j] = -eigen->eigenvectors[i].v[j];
                }
            }
            mag_real_t norm = 0.0;
            for (uint8_t j = 0; j < 6; j++) {
                norm += eigen->eigenvectors[i].v[j] * eigen->eigenvectors[i].v[j];
            }
            norm = sqrt(norm);
            if (fabs(norm) < EIGEN_SINGULARITY_TOLERANCE) {
                return EIGEN_SINGULARITY;
            }
            for (uint8_t j = 0; j < 6; j++) {
                eigen->eigenvectors[i].v[j] /= norm;
            }
            j++;
        } while (!vec_equal_6(&eigen->eigenvectors[i], &previous, EIGEN_TOLERANCE) && (j < EIGEN_MAX_ITER));

        if (j == EIGEN_MAX_ITER) {
            return EIGEN_MAX_ITERATION;
        }
    }

    return EIGEN_SUCCESS;
}

eigen_result_t eig_solve_6(matrix6_t *matA, eigen6_t *eigen) {
    assert(matA != NULL);
    assert(eigen != NULL);
    eigen_result_t res;
    eigen->largest_eigenvalue_index = -1;

    res = eig_solve_eigenvalues_6(matA, eigen);

    if (res != EIGEN_SUCCESS && res != EIGEN_MAX_ITERATION) {
        return res;
    }

    res = eig_solve_eigenvectors_6(matA, eigen);
    return res;
}

static inline void sort3_desc(mag_real_t *a, mag_real_t *b, mag_real_t *c) {
    mag_real_t t;
    if (*a < *b) {
        t  = *a;
        *a = *b;
        *b = t;
    }
    if (*a < *c) {
        t  = *a;
        *a = *c;
        *c = t;
    }
    if (*b < *c) {
        t  = *b;
        *b = *c;
        *c = t;
    }
}

void eig_fast_eigenvalues_3(matrix3_t *matA, eigen3_t *eigen) {
    const mag_real_t a00 = matA->m[0][0];
    const mag_real_t a01 = matA->m[0][1];
    const mag_real_t a02 = matA->m[0][2];
    const mag_real_t a11 = matA->m[1][1];
    const mag_real_t a12 = matA->m[1][2];
    const mag_real_t a22 = matA->m[2][2];
    mag_real_t p1        = a01 * a01 + a02 * a02 + a12 * a12;

    if (p1 == 0.0) {
        eigen->eigenvalues[0] = a00;
        eigen->eigenvalues[1] = a11;
        eigen->eigenvalues[2] = a22;
        sort3_desc(&eigen->eigenvalues[0], &eigen->eigenvalues[1], &eigen->eigenvalues[2]);
        return;
    }
    const mag_real_t q   = (a00 + a11 + a22) / 3.0;
    const mag_real_t b00 = a00 - q;
    const mag_real_t b11 = a11 - q;
    const mag_real_t b22 = a22 - q;
    const mag_real_t p2  = b00 * b00 + b11 * b11 + b22 * b22 + 2.0 * p1;
    const mag_real_t p   = sqrt(p2 / 6.0);
    const mag_real_t detB
      = b00 * (b11 * b22 - a12 * a12) - a01 * (a01 * b22 - a02 * a12) + a02 * (a01 * a12 - b11 * a02);
    mag_real_t r = detB / (2.0 * p * p * p);
    if (r <= -1.0)
        r = -1.0;
    if (r >= 1.0)
        r = 1.0;
    const mag_real_t phi = acos(r) / 3.0;
    mag_real_t eig1      = q + 2.0 * p * cos(phi);
    mag_real_t eig2      = q + 2.0 * p * cos(phi + (2.0 * M_PI / 3.0));
    mag_real_t eig3      = q + 2.0 * p * cos(phi + (4.0 * M_PI / 3.0));

    sort3_desc(&eig1, &eig2, &eig3);
    eigen->eigenvalues[0] = eig1;
    eigen->eigenvalues[1] = eig2;
    eigen->eigenvalues[2] = eig3;
}

void eig_fast_eigenvectors_3(matrix3_t *matA, eigen3_t *eigen) {
    for (uint8_t i = 0; i < 3; i++) {
        const mag_real_t lambda = eigen->eigenvalues[i];
        if (eigen->largest_eigenvalue_index == -1 || lambda > eigen->eigenvalues[eigen->largest_eigenvalue_index]) {
            eigen->largest_eigenvalue_index = i;
        }
        const mag_real_t r0[3] = {matA->m[0][0] - lambda, matA->m[0][1], matA->m[0][2]};
        const mag_real_t r1[3] = {matA->m[1][0], matA->m[1][1] - lambda, matA->m[1][2]};
        const mag_real_t r2[3] = {matA->m[2][0], matA->m[2][1], matA->m[2][2] - lambda};

        const mag_real_t v0[3]
          = {r0[1] * r1[2] - r0[2] * r1[1], r0[2] * r1[0] - r0[0] * r1[2], r0[0] * r1[1] - r0[1] * r1[0]};  // r0 x r1
        const mag_real_t v1[3]
          = {r1[1] * r2[2] - r1[2] * r2[1], r1[2] * r2[0] - r1[0] * r2[2], r1[0] * r2[1] - r1[1] * r2[0]};  // r1 x r2
        const mag_real_t v2[3]
          = {r2[1] * r0[2] - r2[2] * r0[1], r2[2] * r0[0] - r2[0] * r0[2], r2[0] * r0[1] - r2[1] * r0[0]};  // r2 x r0

        const mag_real_t mag0 = v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2];
        const mag_real_t mag1 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
        const mag_real_t mag2 = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];

        const mag_real_t *best_v = v0;
        mag_real_t max_mag       = mag0;

        if (mag1 > max_mag) {
            best_v  = v1;
            max_mag = mag1;
        }
        if (mag2 > max_mag) {
            best_v  = v2;
            max_mag = mag2;
        }

        mag_real_t inv_mag = 1.0 / sqrt(max_mag);
        if (best_v[0] > 0) {
            eigen->eigenvectors[i].v[0] = best_v[0] * inv_mag;
            eigen->eigenvectors[i].v[1] = best_v[1] * inv_mag;
            eigen->eigenvectors[i].v[2] = best_v[2] * inv_mag;
        } else {
            eigen->eigenvectors[i].v[0] = -best_v[0] * inv_mag;
            eigen->eigenvectors[i].v[1] = -best_v[1] * inv_mag;
            eigen->eigenvectors[i].v[2] = -best_v[2] * inv_mag;
        }
    }
}

eigen_result_t eig_solve_3(matrix3_t *matA, eigen3_t *eigen) {
    assert(matA != NULL);
    assert(eigen != NULL);
    eigen->largest_eigenvalue_index = -1;
    eig_fast_eigenvalues_3(matA, eigen);
    eig_fast_eigenvectors_3(matA, eigen);
    return EIGEN_SUCCESS;
}

void eig_orthogonalize_3(eigen3_t *eigen) {
    assert(eigen != NULL);
    // Modified Gram–Schmidt (MGS) orthogonalize eigenvectors for better results
    // For best result start with vector with largest eigenvalue
    switch (eigen->largest_eigenvalue_index) {
    case 0:
        vec_ortogonalize_onto(&eigen->eigenvectors[1], &eigen->eigenvectors[0]);
        vec_ortogonalize_onto(&eigen->eigenvectors[2], &eigen->eigenvectors[0]);
        vec_ortogonalize_onto(&eigen->eigenvectors[1], &eigen->eigenvectors[2]);
        break;
    case 1:
        vec_ortogonalize_onto(&eigen->eigenvectors[0], &eigen->eigenvectors[1]);
        vec_ortogonalize_onto(&eigen->eigenvectors[2], &eigen->eigenvectors[1]);
        vec_ortogonalize_onto(&eigen->eigenvectors[0], &eigen->eigenvectors[2]);
        break;
    case 2:
        vec_ortogonalize_onto(&eigen->eigenvectors[0], &eigen->eigenvectors[2]);
        vec_ortogonalize_onto(&eigen->eigenvectors[1], &eigen->eigenvectors[2]);
        vec_ortogonalize_onto(&eigen->eigenvectors[0], &eigen->eigenvectors[1]);
        break;
    }
    // eigenvector 0 does not need normalization as it is already normalized
    vec_normalize(&eigen->eigenvectors[1]);
    vec_normalize(&eigen->eigenvectors[2]);
}
