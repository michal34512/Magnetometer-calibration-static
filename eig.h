#ifndef EIG_H
#define EIG_H

#include "matrix.h"
#include "vector.h"
#include <stdint.h>

typedef enum {
    EIGEN_SUCCESS = 0,
    EIGEN_MAX_ITERATION,
    EIGEN_SINGULARITY,
} eigen_result_t;

typedef struct {
    mag_real_t eigenvalues[6];
    vector6_t eigenvectors[6];
    int8_t largest_eigenvalue_index;
} eigen6_t;

typedef struct {
    mag_real_t eigenvalues[3];
    vector3_t eigenvectors[3];
    int8_t largest_eigenvalue_index;
} eigen3_t;

eigen_result_t eig_solve_6(matrix6_t *matA, eigen6_t *eigen);
eigen_result_t eig_solve_3(matrix3_t *matA, eigen3_t *eigen);
void eig_orthogonalize_3(eigen3_t *eigen);

#endif  // EIG_H