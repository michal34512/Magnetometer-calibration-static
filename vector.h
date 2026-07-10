#ifndef VECTOR_H
#define VECTOR_H

#include "magcalib_defs.h"

typedef struct {
    mag_real_t v[3];
} vector3_t;

typedef struct {
    mag_real_t v[6];
} vector6_t;

void vec_rotate_xyz(vector3_t *vec, mag_real_t angle_x, mag_real_t angle_y, mag_real_t angle_z);

void vec_mul_scalar(vector3_t *vec, mag_real_t scalar);

void vec_normalize(vector3_t *vec);

void vec_ortogonalize_onto(vector3_t *vec_a, vector3_t *vec_b);

void vec_print_3(vector3_t *vec_a);

#endif  // VECTOR_H