#ifndef VECTOR_H
#define VECTOR_H

#include "magcalib.h"

typedef struct {
    mag_real_t x;
    mag_real_t y;
    mag_real_t z;
} vector3_t;

typedef struct {
    mag_real_t v[6];
} vector6_t;

void vec_rotate_xyz(vector3_t *vec, mag_real_t angle_x, mag_real_t angle_y, mag_real_t angle_z);

void vec_mul_scalar(vector3_t *vec, mag_real_t scalar);

void vec_normalize(vector3_t *vec);

#endif  // VECTOR_H