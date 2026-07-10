#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "vector.h"

void vec_rotate_x(vector3_t *vec, mag_real_t angle) {
    assert(vec != NULL);
    mag_real_t y = vec->v[1];
    mag_real_t z = vec->v[2];
    vec->v[1]    = y * cos(angle) - z * sin(angle);
    vec->v[2]    = y * sin(angle) + z * cos(angle);
}

void vec_rotate_y(vector3_t *vec, mag_real_t angle) {
    assert(vec != NULL);
    mag_real_t x = vec->v[0];
    mag_real_t z = vec->v[2];
    vec->v[0]    = x * cos(angle) + z * sin(angle);
    vec->v[2]    = -x * sin(angle) + z * cos(angle);
}

void vec_rotate_z(vector3_t *vec, mag_real_t angle) {
    assert(vec != NULL);
    mag_real_t x = vec->v[0];
    mag_real_t y = vec->v[1];
    vec->v[0]    = x * cos(angle) - y * sin(angle);
    vec->v[1]    = x * sin(angle) + y * cos(angle);
}

void vec_rotate_xyz(vector3_t *vec, mag_real_t angle_x, mag_real_t angle_y, mag_real_t angle_z) {
    assert(vec != NULL);
    vec_rotate_x(vec, angle_x);
    vec_rotate_y(vec, angle_y);
    vec_rotate_z(vec, angle_z);
}

void vec_mul_scalar(vector3_t *vec, mag_real_t scalar) {
    assert(vec != NULL);
    vec->v[0] *= scalar;
    vec->v[1] *= scalar;
    vec->v[2] *= scalar;
}

void vec_normalize(vector3_t *vec) {
    assert(vec != NULL);
    mag_real_t length = sqrt(vec->v[0] * vec->v[0] + vec->v[1] * vec->v[1] + vec->v[2] * vec->v[2]);
    if (length > 0) {
        vec->v[0] /= length;
        vec->v[1] /= length;
        vec->v[2] /= length;
    }
}

void vec_ortogonalize_onto(vector3_t *vec_a, vector3_t *vec_b) {
    assert(vec_a != NULL);
    assert(vec_b != NULL);

    mag_real_t vec_b_len_sq = vec_b->v[0] * vec_b->v[0] + vec_b->v[1] * vec_b->v[1] + vec_b->v[2] * vec_b->v[2];
    if (vec_b_len_sq == 0) {
        vec_b_len_sq = 1;
    }
    mag_real_t dot = vec_a->v[0] * vec_b->v[0] + vec_a->v[1] * vec_b->v[1] + vec_a->v[2] * vec_b->v[2];
    vec_a->v[0]    = vec_a->v[0] - (dot / vec_b_len_sq) * vec_b->v[0];
    vec_a->v[1]    = vec_a->v[1] - (dot / vec_b_len_sq) * vec_b->v[1];
    vec_a->v[2]    = vec_a->v[2] - (dot / vec_b_len_sq) * vec_b->v[2];
}

void vec_print_3(vector3_t *vec_a) {
    assert(vec_a != NULL);
    printf("Vector 3: [%.6f, %.6f, %.6f]\n", vec_a->v[0], vec_a->v[1], vec_a->v[2]);
}