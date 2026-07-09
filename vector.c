#include "vector.h"

void vec_rotate_x(vector3_t *vec, mag_real_t angle) {
    mag_real_t y = vec->y;
    mag_real_t z = vec->z;
    vec->y       = y * cos(angle) - z * sin(angle);
    vec->z       = y * sin(angle) + z * cos(angle);
}

void vec_rotate_y(vector3_t *vec, mag_real_t angle) {
    mag_real_t x = vec->x;
    mag_real_t z = vec->z;
    vec->x       = x * cos(angle) + z * sin(angle);
    vec->z       = -x * sin(angle) + z * cos(angle);
}

void vec_rotate_z(vector3_t *vec, mag_real_t angle) {
    mag_real_t x = vec->x;
    mag_real_t y = vec->y;
    vec->x       = x * cos(angle) - y * sin(angle);
    vec->y       = x * sin(angle) + y * cos(angle);
}

void vec_rotate_xyz(vector3_t *vec, mag_real_t angle_x, mag_real_t angle_y, mag_real_t angle_z) {
    vec_rotate_x(vec, angle_x);
    vec_rotate_y(vec, angle_y);
    vec_rotate_z(vec, angle_z);
}

void vec_mul_scalar(vector3_t *vec, mag_real_t scalar) {
    vec->x *= scalar;
    vec->y *= scalar;
    vec->z *= scalar;
}

void vec_normalize(vector3_t *vec) {
    mag_real_t length = sqrt(vec->x * vec->x + vec->y * vec->y + vec->z * vec->z);
    if (length > 0) {
        vec->x /= length;
        vec->y /= length;
        vec->z /= length;
    }
}