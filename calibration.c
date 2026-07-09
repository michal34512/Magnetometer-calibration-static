#include "calibration.h"
#include "eig.h"
#include "magcalib.h"
#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#define ALPHA 16

matrix6_t c_inv = {
    .m = {
        {(mag_real_t)(ALPHA - 4) / (-ALPHA * ALPHA + 3 * ALPHA), -(mag_real_t)(ALPHA - 2) / (-ALPHA * ALPHA + 3 * ALPHA), -(mag_real_t)(ALPHA - 2) / (-ALPHA * ALPHA + 3 * ALPHA), 0.f, 0.f, 0.f},
        {-(mag_real_t)(ALPHA - 2) / (-ALPHA * ALPHA + 3 * ALPHA), (mag_real_t)(ALPHA - 4) / (-ALPHA * ALPHA + 3 * ALPHA), -(mag_real_t)(ALPHA - 2) / (-ALPHA * ALPHA + 3 * ALPHA), 0.f, 0.f, 0.f},
        {-(mag_real_t)(ALPHA - 2) / (-ALPHA * ALPHA + 3 * ALPHA), -(mag_real_t)(ALPHA - 2) / (-ALPHA * ALPHA + 3 * ALPHA), (mag_real_t)(ALPHA - 4) / (-ALPHA * ALPHA + 3 * ALPHA)   , 0.f, 0.f, 0.f},
        {0.f, 0.f, 0.f, -1.f / ALPHA, 0.f, 0.f},
        {0.f, 0.f, 0.f, 0.f, -1.f / ALPHA, 0.f},
        {0.f, 0.f, 0.f, 0.f, 0.f, -1.f / ALPHA},
    },
};

calibration_result_t calib_fit_ellipsoid(calibration_result_data_t *calib_res, calibration_object_t *calib_obj) {
    assert(calib_res != NULL);
    assert(calib_obj != NULL);

    // M = C1^-1*(S11-S12*(S22^-1*S21));
    calibration_result_t res;
    matrix4_t s22_inv;
    matrix6x4_t s22_inv_s21;
    matrix6_t s11_s12_s22_inv_s21;
    matrix6_t m;

    if (mat_inv_4x4(&calib_obj->s22, &s22_inv)) {
        // S22^-1*S21
        for (uint8_t r = 0; r < 6; r++) {
            for (uint8_t c = 0; c < 4; c++) {
                s22_inv_s21.m[r][c] = 0.0f;
                for (uint8_t k = 0; k < 4; k++) {
                    s22_inv_s21.m[r][c] += calib_obj->s12.m[r][k] * s22_inv.m[k][c];
                }
            }
        }
        // S11 - S12*(S22^-1*S21))
        for (uint8_t r = 0; r < 6; r++) {
            for (uint8_t c = 0; c < 6; c++) {
                s11_s12_s22_inv_s21.m[r][c] = calib_obj->s11.m[r][c];
                for (uint8_t k = 0; k < 4; k++) {
                    s11_s12_s22_inv_s21.m[r][c] -= calib_obj->s12.m[r][k] * s22_inv_s21.m[k][c];
                }
            }
        }
        // M (Reduced System Matrix)
        for (uint8_t r = 0; r < 6; r++) {
            for (uint8_t c = 0; c < 6; c++) {
                m.m[r][c] = 0.0f;
                for (uint8_t k = 0; k < 6; k++) {
                    m.m[r][c] += c_inv.m[r][k] * s11_s12_s22_inv_s21.m[k][c];
                }
            }
        }
        // eigen decomposition of M
        eigen6_t eigen;
        eigen_result_t eig_res = eig_solve_6(&m, &eigen);

        // Ellipsoid coefficients
        vector6_t coeffA = eigen.eigenvectors[eigen.largest_eigenvalue_index];
        vector6_t coeffB = eigen.eigenvectors[eigen.largest_eigenvalue_index];

        res = CALIBRATION_SUCCESS;
    } else {
        res = CALIBRATION_ERROR_S22_NOT_INVERTIBLE;
    }

    return res;
}

calibration_result_t calib_calibrate(calibration_result_data_t *calib_res, calibration_object_t *calib_obj) {
    assert(calib_res != NULL);
    assert(calib_obj != NULL);

    calibration_result_t res = calib_fit_ellipsoid(calib_res, calib_obj);

    return res;
}

void calib_capture_vec(calibration_object_t *calib_obj, vector3_t *vec) {
    assert(calib_obj != NULL);
    assert(vec != NULL);

    const mag_real_t x = vec->x;
    const mag_real_t y = vec->y;
    const mag_real_t z = vec->z;

    const mag_real_t a[6] = {
      x * x,
      y * y,
      z * z,
      2.0f * y * z,
      2.0f * x * z,
      2.0f * x * y,
    };

    const mag_real_t b[4] = {
      2.0f * x,
      2.0f * y,
      2.0f * z,
      1.0f,
    };

    for (uint8_t r = 0; r < 6; r++) {
        for (uint8_t c = 0; c < 6; c++) {
            calib_obj->s11.m[r][c] += a[r] * a[c];
        }
    }

    for (uint8_t r = 0; r < 6; r++) {
        for (uint8_t c = 0; c < 4; c++) {
            calib_obj->s12.m[r][c] += a[r] * b[c];
        }
    }

    for (uint8_t r = 0; r < 4; r++) {
        for (uint8_t c = 0; c < 4; c++) {
            calib_obj->s22.m[r][c] += b[r] * b[c];
        }
    }
}