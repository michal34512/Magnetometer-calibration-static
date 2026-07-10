#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "calibration.h"
#include "eig.h"
#include "magcalib.h"

#define ALPHA 16

typedef struct {
    matrix3_t m;
    vector3_t pqr;
    mag_real_t d;
} calibration_ctx_t;

const matrix6_t c_inv = {
    .m = {
        {(mag_real_t)(ALPHA - 4) / (-ALPHA * ALPHA + 3 * ALPHA), -(mag_real_t)(ALPHA - 2) / (-ALPHA * ALPHA + 3 * ALPHA), -(mag_real_t)(ALPHA - 2) / (-ALPHA * ALPHA + 3 * ALPHA), 0.f, 0.f, 0.f},
        {-(mag_real_t)(ALPHA - 2) / (-ALPHA * ALPHA + 3 * ALPHA), (mag_real_t)(ALPHA - 4) / (-ALPHA * ALPHA + 3 * ALPHA), -(mag_real_t)(ALPHA - 2) / (-ALPHA * ALPHA + 3 * ALPHA), 0.f, 0.f, 0.f},
        {-(mag_real_t)(ALPHA - 2) / (-ALPHA * ALPHA + 3 * ALPHA), -(mag_real_t)(ALPHA - 2) / (-ALPHA * ALPHA + 3 * ALPHA), (mag_real_t)(ALPHA - 4) / (-ALPHA * ALPHA + 3 * ALPHA)   , 0.f, 0.f, 0.f},
        {0.f, 0.f, 0.f, -1.f / ALPHA, 0.f, 0.f},
        {0.f, 0.f, 0.f, 0.f, -1.f / ALPHA, 0.f},
        {0.f, 0.f, 0.f, 0.f, 0.f, -1.f / ALPHA},
    },
};

static inline void
calib_fit_ellipsoid_decompose_coeff(eigen6_t *eigen, matrix6x4_t *s22_inv_s21, calibration_ctx_t *ctx) {
    const vector6_t coeffA = eigen->eigenvectors[eigen->largest_eigenvalue_index];

    // coefB = - s22_inv_s21 * coeffA
    vector6_t coeffB;
    for (uint8_t r = 0; r < 4; r++) {
        coeffB.v[r] = 0.0f;
        for (uint8_t k = 0; k < 6; k++) {
            coeffB.v[r] -= s22_inv_s21->m[k][r] * coeffA.v[k];
        }
    }

    // coefB, coefA -> ctx->m, ctx->pqr, ctx->d
    ctx->m.m[0][0] = coeffA.v[0];
    ctx->m.m[0][1] = coeffA.v[5];
    ctx->m.m[0][2] = coeffA.v[4];
    ctx->m.m[1][0] = coeffA.v[5];
    ctx->m.m[1][1] = coeffA.v[1];
    ctx->m.m[1][2] = coeffA.v[3];
    ctx->m.m[2][0] = coeffA.v[4];
    ctx->m.m[2][1] = coeffA.v[3];
    ctx->m.m[2][2] = coeffA.v[2];

    ctx->pqr.v[0] = coeffB.v[0];
    ctx->pqr.v[1] = coeffB.v[1];
    ctx->pqr.v[2] = coeffB.v[2];
    ctx->d        = coeffB.v[3];
}

calibration_result_t calib_fit_ellipsoid(calibration_object_t *calib_obj, calibration_ctx_t *ctx) {
    assert(calib_obj != NULL);
    assert(ctx != NULL);

    // M = C1^-1*(S11-S12*(S22^-1*S21));
    calibration_result_t res;
    matrix4_t inv_s22;
    matrix6x4_t inv_s22_s21;
    matrix6_t s11_s12_inv_s22_s21;
    matrix6_t m;

    if (mat_inv_4x4(&calib_obj->s22, &inv_s22)) {
        // S22^-1*S21
        for (uint8_t r = 0; r < 6; r++) {
            for (uint8_t c = 0; c < 4; c++) {
                inv_s22_s21.m[r][c] = 0.0f;
                for (uint8_t k = 0; k < 4; k++) {
                    inv_s22_s21.m[r][c] += calib_obj->s12.m[r][k] * inv_s22.m[k][c];
                }
            }
        }
        // S11 - S12*(S22^-1*S21))
        for (uint8_t r = 0; r < 6; r++) {
            for (uint8_t c = 0; c < 6; c++) {
                s11_s12_inv_s22_s21.m[r][c] = calib_obj->s11.m[r][c];
                for (uint8_t k = 0; k < 4; k++) {
                    s11_s12_inv_s22_s21.m[r][c] -= calib_obj->s12.m[r][k] * inv_s22_s21.m[c][k];
                }
            }
        }
        // M (Reduced System Matrix)
        for (uint8_t r = 0; r < 6; r++) {
            for (uint8_t c = 0; c < 6; c++) {
                m.m[r][c] = 0.0f;
                for (uint8_t k = 0; k < 6; k++) {
                    m.m[r][c] += c_inv.m[r][k] * s11_s12_inv_s22_s21.m[k][c];
                }
            }
        }

        eigen6_t eigen;
        eigen_result_t eig_res = eig_solve_6(&m, &eigen);

        // Ellipsoid coefficients
        switch (eig_res) {
        case EIGEN_MAX_ITERATION:
            res = CALIBRATION_WARNING_EIGEN_MAX_ITERATION;  // decompose the coefficients anyway
            /* fall through */
        case EIGEN_SUCCESS:
            calib_fit_ellipsoid_decompose_coeff(&eigen, &inv_s22_s21, ctx);
            break;
        case EIGEN_SINGULARITY:
            res = CALIBRATION_ERROR_EIGEN_SINGULARITY;
            break;
        default:
            res = CALIBRATION_ERROR_UNKNOWN;
            break;
        }
        res = CALIBRATION_SUCCESS;
    } else {
        res = CALIBRATION_ERROR_S22_NOT_INVERTIBLE;
    }

    return res;
}

static calibration_result_t calib_fit_get_transform(calibration_ctx_t *ctx, calibration_result_data_t *calib_res) {
    assert(calib_res != NULL);
    assert(ctx != NULL);

    // Rotation
    eigen3_t eig;
    eig_solve_3(&ctx->m, &eig);

    // Orthogonalize the eigenvectors to ensure they are orthonormal
    eig_orthogonalize_3(&eig);

    // Calculate semi-axis lengths (gain matrix G)
    vector3_t pqr_r;
    for (uint8_t c = 0; c < 3; c++) {
        pqr_r.v[c] = 0.0f;
        for (uint8_t k = 0; k < 3; k++) {
            pqr_r.v[c] += ctx->pqr.v[k] * eig.eigenvectors[c].v[k];
        }
    }
    const mag_real_t a_ = eig.eigenvalues[0];
    const mag_real_t b_ = eig.eigenvalues[1];
    const mag_real_t c_ = eig.eigenvalues[2];
    const mag_real_t p_ = pqr_r.v[0];
    const mag_real_t q_ = pqr_r.v[1];
    const mag_real_t r_ = pqr_r.v[2];
    const mag_real_t d_ = ctx->d;

    matrix3_t gain = {0};
    // we do not accctually calculate diagonal G to avoid declaring additional matrix
    gain.m[0][0] = 1 / sqrt((p_ * p_) / (a_ * a_) + (q_ * q_) / (a_ * b_) + (r_ * r_) / (a_ * c_) - d_ / a_);
    gain.m[0][1] = 1 / sqrt((p_ * p_) / (a_ * b_) + (q_ * q_) / (b_ * b_) + (r_ * r_) / (b_ * c_) - d_ / b_);
    gain.m[0][2] = 1 / sqrt((p_ * p_) / (a_ * c_) + (q_ * q_) / (b_ * c_) + (r_ * r_) / (c_ * c_) - d_ / c_);

    // Rotate gain matrix to align with the ellipsoid axes (Q^T * G * Q)
    // Q^T * G
    for (uint8_t r = 1; r < 3; r++) {
        for (uint8_t c = 0; c < 3; c++) {
            gain.m[r][c] = gain.m[0][c] * eig.eigenvectors[c].v[r];
        }
    }
    for (uint8_t c = 0; c < 3; c++) {
        gain.m[0][c] *= eig.eigenvectors[c].v[0];
    }
    for (uint8_t r = 0; r < 3; r++) {
        for (uint8_t c = 0; c < 3; c++) {
            calib_res->transform.m[r][c] = 0;
            for (uint8_t k = 0; k < 3; k++) {
                calib_res->transform.m[r][c] += gain.m[r][k] * eig.eigenvectors[k].v[c];
            }
        }
    }
    return CALIBRATION_SUCCESS;
}

static calibration_result_t calib_fit_get_offset(calibration_ctx_t *ctx, calibration_result_data_t *calib_res) {
    if (mat_inv_3x3(&ctx->m)) {
        for (uint8_t r = 0; r < 3; r++) {
            calib_res->offset.v[r] = 0.0f;
            for (uint8_t c = 0; c < 3; c++) {
                calib_res->offset.v[r] -= ctx->m.m[r][c] * ctx->pqr.v[c];
            }
        }
    } else {
        return CALIBRATION_ERROR_ELLIP_MAT_NOT_INVERTIBLE;
    }
    return CALIBRATION_SUCCESS;
}

calibration_result_t calib_calibrate(calibration_result_data_t *calib_res, calibration_object_t *calib_obj) {
    assert(calib_res != NULL);
    assert(calib_obj != NULL);

    calibration_ctx_t ctx;

    calibration_result_t res = calib_fit_ellipsoid(calib_obj, &ctx);
    res                      = calib_fit_get_transform(&ctx, calib_res);
    if (res != CALIBRATION_SUCCESS) {
        return res;
    }
    res = calib_fit_get_offset(&ctx, calib_res);

    return res;
}

void calib_capture_vec(calibration_object_t *calib_obj, vector3_t *vec) {
    assert(calib_obj != NULL);
    assert(vec != NULL);

    const mag_real_t x = vec->v[0];
    const mag_real_t y = vec->v[1];
    const mag_real_t z = vec->v[2];

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

void calib_calibrate_vec(calibration_result_data_t *calib_res, vector3_t *vec) {
    for (uint8_t r = 0; r < 3; r++) {
        vec->v[r] -= calib_res->offset.v[r];
    }

    mag_real_t v0 = vec->v[0];
    mag_real_t v1 = vec->v[1];
    vec->v[0]     = calib_res->transform.m[0][0] * v0 + calib_res->transform.m[0][1] * v1
                + calib_res->transform.m[0][2] * vec->v[2];
    vec->v[1] = calib_res->transform.m[1][0] * v0 + calib_res->transform.m[1][1] * v1
                + calib_res->transform.m[1][2] * vec->v[2];
    vec->v[2] = calib_res->transform.m[2][0] * v0 + calib_res->transform.m[2][1] * v1
                + calib_res->transform.m[2][2] * vec->v[2];
}