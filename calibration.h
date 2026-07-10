#ifndef CALIBRATION_H
#define CALIBRATION_H

#include "matrix.h"
#include "vector.h"

typedef enum {
    /* Success */
    CALIBRATION_SUCCESS = 0,

    /* Warnings*/
    CALIBRATION_WARNING_EIGEN_MAX_ITERATION = 1,

    /* Errors */
    CALIBRATION_ERROR_INVALID_ARGUMENTS        = 100,
    CALIBRATION_ERROR_S22_NOT_INVERTIBLE       = 101,
    CALIBRATION_ERROR_ELLIP_MAT_NOT_INVERTIBLE = 102,
    CALIBRATION_ERROR_EIGEN_SINGULARITY        = 103,
    CALIBRATION_ERROR_UNKNOWN                  = 104,
} calibration_result_t;

typedef struct {
    vector3_t offset;
    matrix3_t transform;
} calibration_result_data_t;

typedef struct {
    matrix6_t s11;
    matrix6x4_t s12;
    matrix4_t s22;
} calibration_object_t;

calibration_result_t calib_calibrate(calibration_result_data_t *calib_res, calibration_object_t *calib_obj);

void calib_capture_vec(calibration_object_t *calib_obj, vector3_t *vec);

void calib_calibrate_vec(calibration_result_data_t *calib_res, vector3_t *vec);

#endif  // CALIBRATION_H