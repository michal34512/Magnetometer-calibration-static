#include "matrix.h"
#include "vector.h"

typedef enum {
    CALIBRATION_SUCCESS = 0,
    CALIBRATION_ERROR_INVALID_ARGUMENTS,
    CALIBRATION_ERROR_S22_NOT_INVERTIBLE
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