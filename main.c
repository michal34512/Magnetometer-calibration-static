#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "magcalib.h"

#define MIN_GAIN   40
#define MAX_GAIN   1
#define MIN_OFFSET -10
#define MAX_OFFSET 10
#define MAX_NOISE  0.1f
#define DATA_LEN   1000

vector3_t generated[DATA_LEN] = {0};

void generate_mag_data(int n) {
    // Random rot ation matrix
    float angleX = ((float)rand() / (float)RAND_MAX) * 2 * M_PI;
    float angleY = ((float)rand() / (float)RAND_MAX) * 2 * M_PI;
    float angleZ = ((float)rand() / (float)RAND_MAX) * 2 * M_PI;

    float gainX = ((float)rand() / (float)RAND_MAX) * (MAX_GAIN - MIN_GAIN) + MIN_GAIN;
    float gainY = ((float)rand() / (float)RAND_MAX) * (MAX_GAIN - MIN_GAIN) + MIN_GAIN;
    float gainZ = ((float)rand() / (float)RAND_MAX) * (MAX_GAIN - MIN_GAIN) + MIN_GAIN;

    float offsetX = ((float)rand() / (float)RAND_MAX) * (MAX_OFFSET - MIN_OFFSET) + MIN_OFFSET;
    float offsetY = ((float)rand() / (float)RAND_MAX) * (MAX_OFFSET - MIN_OFFSET) + MIN_OFFSET;
    float offsetZ = ((float)rand() / (float)RAND_MAX) * (MAX_OFFSET - MIN_OFFSET) + MIN_OFFSET;

    // Generate data
    for (uint16_t i = 0; i < n; i++) {
        generated[i].v[0] = ((float)rand() / (float)RAND_MAX) * 2 - 1;
        generated[i].v[1] = ((float)rand() / (float)RAND_MAX) * 2 - 1;
        generated[i].v[2] = ((float)rand() / (float)RAND_MAX) * 2 - 1;
        vec_normalize(&generated[i]);

        // Add noise
        float scalar = (((float)rand() / (float)RAND_MAX) * 2 - 1) * MAX_NOISE + 1;
        vec_mul_scalar(&generated[i], scalar);

        // Apply ellipsoid transformation (rotation, scaling, offset)
        generated[i].v[0] *= gainX;
        generated[i].v[1] *= gainY;
        generated[i].v[2] *= gainZ;

        vec_rotate_xyz(&generated[i], angleX, angleY, angleZ);

        generated[i].v[0] += offsetX;
        generated[i].v[1] += offsetY;
        generated[i].v[2] += offsetZ;
    }
}

void generate_mag_data_single_axis_rotation(int n) {
    // Choose random rotation axis (0=X, 1=Y, 2=Z)
    int axis = rand() % 3;

    // Generate rotation angle around the chosen axis only
    float angle = ((float)rand() / (float)RAND_MAX) * 2 * M_PI;

    float gainX = ((float)rand() / (float)RAND_MAX) * (MAX_GAIN - MIN_GAIN) + MIN_GAIN;
    float gainY = ((float)rand() / (float)RAND_MAX) * (MAX_GAIN - MIN_GAIN) + MIN_GAIN;
    float gainZ = ((float)rand() / (float)RAND_MAX) * (MAX_GAIN - MIN_GAIN) + MIN_GAIN;

    float offsetX = ((float)rand() / (float)RAND_MAX) * (MAX_OFFSET - MIN_OFFSET) + MIN_OFFSET;
    float offsetY = ((float)rand() / (float)RAND_MAX) * (MAX_OFFSET - MIN_OFFSET) + MIN_OFFSET;
    float offsetZ = ((float)rand() / (float)RAND_MAX) * (MAX_OFFSET - MIN_OFFSET) + MIN_OFFSET;

    // Generate data
    for (uint16_t i = 0; i < n; i++) {
        generated[i].v[0] = ((float)rand() / (float)RAND_MAX) * 2 - 1;
        generated[i].v[1] = ((float)rand() / (float)RAND_MAX) * 2 - 1;
        generated[i].v[2] = ((float)rand() / (float)RAND_MAX) * 2 - 1;
        vec_normalize(&generated[i]);

        // Add noise
        float scalar = (((float)rand() / (float)RAND_MAX) * 2 - 1) * MAX_NOISE + 1;
        vec_mul_scalar(&generated[i], scalar);

        // Apply scaling (gains) - creates ellipsoid shape
        generated[i].v[0] *= gainX;
        generated[i].v[1] *= gainY;
        generated[i].v[2] *= gainZ;

        // Apply rotation around single axis
        if (axis == 0) {
            // Rotation around X axis
            vec_rotate_xyz(&generated[i], angle, 0, 0);
        } else if (axis == 1) {
            // Rotation around Y axis
            vec_rotate_xyz(&generated[i], 0, angle, 0);
        } else {
            // Rotation around Z axis
            vec_rotate_xyz(&generated[i], 0, 0, angle);
        }

        // Apply offset
        generated[i].v[0] += offsetX;
        generated[i].v[1] += offsetY;
        generated[i].v[2] += offsetZ;
    }
}

static void save_generated_data_to_csv(char *filename, vector3_t *data, size_t data_len) {
    FILE *out = fopen(filename, "w");
    if (out == NULL) {
        perror("Failed to open file for writing");
        return;
    }

    fprintf(out, "i,x,y,z\n");
    for (uint16_t i = 0; i < data_len; i++) {
        fprintf(out, "%u,%.9f,%.9f,%.9f\n", (unsigned)i, data[i].v[0], data[i].v[1], data[i].v[2]);
    }
    fclose(out);
}

static void save_ellipsoid(char *filename, calibration_result_data_t *calib_res) {
    FILE *out = fopen(filename, "w");
    if (out == NULL) {
        perror("Failed to open file for writing");
        return;
    }

    fprintf(out, "# Ellipsoid parameters used in calibration\n");
    fprintf(out, "offset,%.9f,%.9f,%.9f\n", calib_res->offset.v[0], calib_res->offset.v[1], calib_res->offset.v[2]);
    fprintf(out,
            "transform_row0,%.9f,%.9f,%.9f\n",
            calib_res->transform.m[0][0],
            calib_res->transform.m[0][1],
            calib_res->transform.m[0][2]);
    fprintf(out,
            "transform_row1,%.9f,%.9f,%.9f\n",
            calib_res->transform.m[1][0],
            calib_res->transform.m[1][1],
            calib_res->transform.m[1][2]);
    fprintf(out,
            "transform_row2,%.9f,%.9f,%.9f\n",
            calib_res->transform.m[2][0],
            calib_res->transform.m[2][1],
            calib_res->transform.m[2][2]);

    fclose(out);
}

int main() {
    srand(42);

    // Simulate random magnetometer data
    generate_mag_data(DATA_LEN);
    save_generated_data_to_csv("scripts/samples.csv", generated, DATA_LEN);

    // Callibration
    calibration_object_t calib_obj = {0};
    calibration_result_data_t calib_res;
    for (uint16_t i = 0; i < DATA_LEN; i++) {
        calib_capture_vec(&calib_obj, &generated[i]);
    }
    calib_calibrate(&calib_res, &calib_obj);
    save_ellipsoid("scripts/ellipsoid.csv", &calib_res);

    // Calibrate data
    mat_print_3(&calib_res.transform);
    vec_print_3(&calib_res.offset);
    for (uint16_t i = 0; i < DATA_LEN; i++) {
        calib_calibrate_vec(&calib_res, &generated[i]);
    }
    save_generated_data_to_csv("scripts/calibrated.csv", generated, DATA_LEN);
    return 0;
}
