#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "magcalib.h"

#define MIN_GAIN   0.5f
#define MAX_GAIN   10
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
        generated[i].x = ((float)rand() / (float)RAND_MAX) * 2 - 1;
        generated[i].y = ((float)rand() / (float)RAND_MAX) * 2 - 1;
        generated[i].z = ((float)rand() / (float)RAND_MAX) * 2 - 1;
        vec_normalize(&generated[i]);

        // Add noise
        float scalar = (((float)rand() / (float)RAND_MAX) * 2 - 1) * MAX_NOISE + 1;
        vec_mul_scalar(&generated[i], scalar);

        // Apply ellipsoid transformation (rotation, scaling, offset)
        generated[i].x *= gainX;
        generated[i].y *= gainY;
        generated[i].z *= gainZ;

        vec_rotate_xyz(&generated[i], angleX, angleY, angleZ);

        generated[i].x += offsetX;
        generated[i].y += offsetY;
        generated[i].z += offsetZ;
    }
}

int main() {
    srand(42);

    // Simulate random magnetometer data
    generate_mag_data(DATA_LEN);
    FILE *out = fopen("scripts/samples.csv", "w");
    if (out == NULL) {
        perror("Failed to open scripts/samples.csv");
        return 1;
    }
    fprintf(out, "i,x,y,z\n");
    for (uint16_t i = 0; i < DATA_LEN; i++) {
        fprintf(out, "%u,%.9f,%.9f,%.9f\n", (unsigned)i, generated[i].x, generated[i].y, generated[i].z);
    }

    // Callibration
    calibration_object_t calib_obj = {0};
    calibration_result_data_t calib_res;
    for (uint16_t i = 0; i < DATA_LEN; i++) {
        calib_capture_vec(&calib_obj, &generated[i]);
    }
    calib_calibrate(&calib_res, &calib_obj);

    fclose(out);
    return 0;
}
