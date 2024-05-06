#ifndef __CAMERA_H__
#define CAMERA_H__

#ifdef __cplusplus
extern "C" {
#endif
#include <libheif/heif.h>
#include <stdbool.h>
int *initialize_camera(const char *yamlFile, unsigned int nMarkers);
int start_video(int *camera, const int cameraId);
int capture_image(int *camera);
void compute_image(int *camera, bool doDisplay);
unsigned int get_values(int *camera, double *transVals, double *rotVals);
void delete_camera(int *camera);
void convertToMat(int *camera, unsigned char *img, int width, int height);

#ifdef __cplusplus
}
#endif

#endif