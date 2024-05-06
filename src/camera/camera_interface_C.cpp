#include "camera.h"
#include "camera_vision.hpp"
#include <libheif/heif.h>


    
int *initialize_camera(const char *yamlFile, unsigned int nMarkers) {
    Camera *camera = new Camera(yamlFile, nMarkers); 
    return (int *) camera;
}

int start_video(int *camera, const int cameraId) {
    return ((Camera *) camera)->start_video(cameraId);
}

int capture_image(int *camera) {
    return ((Camera *) camera)->capture_image();
}

void compute_image(int *camera, bool doDisplay) {
    return ((Camera *) camera)->compute_image(doDisplay);
}

unsigned int get_values(int *_camera, double *transVals, double *rotVals) {
    Camera *camera = (Camera *) _camera;
    for (long unsigned int i=0; i<camera->tvecs.size(); i++) {
        transVals[i*3+0] = camera->tvecs[i][0];
        transVals[i*3+1] = camera->tvecs[i][1];
        transVals[i*3+2] = camera->tvecs[i][2];

        rotVals[i*3+0] = camera->rvecs[i][0];
        rotVals[i*3+1] = camera->rvecs[i][1];
        rotVals[i*3+2] = camera->rvecs[i][2];
    }
    return camera->nbMarkersLastDetected;
}

void delete_camera(int *camera) {
    delete ((Camera *) camera);
}

void displayMat(const cv::Mat& mat) {
    // Create a window to display the image
    cv::namedWindow("Display Window", cv::WINDOW_NORMAL);

    // Show the image in the window
    cv::imshow("Display Window", mat);

    // Wait for a key press to close the window
    cv::waitKey(1);
}


void convertToMat(int *_camera, unsigned char *img, int width, int height) {
    Camera *camera = (Camera *) _camera;
    // Create a cv::Mat object of appropriate size and type
    camera->imageCapture = cv::Mat(height, width, CV_8UC3);

    // Iterate over each row and each column of the image
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // Access the pixel in the input image
            unsigned char *pixelPtr = img + (y * width + x) * 3;

            // Assign pixel values to the corresponding location in the cv::Mat
            camera->imageCapture.at<cv::Vec3b>(y, x)[0] = pixelPtr[2];  // Blue
            camera->imageCapture.at<cv::Vec3b>(y, x)[1] = pixelPtr[1];  // Green
            camera->imageCapture.at<cv::Vec3b>(y, x)[2] = pixelPtr[0];  // Red
        }
    }

}

