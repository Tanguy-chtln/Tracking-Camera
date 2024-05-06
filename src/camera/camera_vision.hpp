#ifndef __CAMERA_VISION_H__
#define __CAMERA_VISION_H__

#include <opencv2/aruco.hpp>
#include <opencv2/opencv.hpp>

class Camera {
    private:
        cv::Ptr<cv::aruco::Dictionary> dictionary;
        cv::Mat objPoints;
        cv::Mat cameraMatrix;
        cv::Mat distortionCoefficients;
        cv::VideoCapture inputVideo;
        int isVideoStarted;

    public: 
        int nbMarkersLastDetected;
        std::vector<cv::Vec3d> rvecs;
        std::vector<cv::Vec3d> tvecs;
        cv::Mat imageCapture;
        cv::Mat imageCpy;
        unsigned int nMarkers;

        Camera(const char *yamlFile, unsigned int nMarkers);
        int start_video(const int cameraId);
        void close_video();
        void compute_image(bool doDisplay);
        int capture_image();

};

#endif