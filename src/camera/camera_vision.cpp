
#include <yaml-cpp/yaml.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "camera_vision.hpp"

using namespace std;

void get_camera_matrix(std::string yamlFile, cv::Mat cameraMatrix) {
    std::ifstream fin(yamlFile);

    YAML::Node data = YAML::Load(fin);
    for (int i = 0; i < 3; ++i) 
        for (int j = 0; j < 3; ++j) 
            cameraMatrix.at<double>(i, j) = data["camera_matrix"][i][j].as<double>();
    fin.close();
}

void get_dist_coeff(std::string yamlFile, cv::Mat distortionCoefficients) {
    std::ifstream fin(yamlFile);
    YAML::Node data = YAML::Load(fin);

    for (int i = 0; i < 5; ++i) {
        distortionCoefficients.at<double>(i) = data["dist_coeff"][0][i].as<double>();
    }
    fin.close();
}


Camera::Camera(const char *yamlFile, unsigned int nMarkers):
    dictionary(new cv::aruco::Dictionary()),
    objPoints(cv::Mat(4, 1, CV_32FC3)), 
    cameraMatrix(cv::Mat(3, 3, CV_64F)), 
    distortionCoefficients(cv::Mat(5, 1, CV_64F)),
    isVideoStarted(0),
    nbMarkersLastDetected(0),
    rvecs(std::vector<cv::Vec3d>(nMarkers)),
    tvecs(std::vector<cv::Vec3d>(nMarkers)),
    imageCapture(cv::Mat()),
    imageCpy(cv::Mat()),
    nMarkers(nMarkers)        
    {

    float markerLength = 0.16;

    get_camera_matrix(yamlFile, this->cameraMatrix);
    get_dist_coeff(yamlFile, this->distortionCoefficients);

    this->objPoints.ptr<cv::Vec3f>(0)[0] = cv::Vec3f(-markerLength/2.f, markerLength/2.f, 0);
    this->objPoints.ptr<cv::Vec3f>(0)[1] = cv::Vec3f(markerLength/2.f, markerLength/2.f, 0);
    this->objPoints.ptr<cv::Vec3f>(0)[2] = cv::Vec3f(markerLength/2.f, -markerLength/2.f, 0);
    this->objPoints.ptr<cv::Vec3f>(0)[3] = cv::Vec3f(-markerLength/2.f, -markerLength/2.f, 0);
    //std::vector<cv::Vec3d> rvecs(nMarkers), tvecs(nMarkers);

    // cv::aruco::DetectorParameters detectorParams = cv::aruco::DetectorParameters();
    *(this->dictionary) = cv::aruco::getPredefinedDictionary(cv::aruco::DICT_6X6_50);
    // cv::Ptr<cv::aruco::Dictionary> dictionary= &_dictionary;
}


int Camera::start_video(const int cameraId) {
    std::cout << "Camera autofocus enable succeded ? -> " << this->inputVideo.set(cv::CAP_PROP_AUTOFOCUS, 0)  << std::endl;
    if (this->inputVideo.open(cameraId)) {
        this->isVideoStarted = 1;
        return 1;
    }
    return 0;
}

void Camera::close_video() {
    if (this->isVideoStarted) {
        this->inputVideo.release();
        this->isVideoStarted = 0;
    } 
}

int Camera::capture_image() {
    if (this->isVideoStarted && this->inputVideo.grab()) {
        this->inputVideo.retrieve(this->imageCapture);
        return 1;
    }
    return 0;
}


void Camera::compute_image(bool doDisplay) {

    unsigned int i = 0;
    //cv::aruco::ArucoDetector detector(dictionary, detectorParams);

    this->imageCapture.copyTo(this->imageCpy);
    std::vector<int> ids;
    std::vector<std::vector<cv::Point2f>> corners, rejected;

    // cv::Ptr<cv::aruco::Dictionary> _dict = &(this->dictionary);
    cv::aruco::detectMarkers(this->imageCapture, this->dictionary, corners, ids);
    if (ids.size() > 0) {
        solvePnP(this->objPoints, corners.at(i), this->cameraMatrix, this->distortionCoefficients, this->rvecs.at(i), this->tvecs.at(i));
        cv::aruco::drawDetectedMarkers(this->imageCpy, corners, ids);
        cv::drawFrameAxes(this->imageCpy, this->cameraMatrix, this->distortionCoefficients, this->rvecs.at(i), this->tvecs.at(i), 0.1);
        i++;
        this->nbMarkersLastDetected++;
    } else {
        this->nbMarkersLastDetected = 0;
    }
    if (doDisplay) {
        cv::imshow("out", this->imageCpy);
        int key = cv::waitKey(1);
        if (key == 27) {
            printf("Exiting trough openCV");
            close_video();
            exit(0);
        }
        // cv::destroyAllWindows();
    }
}


// void gen() {
//     cv::Mat markerImage;
//     cv::aruco::Dictionary dictionary = cv::aruco::getPredefinedDictionary(cv::aruco::DICT_6X6_250);
//     cv::aruco::drawMarker(&dictionary, 23, 200, markerImage, 1);
//     cv::imwrite("marker23.png", markerImage);
// }
