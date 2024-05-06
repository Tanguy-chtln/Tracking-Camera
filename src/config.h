#ifndef __CONFIG_H__
#define __CONFIG_H__
#include <stdbool.h>
#include <stdint.h>
#include "SerialQueue.h"
#include "filter_in_use.h"

#define READ_BUFFER_SIZE 20
#define CMD_LINE_MAX_LENGTH 20
#define SERIAL_FREQ B19200
#define SERIAL_PORT "/dev/ttyACM0"
#define SOCKET_READ_ENABLE true
#define SOCKET_WRITE_ENABLE true

#define CAMERA_ENABLE true
#define CAMERA_ID 2
#define CAMERA_CONFIG_FILE "config/camera.yaml" 
#define MOTOR_STEP_ANGLE (360./(2*NB_MOTOR_PERIODS_PER_REVOLUTION))


struct env_t {
    bool isRunning;
    struct buffer_queue_t *msgToSendQueue;
    bool isCameraRunning;
    double motorAngleSide;
    double motorAngleUp;
};

#endif