
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "config.h"
#include "SerialQueue.h"
#include "camera.h"
#include "timer.h"

void *aruko_finder(void *arg) {
    struct env_t *env = (struct env_t *) arg;
    int sprintfOutput = 3;
    char buf[255] = {"get\0"};
    int * camera;
    double translationArray[3];
    double rotationArray[3];
    unsigned int dataReceived=0;

    double motorAngleCpySide = 0.;
    double relAngleSide;
    double absAngleSide;

    double relAngleUp;
    double absAngleUp;
    double motorAngleCpyUp = 0.;

    printf("Starting camera...\n");
    camera = initialize_camera(CAMERA_CONFIG_FILE, 1);
    start_video(camera, CAMERA_ID);
    printf("Camera started\n");
    if ((capture_image(camera) && env->isRunning)) 
        env->isCameraRunning = true;
    
    while (capture_image(camera) && env->isRunning && env->isCameraRunning) {
        compute_image(camera, true);
        dataReceived = get_values(camera, translationArray, rotationArray);
        
        relAngleSide =-45 * (atan(translationArray[0]/translationArray[2]) - (.45/2)) / (.45/2);
        absAngleSide = relAngleSide + motorAngleCpySide;
        relAngleUp =-20 * (atan(translationArray[1]/translationArray[2]) - (.20/2)) / (.20/2);
        absAngleUp = relAngleUp + motorAngleCpyUp;
        // printf("%f %f %f %f\n", translationArray[0], translationArray[1], translationArray[2], absAngleUp);
        if (dataReceived > 0) {
            sprintfOutput = sprintf(buf, "0%lf\n", absAngleSide);
            buffer_queue_write(env->msgToSendQueue, buf, sprintfOutput);
            sprintfOutput = sprintf(buf, "1%lf\n", absAngleUp);
            buffer_queue_write(env->msgToSendQueue, buf, sprintfOutput);
        } 
        motorAngleCpySide = env->motorAngleSide;
        motorAngleCpyUp = env->motorAngleUp;
    }
    env->isRunning = false;
    printf("Closing camera\n");
    delete_camera(camera);
    printf("Camera closed\n");
    return NULL;
}