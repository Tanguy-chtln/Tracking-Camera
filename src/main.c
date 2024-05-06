#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <stdbool.h>
#include <unistd.h>

#include "SerialQueue.h"
#include "filter_in_use.h"
#include "config.h"
#include "socket.h"
#include "aruko_detect.h"
#include "timer.h"



int main() {
    struct env_t env = {false, NULL, false, 0., 0.};
    char inputStream[CMD_LINE_MAX_LENGTH+1];
    int inputAngle;
    int orderMsgSize;
    pthread_t socketThread;
    pthread_t cameraThread;
    timer_t filterTimer;
    env.msgToSendQueue = buffer_queue_create(QUEUE);

    env.isRunning =true;

    pthread_create(&socketThread, NULL, socket_manager, (void *) &env);
    printf("Beginning start routine\n");
    buffer_queue_write(env.msgToSendQueue, "045\n", 4);
    buffer_queue_write(env.msgToSendQueue, "130\n", 4);
    sleep(1);
    if (env.isRunning) {
        sleep(3);
        buffer_queue_write(env.msgToSendQueue, "00\n", 3);
        buffer_queue_write(env.msgToSendQueue, "10\n", 3);
        sleep(4);
        printf("End of start routine\n");
        filterTimer = start_filter(&env);

        if (CAMERA_ENABLE) {
            pthread_create(&cameraThread, NULL, aruko_finder, (void *) &env);
            while (!env.isCameraRunning && env.isRunning);
        }

        while (env.isRunning) { // Tracker shuts down when stop is written in the terminal
            fgets(inputStream, CMD_LINE_MAX_LENGTH, stdin);
            if (!strcmp(inputStream, "stop\n")) {
                printf("Stopping...\n");
                env.isRunning = false;
            } else if (!strcmp(inputStream, "cam stop\n")) {
                printf("Stopping camera");
                env.isCameraRunning = false;
            }
             else if ( strlen(inputStream) == 1 || (strlen(inputStream) > 1 && (inputAngle = atoi(inputStream+1)) == 0 && inputStream[1] != '0')) {
                printf("Wrong input : Enter either a angle (integer) or 'stop'\n");
            } else if(inputStream[0] == 'a') {
                orderMsgSize = sprintf(inputStream, "0%d\n", inputAngle);
                buffer_queue_write(env.msgToSendQueue, inputStream, orderMsgSize);
            } else if(inputStream[0] == 'b') {
                orderMsgSize = sprintf(inputStream, "1%d\n", inputAngle);
                buffer_queue_write(env.msgToSendQueue, inputStream, orderMsgSize);
            }       
        }
        if (CAMERA_ENABLE) 
            pthread_join(cameraThread, NULL);
        
        stop_filter(filterTimer);
    }
    pthread_join(socketThread, NULL);
    
    buffer_queue_destroy(env.msgToSendQueue);
    return EXIT_SUCCESS;
}
