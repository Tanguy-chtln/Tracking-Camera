#include <stdio.h>
#include <stdlib.h>
#include <termios.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>

#include "SerialQueue.h"
#include "config.h"
#include "timer.h"

char *get_string_cut(char *buff, unsigned int size) {
    const char pattern_str[2] = {'\n', '\n'};
    const short pattern = *(short *) pattern_str;
    for (char *ptr = buff; ptr<buff+size-1; ptr++) {
        if (*(short *) (ptr) == pattern)
            return ptr;
    }
    return NULL;
}

void *socket_manager(void *arg) {
    struct env_t *env = (struct env_t *) arg;
    ssize_t msgReadLength;
    ssize_t msgWriteLength;
    struct string_t printString = STRING_INITIALIZER;
    char readBuffer[READ_BUFFER_SIZE+1];
    unsigned int bufferFillUp = 0;
    char *subReadBuffEndPtr = NULL;
    struct string_t writeString = STRING_INITIALIZER;

    int fd = open(SERIAL_PORT, O_RDWR | O_NOCTTY | O_NDELAY);
    if (fd == -1) {
        perror("Can't open Serial Port");
        env->isRunning = false;
        return NULL;
    }

    struct termios options;
    tcgetattr(fd, &options);
    cfsetispeed(&options, SERIAL_FREQ); // Vitesse de transmission entrante
    cfsetospeed(&options, SERIAL_FREQ); // Vitesse de transmission sortante
    options.c_cflag |= (CLOCAL | CREAD); // Enable receiver, ignore modem control lines
    options.c_cflag &= ~CSIZE;
    options.c_cflag |= CS8; // 8-bit data
    options.c_cflag &= ~PARENB; // No parity bit
    options.c_cflag &= ~CSTOPB; // 1 stop bit
    options.c_cflag &= ~CRTSCTS; // No hardware flow control

    // Set input mode (non-canonical, no echo)
    options.c_lflag &= ~(ICANON | ECHO | ECHOE | ISIG);

    // Set output mode (raw output)
    options.c_oflag &= ~OPOST;

    // Set read timeouts
    options.c_cc[VMIN] = 0;
    options.c_cc[VTIME] = 10; // 1 second timeout
    if (tcsetattr(fd, TCSANOW, &options) != 0) {
        perror("Error from tcsetattr");
        close(fd);
        env->isRunning = false;
        return NULL;
    }


    while (env->isRunning) {
        if (SOCKET_READ_ENABLE) {
            msgReadLength = read(fd, readBuffer + bufferFillUp, READ_BUFFER_SIZE);
            if (msgReadLength == -1) {
                
                //perror("Error while reading in socket");
                //exit(EXIT_FAILURE);
            } else if (msgReadLength > 0) {
                bufferFillUp += msgReadLength;
                while ((subReadBuffEndPtr = get_string_cut(readBuffer, bufferFillUp)) != NULL) {
                    string_write(&printString, readBuffer, (unsigned int) (subReadBuffEndPtr - readBuffer));
                    printf("Message retrieved : '%f'\n", atoi(printString.text)*MOTOR_STEP_ANGLE);
                    printf("Comparaison : %f %f\n", env->motorAngleSide, env->motorAngleUp);
                    env->motorAngleSide = atoi(printString.text) * MOTOR_STEP_ANGLE;
                    memcpy(readBuffer, subReadBuffEndPtr+2, bufferFillUp - (unsigned int) (subReadBuffEndPtr - readBuffer) - 2);
                    bufferFillUp -= (unsigned int) (subReadBuffEndPtr - readBuffer) + 2;
                }
            }
        }
        
        if (SOCKET_WRITE_ENABLE && buffer_queue_read(env->msgToSendQueue, &writeString)) {
            msgWriteLength = write(fd, writeString.text, writeString.textLength);
            if (msgWriteLength == -1) {
                perror("Error while writing in socket");
                exit(EXIT_FAILURE);
            } else if (msgWriteLength != writeString.textLength) {
                printf("Message was not written entierly in socket");
            } else if (writeString.text[0] == '0') {
                update_command_side(atof(writeString.text+1));
                printf("Message sent : %s\n", writeString.text);
            } else if (writeString.text[0] == '1') {
                update_command_up(atof(writeString.text+1));
                printf("Message sent : %s\n", writeString.text);
            }
        }
    }
    close(fd);
    string_destroy(&printString);
    string_destroy(&writeString);
    return NULL;
}
