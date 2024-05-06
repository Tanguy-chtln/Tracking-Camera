#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <asm/unistd.h>
#include <bits/types/timer_t.h>

#include "filter_in_use.h"
#include "config.h"
#include "filter.h"


int filterSide;
int filterUp;

double lastCommandSide = 0;
double lastCommandUp = 0;

struct env_t *glbEnv = NULL;

void update_command_side(double command) {
    lastCommandSide = command;
}

void update_command_up(double command) {
    lastCommandUp = command;
}

double convert_command(int command) {
    return command / 360.;
}

void filter_update(int) {
    filter_append(filterSide, convert_command(lastCommandSide));
    glbEnv->motorAngleSide = filter_get_value(filterSide) * 360.;
    filter_append(filterUp, convert_command(lastCommandUp));
    glbEnv->motorAngleUp = filter_get_value(filterUp) * 360.;
}

timer_t start_timer() {
    struct sigevent sev;
    timer_t timerid;
    struct itimerspec its;

    // Créer une minuterie POSIX
    sev.sigev_notify = SIGEV_SIGNAL;
    sev.sigev_signo = SIGALRM;
    sev.sigev_value.sival_ptr = &timerid;
    if (timer_create(CLOCK_REALTIME, &sev, &timerid) == -1) {
        perror("timer_create");
        exit(EXIT_FAILURE);
    }

    // Configurer la période de la minuterie
    its.it_value.tv_sec = 0;
    its.it_value.tv_nsec = 1000000000/SAMPLE_FREQUENCY;
    its.it_interval.tv_sec = 0;
    its.it_interval.tv_nsec = 1000000000/SAMPLE_FREQUENCY;

    // Démarrer la minuterie
    if (timer_settime(timerid, 0, &its, NULL) == -1) {
        perror("timer_settime");
        exit(EXIT_FAILURE);
    }

    // Enregistrer le gestionnaire de signal
    signal(SIGALRM, filter_update);
    return timerid;
}

timer_t start_filter(struct env_t *arg) {
    glbEnv = (struct env_t *) arg;
    FLOATING_TYPE numerator[FILTER_ORDER+1] = FILTER_NUMERATOR;
    FLOATING_TYPE denominator[FILTER_ORDER+1] = FILTER_DENOMINATOR;

    filterSide = filter_create(FILTER_ORDER, numerator, denominator);
    filterUp = filter_create(FILTER_ORDER, numerator, denominator);

    return start_timer();
}


void stop_filter(timer_t id) {
    if (timer_delete(id) == -1) {
        perror("timer_delete");
        exit(EXIT_FAILURE);
    }
    sleep(1);
    signal(SIGALRM, SIG_DFL);
    filter_destroy(filterUp);
    filter_destroy(filterSide);
}
