#ifndef __TIMER_H__
#define __TIMER_H__
#include <bits/types/timer_t.h>

void update_command_side(int command);
void update_command_up(double command);
timer_t start_filter(struct env_t *arg);
void stop_filter(timer_t id);

#endif