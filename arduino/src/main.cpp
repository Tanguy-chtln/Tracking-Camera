// pio run -t upload
// cu -s 115200 -l /dev/ttyACM0

#include "filter.hpp"
#include <Arduino.h>
#include <FspTimer.h>
#include "filter_in_use.h"
#define FLOATING_TYPE double

FspTimer timer_for_filter;

const int ena_pin[2] = {7, 13}, dir_pin[2] = {4, 12}, pul_pin[2] = {2,8}; int pul_value[2] = {0,0};
// const int SAMPLE_PIN = 7;
// const int PRINT_PIN = 8;
void motor_on(bool on, int motorId){ digitalWrite(ena_pin[motorId], on ? LOW: HIGH); }
void free_wheeling(int motorId){ motor_on(false, motorId); }
void move_clockwise(bool cw, int motorId){ digitalWrite(dir_pin[motorId], cw ? LOW: HIGH); }
bool is_end_of_line(char c) { return c == '\r' or c == '\n'; }
const int MAX_LENGTH = 20;
char incoming_Byte = '\r';
String text("");
const FLOATING_TYPE STEP = (1./(2*NB_MOTOR_PERIODS_PER_REVOLUTION));
const unsigned int sample_frequence = SAMPLE_FREQUENCY; // In Hz

FLOATING_TYPE command_angle[2] = {0};
FLOATING_TYPE command_norm_angle[2] = {0.};
FLOATING_TYPE motor_pwm[2] = {0.};
FLOATING_TYPE last_motor_pwm[2] = {motor_pwm[0], motor_pwm[1]};
long int nbSteps[2] = {0}; 
long int nbStepsCopy[2] = {nbSteps[0], nbSteps[1]}; 


const int order = 4;
RII_filter<FLOATING_TYPE, FILTER_ORDER> filter1{
    FILTER_NUMERATOR, // b : Numerator
    FILTER_DENOMINATOR // a : Denominator
};

const int order2 = 4;
RII_filter<FLOATING_TYPE, FILTER_ORDER> filter2{
    FILTER_NUMERATOR, // b : Numerator
    FILTER_DENOMINATOR // a : Denominator
};



bool read_line()
{
  if (is_end_of_line(incoming_Byte) or text.length() == MAX_LENGTH)
  {
    text = "";
    incoming_Byte = 'a';
  }
  while (Serial.available() and text.length() < MAX_LENGTH)
  {
    incoming_Byte = Serial.read();
    if (is_end_of_line(incoming_Byte))
    {
      break;
    }
    text += incoming_Byte;
  }
  return is_end_of_line(incoming_Byte) or text.length() == MAX_LENGTH;
}



void motor_move(int motorId) {
  if (motorId == 0) {
    filter1.append(command_norm_angle[0]);
    motor_pwm[0] = filter1.get_value();
  } else if (motorId == 1) {
    filter2.append(command_norm_angle[1]);
    motor_pwm[1] = filter2.get_value();
  }


  if ( motor_pwm[motorId] - last_motor_pwm[motorId] > STEP ) {
    nbSteps[motorId] += 1;
    last_motor_pwm[motorId] = nbSteps[motorId] * STEP;
    move_clockwise(false, motorId);
    pul_value[motorId] = 1-pul_value[motorId];
    digitalWrite(pul_pin[motorId], pul_value[motorId]);
  } else if (motor_pwm[motorId] - last_motor_pwm[motorId] < -STEP) {
    nbSteps[motorId] -= 1;
    last_motor_pwm[motorId] = nbSteps[motorId] * STEP;
    move_clockwise(true, motorId);
    pul_value[motorId] = 1-pul_value[motorId];
    digitalWrite(pul_pin[motorId], pul_value[motorId]);
  }
}

void make_a_sample(timer_callback_args_t __attribute((unused)) *p_args) {
    //digitalWrite(SAMPLE_PIN, HIGH);
    motor_move(0);
    motor_move(1);
    //digitalWrite(SAMPLE_PIN, LOW);
}


bool beginTimer(float rate) {
  uint8_t timer_type = GPT_TIMER;
  int8_t tindex = FspTimer::get_available_timer(timer_type);
  if (tindex < 0){
    tindex = FspTimer::get_available_timer(timer_type, true);
  }
  if (tindex < 0) return false;
  FspTimer::force_use_of_pwm_reserved_timer();
  if( !timer_for_filter.begin(
      TIMER_MODE_PERIODIC, timer_type, tindex, rate, 0.0f, make_a_sample
  ) ) return false;
  if (!timer_for_filter.setup_overflow_irq()) return false;
  if (!timer_for_filter.open()) return false;
  if (!timer_for_filter.start()) return false;
  return true;
}



void setup()
{

  Serial.begin(19200);
  if( !beginTimer(sample_frequence) ){
    Serial.println("Failed to init timer and interruption for filter.");
  }

  filter1.reset();
  filter2.reset();
  pinMode(ena_pin[0], OUTPUT); pinMode(pul_pin[0], OUTPUT); pinMode(dir_pin[0], OUTPUT);
  pinMode(ena_pin[1], OUTPUT); pinMode(pul_pin[1], OUTPUT); pinMode(dir_pin[1], OUTPUT);
  // pinMode(SAMPLE_PIN, OUTPUT);
  // pinMode(PRINT_PIN, OUTPUT);
  free_wheeling(0); move_clockwise(false, 0); motor_on(true, 0); pul_value[0]=0;
  free_wheeling(1); move_clockwise(false, 1); motor_on(true, 1); pul_value[1]=0;
}

FLOATING_TYPE command_norm_angle_tmp[2] = {0};
// uint32_t STEP_NB_BY_REV = 200, MSTEP_NB_BY_STEP = 2; // mode "2/B"
// uint32_t TIME_BY_REV = 1000000; // in micro seconds, tour d’une seconde
// uint32_t micro_step_delay = TIME_BY_REV/(MSTEP_NB_BY_STEP*STEP_NB_BY_REV);

void loop()
{
  while (read_line())
  {
    // digitalWrite(PRINT_PIN, HIGH);
    if (text[0] == '0') {
      command_angle[0] = text.substring(1).toFloat();
      command_norm_angle_tmp[0] = (command_angle[0]) / 360.;
      noInterrupts();
        command_norm_angle[0] = command_norm_angle_tmp[0];
        nbStepsCopy[0] = nbSteps[0];
        interrupts();
        Serial.println(nbStepsCopy[0]);
    } else if (text[0] == '1') {
      command_angle[1] = text.substring(1).toFloat();
      command_norm_angle_tmp[1] = (command_angle[1]) / 360.;
      noInterrupts();
        command_norm_angle[1] = command_norm_angle_tmp[1];
        nbStepsCopy[1] = nbSteps[1];
        interrupts();
        Serial.println(nbStepsCopy[1]);
    }
    
    // digitalWrite(PRINT_PIN, LOW);
  }  
}

