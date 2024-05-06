#ifndef __FILTER_H__
#define __FILTER_H__

#ifdef __cplusplus
extern "C" {
#endif
#include "filter_in_use.h"

int filter_create(unsigned int order, FLOATING_TYPE numerator[], FLOATING_TYPE denominator[]);
void filter_append(int filter, FLOATING_TYPE value);
FLOATING_TYPE filter_get_value(int filter);
FLOATING_TYPE filter_get_original_value(int filter);
void filter_destroy(int filter);


#ifdef __cplusplus
}
#endif

#endif