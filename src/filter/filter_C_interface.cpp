#include "filter.h"
#include "filter.hpp"
#include <iostream>

int **filters = NULL;
int *filtersOrders = NULL;
int nbFilters = 0;
int nbFiltersMax = 0;

int filter_create(unsigned int order, FLOATING_TYPE numerator[], FLOATING_TYPE denominator[]) {
    if (filters == NULL) {
        nbFiltersMax = 2;
        filters = (int **) calloc(nbFiltersMax, sizeof(int *));
        filtersOrders = (int *) calloc(nbFiltersMax, sizeof(int));
    } else if (nbFiltersMax >= nbFilters) {
        nbFiltersMax *= 2;
        filters = (int **) realloc(filters, sizeof(int *) * nbFiltersMax);
        filtersOrders = (int *) realloc(filtersOrders, sizeof(int) * nbFiltersMax);
    }

    int *filterOrderSlot = filtersOrders+nbFilters;
     
    switch (order) {
        case 1:
             filters[nbFilters] = (int *) new RII_filter<FLOATING_TYPE, 1>{
                {numerator[0], numerator[1]}, // b : Numerator
                {denominator[0], denominator[1]}}; // a : Denominator
            *filterOrderSlot = 1;
            break;

        case 2:
            filters[nbFilters] = (int *) new RII_filter<FLOATING_TYPE, 2>{
                {numerator[0], numerator[1], numerator[2] }, // b : Numerator
                {denominator[0], denominator[1], denominator[2]}}; // a : Denominator
            *filterOrderSlot = 2;
            break;

        case 3:
            filters[nbFilters] = (int *) new RII_filter<FLOATING_TYPE, 3>{
                {numerator[0], numerator[1], numerator[2], numerator[3] }, // b : Numerator
                {denominator[0], denominator[1], denominator[2], denominator[3]}}; // a : Denominator
            *filterOrderSlot = 3;
            break;
        
        case 4:
            filters[nbFilters] = (int *) new RII_filter<FLOATING_TYPE, 4>{
                {numerator[0], numerator[1], numerator[2], numerator[3],numerator[4] }, // b : Numerator
                {denominator[0], denominator[1], denominator[2], denominator[3], denominator[4]}}; // a : Denominator
            *filterOrderSlot = 4;
            break;
        
        case 5:
            filters[nbFilters] = (int *) new RII_filter<FLOATING_TYPE, 5>{
                {numerator[0], numerator[1], numerator[2], numerator[3],numerator[4], numerator[5] }, // b : Numerator
                {denominator[0], denominator[1], denominator[2], denominator[3], denominator[4], denominator[5]}}; // a : Denominator
            *filterOrderSlot = 5;
            break;

        case 6:
            filters[nbFilters] = (int *) new RII_filter<FLOATING_TYPE, 6>{
                {numerator[0], numerator[1], numerator[2], numerator[3], numerator[4], numerator[5], numerator[6] }, // b : Numerator
                {denominator[0], denominator[1], denominator[2], denominator[3], denominator[4], denominator[5], denominator[6]}}; // a : Denominator
            *filterOrderSlot = 6;
            break;

        default:
            std::cout << "Wrong filter order given" << std::endl;
            exit(EXIT_FAILURE);
    }
    return nbFilters++;
}


void filter_append(int filter, FLOATING_TYPE value) {
    int filterOrder = filtersOrders[filter];
    if (filter >= nbFilters) {
        std::cout << "Wrong filter given" << std::endl;
        exit(EXIT_FAILURE);
    }

    switch(filterOrder) {
        case 1:
            ((RII_filter<FLOATING_TYPE, 1> *) filters[filter])->append(value);
            break;
        
        case 2:
            ((RII_filter<FLOATING_TYPE, 2> *) filters[filter])->append(value);
            break;
        
        case 3:
            ((RII_filter<FLOATING_TYPE, 3> *) filters[filter])->append(value);
            break;
        
        case 4:
            ((RII_filter<FLOATING_TYPE, 4> *) filters[filter])->append(value);
            break;
        
        case 5:
            ((RII_filter<FLOATING_TYPE, 5> *) filters[filter])->append(value);
            break;
        
        case 6:
            ((RII_filter<FLOATING_TYPE, 6> *) filters[filter])->append(value);
            break;
    }
}


FLOATING_TYPE filter_get_value(int filter) {
    int filterOrder = filtersOrders[filter];
    if (filter >= nbFilters) {
        std::cout << "Wrong filter given" << std::endl;
        exit(EXIT_FAILURE);
    }

    switch(filterOrder) {
        case 1:
            return ((RII_filter<FLOATING_TYPE, 1> *) filters[filter])->get_value();
        
        case 2:
            return ((RII_filter<FLOATING_TYPE, 2> *) filters[filter])->get_value();
        
        case 3:
            return ((RII_filter<FLOATING_TYPE, 3> *) filters[filter])->get_value();
        
        case 4:
            return ((RII_filter<FLOATING_TYPE, 4> *) filters[filter])->get_value();
        
        case 5:
            return ((RII_filter<FLOATING_TYPE, 5> *) filters[filter])->get_value();
        
        case 6:
            return ((RII_filter<FLOATING_TYPE, 6> *) filters[filter])->get_value();
    }
    std::cout << "Something went wrong\n" << std::endl;
    exit(EXIT_FAILURE);
}


FLOATING_TYPE filter_get_original_value(int filter) {
    int filterOrder = filtersOrders[filter];
    if (filter >= nbFilters) {
        std::cout << "Wrong filter given" << std::endl;
        exit(EXIT_FAILURE);
    }

    switch(filterOrder) {
        case 1:
            return ((RII_filter<FLOATING_TYPE, 1> *) filters[filter])->get_original_value();
        
        case 2:
            return ((RII_filter<FLOATING_TYPE, 2> *) filters[filter])->get_original_value();
        
        case 3:
            return ((RII_filter<FLOATING_TYPE, 3> *) filters[filter])->get_original_value();
        
        case 4:
            return ((RII_filter<FLOATING_TYPE, 4> *) filters[filter])->get_original_value();
        
        case 5:
            return ((RII_filter<FLOATING_TYPE, 5> *) filters[filter])->get_original_value();
        
        case 6:
            return ((RII_filter<FLOATING_TYPE, 6> *) filters[filter])->get_original_value();
    }
    std::cout << "Something went wrong\n" << std::endl;
    exit(EXIT_FAILURE);
}


void filter_destroy(int filter) {
    if (filter >= nbFilters) {
        std::cout << "Wrong filter given" << std::endl;
        exit(EXIT_FAILURE);
    }
    int filterOrder = filtersOrders[filter];
    int *filterTr = filters[filter];

    switch(filterOrder) {
        case 1:
            delete ((RII_filter<FLOATING_TYPE, 1> *) filterTr);
            break;
        
        case 2:
            delete ((RII_filter<FLOATING_TYPE, 2> *) filterTr);
            break;
        
        case 3:
            delete ((RII_filter<FLOATING_TYPE, 3> *) filterTr);
            break;
        
        case 4:
            delete ((RII_filter<FLOATING_TYPE, 4> *) filterTr);
            break;
        
        case 5:
            delete ((RII_filter<FLOATING_TYPE, 5> *) filterTr);
            break;
        
        case 6:
            delete ((RII_filter<FLOATING_TYPE, 6> *) filterTr);
            break;
    }
    for (int **curr = filters+filterOrder; curr < filters+nbFilters-1; curr++ ) {
        curr[0] = curr[1];
    }
    for (int *curr = filtersOrders+filterOrder; curr < filtersOrders+nbFilters-1; curr++ ) {
        curr[0] = curr[1];
    }
    nbFilters--;
    if (nbFilters == 0) {
        free(filters);
        free(filtersOrders);
    }
}