/*****
 *  Code for examination project in High Performance Computing and Programming
 * 
 *  common.h common definitions header file
 *
 *  Author: Marcus Holm
 *  Modified by: Elias Rudberg
 *
 **/

#ifndef _COMMON_H
#define _COMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define float_t float 
#define NUM_HIST_BOXES 10

//changed this shit
typedef struct star{
	  struct pos{
    float_t x, y, z;          // x & y random in (-1e5, 1e5), z random in (-3e3, 3e3)
  } position;
  int index;                   // counting index
  unsigned short subType;      // random: 0-9
  char spectralType;           // random: O, B, A, F, G, K, M, L, T
  char designation[8]; 	       // sprintf("%c%d.%d", spectralType, subType, index)
  float_t magnitude;           // random: (-10, +20)
} star_t;

typedef struct hist_params{
  int hist_size;
  float_t min, max, bin_size;
} hist_param_t;

#endif
