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

typedef struct star1{
  unsigned short subType;      // random: 0-9
  char spectralType;           // random: O, B, A, F, G, K, M, L, T
  char designation[8]; 	       // sprintf("%c%d.%d", spectralType, subType, index)
  float_t magnitude;  
	  struct pos1{
    float_t x, y, z;          // x & y random in (-1e5, 1e5), z random in (-3e3, 3e3)
  } position;	// random: (-10, +20)	
  int index;                   // counting index
} star_t1;

typedef struct star2{
	  struct pos2{
    float_t x, y, z;          // x & y random in (-1e5, 1e5), z random in (-3e3, 3e3)
  } position;
  int index;                   // counting index
  unsigned short subType;      // random: 0-9
  char spectralType;           // random: O, B, A, F, G, K, M, L, T
  char designation[9]; 	       // sprintf("%c%d.%d", spectralType, subType, index)
  float_t magnitude;           // random: (-10, +20)
} star_t2;

typedef struct star3{
	  struct pos3{
    float_t x, y, z;          // x & y random in (-1e5, 1e5), z random in (-3e3, 3e3)
  } position;
  char spectralType;           // random: O, B, A, F, G, K, M, L, T
  char designation[8]; 	       // sprintf("%c%d.%d", spectralType, subType, index)
  float_t magnitude;           // random: (-10, +20)7
	int index;                   // counting index
	unsigned short subType;      // random: 0-9
} star_t3;



typedef struct hist_params{
  int hist_size;
  float_t min, max, bin_size;
} hist_param_t;

#endif
