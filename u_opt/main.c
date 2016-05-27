/*****
 *  Code for examination project in High Performance Computing and Programming
 * 
 *  main.c main implementation file 
 *  
 *  Author: Marcus Holm
 *  Modified by: Elias Rudberg
 *
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "funcs.h"

void printtime(clock_t s, clock_t e)
{
  printf("time: %f\n", (double)(e-s)/CLOCKS_PER_SEC);    
}

int main(int argc, char **argv)
{
  int N, i;
  clock_t start, end;
  if(argc < 2)
    {
      printf("usage: ./a.out N\n");
      return 0;
    }
  N = atoi(argv[1]);
  star_t *stars;
  stars = (star_t *) malloc(N*sizeof(star_t));
   
  printf("creating random stars: \t");
  start = clock();
   
  create_random_array(stars, N);
   
  end = clock();
  printtime(start, end);
  // print_stars(stars, N);
         
  printf("sorting stars:    \t");
  start = clock();

  sort(stars, N);
   
  end = clock();
  printtime(start, end);
  // print_stars(stars, N);
   
  printf("allocating matrix: \t");
  start = clock();
  float_t **matrix = (float_t**)malloc(N * sizeof(float_t*));
	for(i = 0; i < N; i++) 
		matrix[i] = (float_t *)malloc(N * sizeof(float_t));
   
  end = clock();
  printtime(start, end);
   
  printf("filling matrix: \t");
  start = clock();
  fill_matrix(stars, matrix, N);
   
  end = clock();
  printtime(start, end);
  // print_matrix(matrix, N);
   
  printf("generating histogram: \t");
  start = clock();
  int *histogram = (int *)calloc(NUM_HIST_BOXES,sizeof(int));
  hist_param_t histparams = generate_histogram(matrix, histogram, N, NUM_HIST_BOXES);
  end = clock();
  printtime(start, end);
	for(i = 0; i < N; i++) 
		free(matrix[i]);
	free(matrix);
  display_histogram(histogram, histparams);
}
