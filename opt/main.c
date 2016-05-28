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
  int N, i ,j, k, s;
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
  //print_stars(stars, N);   // <--- This has output
         
  printf("sorting stars:    \t");
  start = clock();

  sort(stars, N);
   
  end = clock();
  printtime(start, end);
  //print_stars(stars, N);   // <--- This has output
   
  printf("allocating matrix: \t");
  start = clock();
  float_t *matrix = malloc((N * N / 2 + N/2) * sizeof(float_t));
  end = clock();
  printtime(start, end);
   
  printf("filling matrix: \t");
  start = clock();
  fill_matrix(stars, matrix, N);
	printf("%.2e\n",matrix[19]);
	float_t **ptrs[N*N];
	k = 0;
	s = 0;
	if(N<500)
		for ( i = 0; i < N; i++) {
			//ptrs[i*N + i] = &matrix[s];
			printf("%i   ---   %i\n", i + i*N, s);
			for (j = i; j < N ; j++){
				printf("%i   ---   %i\n", i*N+ j, s+j-i);
				printf("%i   ---   %i\n", i + j*N, s+j-i);
				ptrs[i*N + j] = &matrix[s + j - i];
				ptrs[i + j*N] = &matrix[s + j - i]; 
			}
				k ++;
				s +=  N-i; 
		}
  end = clock();
  printtime(start, end);
	if(N < 20){
		print_matrix(ptrs, N);  // <--- This has output
	}
  printf("generating histogram: \t");
  start = clock();
  //int *histogram = (int *)calloc(NUM_HIST_BOXES+1,sizeof(int));
  //hist_param_t histparams = generate_histogram(matrix, histogram, N, NUM_HIST_BOXES);
  end = clock();
  printtime(start, end);
	free(matrix);
	free(stars);
  //display_histogram(histogram, histparams);
  //free(histogram);
}
