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
#include "ref_input.h"

void printtime(clock_t s, clock_t e)
{
  printf("time: %f\n", (double)(e-s)/CLOCKS_PER_SEC);    
}

int main(int argc, char **argv)
{
  clock_t start, end;
  if(argc < 2)
    {
      printf("usage: ./a.out N\n");
      return 0;
    }
  const int N = atoi(argv[1]);
  star_t *stars;
  stars = (star_t *) malloc(N*sizeof(star_t));
   
  printf("creating random stars: \t");
  start = clock(); 
   
  create_random_array(stars, N);
  //create_ref_star_array(stars, N);
   
  end = clock();
  printtime(start, end);
  //print_stars(stars, N);   // <--- This has output
         
  printf("sorting stars:    \t");
  start = clock();

  sort(stars, N);
   
  end = clock();
  printtime(start, end);
	
	//checksorted(stars,N);
  //print_stars(stars, N);   // <--- This has output
   
  printf("allocating matrix: \t");
  float_t *matrix;
  start = clock();
  if(N%2){
  	matrix = malloc((N * N / 2 + N/2 + 1) * sizeof(float_t));
  	}
  else{
  	matrix = malloc((N * N / 2 + N/2) * sizeof(float_t));
  	}
  end = clock();
  printtime(start, end);
   
  printf("filling matrix: \t");
  start = clock();
	fill_matrix(stars, matrix, N);
  end = clock();
  printtime(start, end);
	
	if(N<14){   //<- this is tight
		int i ,j, s;
		float_t *ptrs[N*N] ;
		s = 0;
		for ( i = 0; i < N; i++) {
			for (j = i; j < N ; j++){
				ptrs[i*N + j] = &matrix[s + j - i];
				ptrs[(i + j*N)] = &matrix[s + j - i]; 
				//if(s + j -i == 49)
					//printf("%.2e in place %i.. to add to %i and %i\n", // <- proof of concept
					//matrix[s + j - i]  ,s+j-i, i*N + j, (i + j*N)*(1 - i==j));
			}
				s +=  N-i; 
		}
		
	print_matrix(ptrs, N);   // <--- This has output if defined alse helps to prove a point
  }
  printf("generating histogram: \t");
  start = clock();
  int *histogram = (int *)calloc(NUM_HIST_BOXES+1,sizeof(int));
  hist_param_t histparams = generate_histogram(matrix, histogram, N, NUM_HIST_BOXES);
  end = clock();
  printtime(start, end);
	free(matrix);
	free(stars);
  display_histogram(histogram, histparams);
  free(histogram);
}
