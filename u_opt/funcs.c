/*****
 *  Code for examination project in High Performance Computing and Programming
 * 
 *  funcs.c functions implementation file
 *
 *  Author: Marcus Holm
 *  Modified by: Elias Rudberg
 *
 **/

#include "funcs.h"
#include "time.h"
#include "string.h"


void create_random_array(star_t * stars, int size)
{
	srand(time(NULL));
	int i;
	const char STypes[9] = {"OBAFGKMLT"};	 
	star_t* filler = malloc(sizeof(star_t)); 
	for(i = 0; i < size; i++){
		filler->index = i; 
		filler->subType = (unsigned short)(rand() % 10);
		filler->magnitude = 30*(float)rand()/(float)RAND_MAX-10;
		filler->spectralType = STypes[rand() % 9];
		filler->position.x = 200000*(float)rand()/(float)RAND_MAX-1000000;
		filler->position.y = 200000*(float)rand()/(float)RAND_MAX-1000000;
		filler->position.z = 6000*(float)rand()/(float)RAND_MAX-3000;
		stars[i] = filler[0];
	}
	free(filler);
}

void print_stars(star_t* array, int n)
{
  int i;
  printf("\nprint_stars, n = %d:\n", n);
  for(i = 0; i<n; i++)
    printf("%s ",array[i].designation);
  printf("\n");
}


float_t starfunc(star_t a, star_t b)
{
  unsigned short x = a.subType;
  unsigned short y = b.subType;
  return sqrt(x + y + x*y/0.6);
}


void sort(star_t* array, int n)  //bubble sort
{
	  int i, j;
  for(i = 0; i < n-1; i++){
  	float d1 = sqrt(array[i].position.x*array[i].position.x + array[i].position.y*array[i].position.y);
    for(j = 1+i; j < n; j++) {
    float d2 = sqrt(array[j].position.x*array[j].position.x + array										[j].position.y*array[j].position.y);
      if(d1 > d2) {
				star_t tmp = array[i];
				array[i] = array[j];
				array[j] = tmp;
      }
    }
  }
}

void fill_matrix(star_t * array, float_t **matrix, int size)
{
  int i, j; 
  for(i = 0 ; i < size; i++)
      for(j = 0 ; j < size ; j++){
			float d2 = sqrt(pow(array[j].position.x - array[i].position.x,2) + pow(array[i].position.y - array[j].position.y,2));
			matrix[i][j] = d2 + starfunc(array[j],array[i]);
    }
}

void print_matrix(float_t** theMatrix, int n)
{
  int i, j;
  printf("\nprint_matrix, n = %d:\n", n);
  for(i = 0 ; i < n; i++)
    {
      for(j = 0 ; j < n ; j++)
	printf("%.2f " , theMatrix[i][j]);
      putchar('\n');
    }
}

hist_param_t generate_histogram(float_t **matrix, int *histogram, int mat_size, int hist_size)
{
	int i,j;
	float min = 10000;
	float max = 0;
	float tmp, bin_size;
	hist_param_t histparams;
	float_t **matrix_cpy = (float_t**)malloc(mat_size *sizeof(float_t*));
	for(i = 0; i < mat_size; i++) 
		matrix_cpy[i] = (float_t *)malloc(mat_size * sizeof(float_t));
		
	for(i = 0; i < mat_size; i++)
    memcpy(&matrix_cpy[i][0], &matrix[i][0], mat_size * sizeof(float_t));
    
  for(i = 1 ; i < mat_size-1; i++)
      for(j = 1 ; j < mat_size-1 ; j++){
      	tmp = matrix[i][j];
      	matrix[i][j] = abs(tmp - matrix_cpy[i-1][j]) + abs(tmp -matrix_cpy[i][j-1]) + abs(tmp - matrix_cpy[i+1][j])  + abs(tmp - matrix_cpy[i][j+1]); //Input the vector opåerator for dis cpu?
      	matrix[i][j] /= 4; //Input floating point operator << instead
      if(matrix[i][j] > max) max = matrix[i][j];
      else if(matrix[i][j] < min) min = matrix[i][j];
    }
  bin_size = (max-min)/9;
  for(i = 1 ; i < mat_size-1; i++)
      for(j = 1; j < mat_size -1; j++){
      	int a = (int)(matrix[i][j]/(bin_size+2));
  			histogram[a] += 1;
  			if(a > 8)
  				printf("\n argument: %i\n",a);
  		}
  printf("\n argument: %i and before %f\n",(int)(bin_size), bin_size);
  for(i = 0; i < mat_size; i++)
  	free(matrix_cpy[i]);
	free(matrix_cpy);
	histparams.hist_size = hist_size;
	histparams.bin_size = bin_size;
	histparams.min = min;
	histparams.max = max;
	return histparams;
}

void display_histogram(int *histogram, hist_param_t histparams)
{
  printf("\ndisplay_histogram:\n");
  int i,j;
  for(i = 0; i < histparams.hist_size && histparams.bin_size*i < histparams.max; i++)
    {
      printf("%11.3e ", histparams.bin_size * i + histparams.min);
    }
  printf("%11.3e\n", histparams.max);
  for(j = 0; j < i; j++)
    {
      printf("%11d ", histogram[j]);
    }
  printf("\n");
}
