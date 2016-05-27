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
#define QS

void create_random_array(star_t * stars, int size)
{
	srand(time(NULL));
	int i;
	const char STypes[9] = {"OBAFGKMLT"};
	for(i = 0; i < size; i++){
		stars[i].index = i; 
		stars[i].subType = (unsigned short)(rand() % 10);
		stars[i].magnitude = 30*(float)rand()/(float)RAND_MAX-10;
		stars[i].spectralType = STypes[rand() % 9];
		stars[i].position.x = 200000*(float)rand()/(float)RAND_MAX-1000000;
		stars[i].position.y = 200000*(float)rand()/(float)RAND_MAX-1000000;
		stars[i].position.z = 6000*(float)rand()/(float)RAND_MAX-3000;
		}
}

static inline float_t L2norm(star_t star){
    return sqrt(star.position.x*star.position.x+star.position.y*star.position.y+star.position.z*star.position.z);
}

void print_stars(star_t* array, int n)
{
    int i;
    printf("\nprint_stars, n = %d:\n", n);
    for(i = 0; i<n; i++)
        //printf("%s ",array[i].designation);
    printf("\n");
    for (i=0; i<n; i++)
        //printf("%f ",L2norm(array[i]));
    printf("\n");
}

static inline float_t starfunc(star_t a, star_t b)
{
  unsigned short x = a.subType;
  unsigned short y = b.subType;
  return sqrt(x + y + x*y/0.6);
}

void sort(star_t* array, int n)  //bubble sort
{
	  int i, j;
    float p, d1, d2;
		star_t t;
    if (n < 2)
        return;
    p = sqrt(array[n/2].position.x*array[n/2].position.x + array[n/2].position.y*array[n/2].position.y);
    for (i = 0, j = n - 1;; i++, j--) {
		d1 = sqrt(array[i].position.x*array[i].position.x + array[i].position.y*array[i].position.y);
		d2 = sqrt(array[j].position.x*array[j].position.x + array[j].position.y*array[j].position.y);
        while (d1 < p){
					i++;
					d1 = sqrt(array[i].position.x*array[i].position.x + array[i].position.y*array[i].position.y); 
				}
        while (p < d2){
          j--;
					d2 = sqrt(array[j].position.x*array[j].position.x + array[j].position.y*array[j].position.y);
				}
        if (i >= j)
            break;
        t = array[i];
        array[i] = array[j];
        array[j] = t;
    }
    sort(array, i);
    sort(array + i, n - i);
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
    bin_size = (max-min)/10.0;
    for(i = 1 ; i < mat_size-1; i++)
      for(j = 1; j < mat_size -1; j++){
      	int a = (int)(matrix[i][j]/(bin_size));
  			histogram[a] += 1;
  		}
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
