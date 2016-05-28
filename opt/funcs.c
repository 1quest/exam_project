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
		sprintf(stars[i].designation,"%c%d.%d", stars[i].spectralType, stars[i].subType, i);
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
        printf("%s ",array[i].designation);
    printf("\n");
		
		for(i = 0; i<n; i++)
        printf("%f ",L2norm(array[i]));
    printf("\n");
}

static inline float_t starfunc(star_t a, star_t b)
{
  unsigned short x = a.subType;
  unsigned short y = b.subType;
  return sqrt(x + y + x*y/0.6);
}

void sort(star_t* array, int n)
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

void fill_matrix(star_t * array, float_t *matrix, int size)
{
  int i, j, a; 
	a = 0;
	//putchar('\n'); //For proving its correct
	//float_t * temp = (float_t*)malloc(size * sizeof(float_t));
  for(i = 0 ; i < size; i++){
		float x1 = array[i].position.x;
		float y1 = array[i].position.y; 
    for(j = i; j < size; j++){
			float d2 = sqrt(pow(array[j].position.x - x1,2) + pow(y1 - array[j].position.y,2));
			matrix[a+j-i] = (float_t)(d2 + starfunc(array[j],array[i]));
			//printf("%.2e ",matrix[a+j-i]); //For proving its correct
    }
		//putchar('\n');  //For proving its correct
		a+=size-i;
	}
}

void print_matrix(float_t** theMatrix, int n)
{
  int i, j;
  printf("\nprint_matrix, n = %d:\n", n);
  for(i = 0 ; i < n; i++)
    {
      for(j = 0 ; j < n ; j++)
			printf("%.2e " , *theMatrix[i*n+j]); //.2f
      putchar('\n');
    }
}

hist_param_t generate_histogram(float_t *matrix, int *histogram, int mat_size, int hist_size)
{
	int i,j, a, k, crd; //crd = coordinates
	float min = 10000;
	float max = 0;
	float bin_size;
	float_t tmp;
	hist_param_t histparams;
	float_t *matrix_cpy = (float_t*)malloc((mat_size * mat_size / 2 + mat_size/2) * sizeof(float_t));
	k = mat_size - 1;
	a = mat_size - 1;
	memcpy(matrix_cpy, matrix, (mat_size * mat_size / 2 + mat_size/2) * sizeof(float_t));  
    for(i = 1 ; i < mat_size-1; i++){
      for(j = i ; j < mat_size ; j++){
      	crd = a+j-i;
      	tmp = matrix[crd];
      	matrix[crd] = abs(tmp - matrix_cpy[crd - k+1]) +  //up
      	 abs(tmp -matrix_cpy[crd + k - 2 - (2*k-3)*(i==j-1) - (1/i)]) + //down
      	 abs(tmp - matrix_cpy[crd - 1 + 2*(i==j-1)])  + //left
      	 abs(tmp - matrix_cpy[crd + 1]); //Input the vector opåerator for dis cpu?
      	matrix[crd] /= 4; //Input floating point operator << instead
      	//printf("index is: %i", a+j-i );
      if(matrix[crd] > max) max = matrix[a+j-i];
      else if(matrix[crd] < min) min = matrix[a+j-i];
      //if(a+j-i==46)
      	//printf("%.2e in place %i.. up:%i   down:%i   left:%i   right:%i \n", // <- proof of concept
			//matrix[19]  ,a+j-i, a+j-i - k+1 - (1/i)  ,  a+j-i + k-2  - (2*k-3)*(i==j-1) - (1/i) ,a+j-i-1 + 2*(i==j-1), a+j-i+1);  //<<shorten these
    	}
			k = mat_size - i;
			a += mat_size - i;
		}
    bin_size = (max-min)/(hist_size);
		a = mat_size - 1;
    for(i = 1 ; i < mat_size-1; i++){
      for(j = i+1; j < mat_size; j++){
      	k = (int)((matrix[a+j-i]-min)/bin_size); //Do this for four number at a time
      	printf("k is: %i  and index %i.  put in? %i\n",k, a+j-i , (i==j-1)); // <- proof of concept
  			histogram[k] += 2 - (i==(j-1));
  		}
			a+= mat_size - i ;
		}
		histogram[9] += histogram[10];
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
