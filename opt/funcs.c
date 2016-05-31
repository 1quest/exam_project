/*****
 *  
 *  DISCLAIMER: MY INTENTION WAS TO IMPLEMENT SIMD BUT MY PUNY BRAIN 
 *  
 *
 **/

#include "funcs.h"
#include "time.h"
#include "string.h"
#include <immintrin.h>
#define vec_size 8

void create_random_array(star_t * stars, const int size)
{
	srand(time(NULL));
	int i;
	const char STypes[] = {"OBAFGKMLT"};
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
//Change this to take three float_t
static inline float_t L2norm(star_t star){
    return sqrt(star.position.x*star.position.x+star.position.y*star.position.y+star.position.z*star.position.z);
}

void print_stars(star_t* array, const int n)
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
//Change this to ushort in so that inline becomes useful
static inline float_t starfunc(star_t a, star_t b)
{
  unsigned short x = a.subType;
  unsigned short y = b.subType;
  return sqrt(x + y + x*y/0.6);
}

//Do something supersmart about the d1&d2 if you have time
void sort(star_t* array, const int n)
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

int checksorted(star_t* array, const int n)
{
	int i;
	float d1,d2;
	for (i = 0; i < n -1; i++) {
		d1 = sqrt(array[i].position.x*array[i].position.x + array[i].position.y*array[i].position.y);
		d2 = sqrt(array[i+1].position.x*array[i+1].position.x + array[i+1].position.y*array[i+1].position.y);
		if( d1 > d2 ){
			printf("The list hasnt been sorted : i is %i\n", i);
			return 0;
		}
	}
	printf("List was sorted\n");
	return 1;
}

void fill_matrix(star_t * array, float_t *matrix, const int size)
{
  int i, j, a, z; 
	a = 0;
	z = 0;
	//putchar('\n'); //For proving its correct
	//float_t * temp = (float_t*)malloc(size * sizeof(float_t));
  for(i = 0 ; i < size; i++){
		float_t x1 = array[i].position.x;
		float_t y1 = array[i].position.y; 
		float_t z1 = array[i].position.z;
    for(j = i; j < size; j++){ //Add vector addition
			float_t d2 = sqrt(pow(array[j].position.x - x1,2) + pow(y1 - array[j].position.y,2) + pow(z1 - array[j].position.z,2));
			matrix[z] = (float_t)(d2 + starfunc(array[j],array[i]));
			//printf("%.2e ",matrix[a+j-i]); //For proving its correct
			//printf("%.2i in place %i.. which becomes n%i%i\n", z, a+j-i, i,j);  //<<shorten these
			z++;
    }
		//putchar('\n');  //For proving its correct
		a+=size-i;
	}

}

//Hit and miss
void fill_matrixav(star_t * array, float_t *matrix, const int size)
{
  int i;
	float_t *x1 = malloc(size * sizeof(float_t));
	float_t *y1 = malloc(size * sizeof(float_t));
	float_t *z1 = malloc(size * sizeof(float_t));
	float_t *sf = malloc(size * sizeof(float_t));
  for(i = 0 ; i < size; i++){
		x1[i] = array[i].position.x;
		y1[i] = array[i].position.y; 
		z1[i] = array[i].position.z;
		sf[i] = array[i].subType;
	}
	fill_matravx(matrix, size, x1, y1, z1, sf);
	free(x1);
	free(y1);
	free(z1);
	free(sf);
}

//Hit and miss
void fill_matravx(float_t *matrix, int size, float_t *xv, float_t * yv, float_t * zv, float_t * sf){
    int i,j,a;
    __m256 xi, yi, zi, xj, yj, zj, dist, sfi, sfj, sf2, sfm, MUL;
		a = 0;
		float c = 0.6;
    MUL = _mm256_set1_ps(c);
    for(i = 0 ; i < size; i++){
        xi = _mm256_set1_ps(xv[i]);
        yi = _mm256_set1_ps(yv[i]);
        zi = _mm256_set1_ps(zv[i]);
        sfi = _mm256_set1_ps(sf[i]);        
      for(j = i; j < size-vec_size; j+=vec_size){           
            xj = _mm256_loadu_ps(xv+j);
            yj = _mm256_loadu_ps(yv+j);
            zj = _mm256_loadu_ps(zv+j);
            sfj= _mm256_loadu_ps(sf+j);
            
						//Starfunc implementation on vector-form
            sfj= _mm256_mul_ps(sfi,sfj);
            sfm= _mm256_div_ps(sfj,MUL);
            sf2= _mm256_add_ps(sfi,sfj);
            sfj = _mm256_add_ps(sf2,sfm);
            sfj = _mm256_sqrt_ps(sfj);
            
            xj = _mm256_sub_ps(xi,xj);
            yj = _mm256_sub_ps(yi,yj);
            zj = _mm256_sub_ps(zi,zj);
            xj = _mm256_mul_ps(xj,xj);
            yj = _mm256_mul_ps(yj,yj);
            zj = _mm256_mul_ps(zj,zj);
            dist = _mm256_add_ps(xj,yj);
            dist = _mm256_add_ps(dist,zj);
            dist = _mm256_sqrt_ps(dist);
            
            dist = _mm256_add_ps(dist,sfj);
            _mm256_storeu_ps(matrix+a,dist);						
						
        }
        a=a+8;
        printf("place: %i", a);
				printf("\n");
    }
}


void print_matrix(float_t** theMatrix, const int n)
{
  int i, j;
  printf("\nprint_matrix, n = %d:\n", n);
  for(i = 0 ; i < n; i++)
    {
      for(j = 0 ; j < n ; j++)
			printf("%.2f " , *theMatrix[i*n+j]); //.2f
      putchar('\n');
    }
}

hist_param_t generate_histogram(float_t *matrix, int *histogram, const int mat_size, int hist_size)
{
	int i,j, a, k, crd, z, N; //crd = coordinates
	float max = 0;
	float bin_size;
	float_t tmp;
	hist_param_t histparams;
	float_t *matrix_cpy;
	if(mat_size%2){
  	N = ((mat_size-2) * (mat_size-2) / 2 + (mat_size-2)/2 + 1);
  	}
  else{
  	N = ((mat_size-2) * (mat_size-2) / 2 + (mat_size-2)/2);
  	}
  matrix_cpy = (float_t*)malloc(N * sizeof(float_t));
	k = mat_size - 1;
	a = mat_size - 1;
	z = 0;
	
	//Need to do this mayn
	float min = 2 * abs(matrix[mat_size] - matrix[2]) +  //up
      	2 * abs(matrix[mat_size] - matrix[mat_size + 1]); //right
				 
	//Von Neumann calculations
    for(i = 1 ; i < mat_size-1; i++){			
      for(j = i + 1 ; j < mat_size ; j++){
      	crd = a+j-i; //Load in four "top" and four "bottom" cells to coalesce memory better
      	tmp = matrix[crd]; // and then perform vector addition
      	matrix_cpy[z] = abs(tmp - matrix[crd - k+1 - (1/i)]) +  //up
      	 abs(tmp -matrix[crd + k - 2 - (2*k-3)*(i==j-1) + (1/i)] - 2*(z==0)) + //down
      	 abs(tmp - matrix[crd - 1 + 2*(i==j-1)])  + //left
      	 abs(tmp - matrix[crd + 1]); //Input the vector opåerator for dis cpu?
      	 matrix_cpy[z] /= 4; //Input floating point operator << instead
				//printf("index is: %i", a+j-i );
				if(matrix_cpy[z] > max) max = matrix_cpy[z];
				else if(matrix_cpy[z] < min) min = matrix_cpy[z];
				//if(a+j-i==46)
				//printf("%.2e in place %i.. up:%i   down:%i   left:%i   right:%i \n", // <- proof of concept
				//matrix[19]  ,a+j-i, crd - k+1 - (1/i), crd + k - 2 - (2*k-3)*(i==j-1) + (1/i) - 2*(z==0),a+j-i-1 + 2*(i==j-1), crd+1); //<<shorten these
    		z++;
			}
			k = mat_size - i;
			a += mat_size - i;
		}
		//printf("WHATS UP BOSS: %i, N IS:%i\n", z , N);
    bin_size = (max-min)/(hist_size);
		a = 0;
		z = 0;
    for(i = 0 ; i < N; i++){
      	k = (int)((matrix_cpy[i]-min)/bin_size); //Do this for four number at a time
      	//printf("k is: %i  and index %i.  put in? %i\n",k, a+j-i , (i==j-1)); // <- proof of concept
  			histogram[k] += 1 + (i!=a);
				z += (i==a);
				a += (mat_size - 2  - z) * (i==a);
		}
		//printf("WHATS UP BOSS: %i\n", z);
		histogram[9] += histogram[10];
	histparams.hist_size = hist_size;
	histparams.bin_size = bin_size;
	histparams.min = min;
	histparams.max = max;
	free(matrix_cpy);
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
