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

void create_random_array(star_t * stars, const int size)
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
  int i, j, a; 
	a = 0;
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
}


void fill_matravx(float_t *matrix, int size, float_t *xv, float_t * yv, float_t * zv, float_t * sf){
    int i,j,a;
    __m256 xi,yi,zi,xj,yj,zj,x1,y1,z1,x2,y2,z2,dist,dist2,dist3,sfi,sfj,sf1,sf2,sf3,sfm,sfr,res,cDiv;
    float con = 0.6;
		a = 0;
    cDiv = _mm256_set1_ps(con);
    for(i = 1 ; i < size-1; i++){
        xi = _mm256_set1_ps(xv[i]);
        yi = _mm256_set1_ps(yv[i]);
        zi = _mm256_set1_ps(zv[i]);
        sfi = _mm256_set1_ps(sf[i]);        
      for(j = i+1; j < size; j+=8){           
            xj = _mm256_loadu_ps(xv+j);
            yj = _mm256_loadu_ps(yv+j);
            zj = _mm256_loadu_ps(zv+j);
            sfj= _mm256_loadu_ps(sf+j);
            
            sf1= _mm256_mul_ps(sfi,sfj);
            sfm= _mm256_div_ps(sf1,cDiv);
            sf2= _mm256_add_ps(sfi,sfj);
            sf3 = _mm256_add_ps(sf2,sfm);
            sfr = _mm256_sqrt_ps(sf3);
            
            x1 = _mm256_sub_ps(xi,xj);
            y1 = _mm256_sub_ps(yi,yj);
            z1 = _mm256_sub_ps(zi,zj);
            x2 = _mm256_mul_ps(x1,x1);
            y2 = _mm256_mul_ps(y1,y1);
            z2 = _mm256_mul_ps(z1,z1);
            dist2 = _mm256_add_ps(x2,y2);
            dist3 = _mm256_add_ps(dist2,z2);
            dist = _mm256_sqrt_ps(dist3);
            
            res = _mm256_add_ps(dist,sfr);
            _mm256_storeu_ps(matrix+a+j-i,res);
        }
		 a+=size-i;
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
	int i,j, a, k, crd; //crd = coordinates
	float max = 0;
	float bin_size;
	float_t tmp;
	hist_param_t histparams;
	float_t *matrix_cpy = (float_t*)malloc((mat_size * mat_size / 2 + mat_size/2) * sizeof(float_t));
	k = mat_size - 1;
	a = mat_size - 1;
	memcpy(matrix_cpy, matrix, (mat_size * mat_size / 2 + mat_size/2) * sizeof(float_t));  
	
	//Need to do this mayn
	float min = 2 * abs(matrix_cpy[mat_size] - matrix_cpy[2]) +  //up
      	2 * abs(matrix_cpy[mat_size] - matrix_cpy[mat_size + 1]); //right
				 
	//Von Neumann calculations
    for(i = 1 ; i < mat_size-1; i++){
      for(j = i +1 ; j < mat_size ; j++){
      	crd = a+j-i; //Load in four "top" and four "bottom" cells to coalesce memory better
      	tmp = matrix[crd]; // and then perform vector addition
      	matrix[crd] = abs(tmp - matrix_cpy[crd - k+1 - (1/i)]) +  //up
      	 abs(tmp -matrix_cpy[crd + k - 2 - (2*k-3)*(i==j-1) + (1/i)]) + //down
      	 abs(tmp - matrix_cpy[crd - 1 + 2*(i==j-1)])  + //left
      	 abs(tmp - matrix_cpy[crd + 1]); //Input the vector opåerator for dis cpu?
      	matrix[crd] /= 4; //Input floating point operator << instead
      //printf("index is: %i", a+j-i );
      if(matrix[crd] > max) max = matrix[a+j-i];
      else if(matrix[crd] < min) min = matrix[a+j-i];
      //if(a+j-i==46)
      //printf("%.2e in place %i.. up:%i   down:%i   left:%i   right:%i \n", // <- proof of concept
			//matrix[19]  ,a+j-i, crd - k+1 - (1/i), crd + k - 2 - (2*k-3)*(i==j-1) + (1/i) ,a+j-i-1 + 2*(i==j-1), a+j-i+1);  //<<shorten these
    	}
			k = mat_size - i;
			a += mat_size - i;
		}
		
    bin_size = (max-min)/(hist_size);
		a = mat_size - 1;
    for(i = 1 ; i < mat_size-1; i++){
      for(j = i+1; j < mat_size; j++){
      	k = (int)((matrix[a+j-i]-min)/bin_size); //Do this for four number at a time
      	//printf("k is: %i  and index %i.  put in? %i\n",k, a+j-i , (i==j-1)); // <- proof of concept
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
