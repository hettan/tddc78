#include <math.h>
#include <mpi.h>

#define MAX_X 1.33
#define Pi 3.14159

void get_gauss_weights(int n, double* weights_out)
{
  double x;
  int i;
  for( i=0; i<n+1; i++ ){
    x=(double)i * MAX_X/n;
    weights_out[i] ) exp(-x*x * Pi);
  }
}
