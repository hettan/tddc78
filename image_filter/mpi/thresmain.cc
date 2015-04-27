//#include <stdlib.h>
//#include <stdio.h>
//#include <string.h>
//#include <mpi.h>
#include <iostream>
#include <time.h>
//#include "ppmio.h"

#define MAX_PIXELS (3000*3000)
using namespace std;

struct pixel
{
  unsigned char r,g,b;
};

/*
void init_pixel_type(MPI_Datatype& pixel_type)
{
  const int struct_size = 3;
  int blockcounts[] = {1, 1, 1};
  MPI_Aint disp[] = {0, 0, 0};
  MPI_Datatype old_types[struct_size] = {MPI::UNSIGNED_CHAR,
					 MPI::UNSIGNED_CHAR,
					 MPI::UNSIGNED_CHAR};

  MPI_Type_struct(struct_size, blockcounts, disp, old_types, &pixel_type);
  MPI_Type_commit(&pixel_type);
  }*/


void thresfilter(const int xsize, const int ysize, pixel* src){
#define uint unsigned int 
  
  uint sum, i, psum, nump;

  nump = xsize * ysize;

  for(i = 0, sum = 0; i < nump; i++) {
    sum += (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
  }

  sum /= nump;

  for(i = 0; i < nump; i++) {
    psum = (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
    if(sum > psum) {
      src[i].r = src[i].g = src[i].b = 0;
    }
    else {
      src[i].r = src[i].g = src[i].b = 255;
    }
  }
}

int main (int argc, char ** argv) {
  cout << "here" << endl;
  int xsize, ysize, colmax;
  pixel src[MAX_PIXELS];
  struct timespec stime, etime;

  /* Take care of the arguments */
  cout << "hllo" << endl;
  if (argc != 3) {
    cout << "Usage: " << argv[0] << " infile outfile" << endl;
    return 1;
  }

  /* read file */
  /*  if(read_ppm (argv[1], &xsize, &ysize, &colmax, (char *) src) != 0)
    return 1;
  */
  if (colmax > 255) {
    cout << "Too large maximum color-component value" << endl;
    return 1;
  }

  cout << "Has read the image, calling filter" << endl;

  clock_gettime(CLOCK_REALTIME, &stime);

  thresfilter(xsize, ysize, src);

  clock_gettime(CLOCK_REALTIME, &etime);

  const double delta_time = (etime.tv_sec  - stime.tv_sec) + 1e-9*(etime.tv_nsec  - stime.tv_nsec);
  cout << "Filtering took: " << delta_time <<  " secs" << endl  ;

  /* write result */
  cout << "Writing output file" << endl;
    
  /*  if(write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
   return 1;
  */

  return 0;
}
