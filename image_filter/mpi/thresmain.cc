#include <mpi.h>
#include <iostream>
#include <time.h>
#include "ppmio.h"

using namespace std;

struct pixel
{
  unsigned char r,g,b;
};

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
}

//Calculate the sum of all pixels in src
uint get_sum(const int xsize, const int ysize, pixel* src){
#define uint unsigned int
  uint sum, i, nump;
  
  nump = xsize * ysize;
  
  for(i = 0, sum = 0; i < nump; i++)
    sum += (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
   
  return sum;
}
  
void thresfilter(const int xsize, const int ysize, pixel* src, const int sum){  
#define uint unsigned int
  uint psum, nump, i;
  nump = xsize * ysize;

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
  
  MPI::Status status;
  MPI_Init( &argc, &argv );
  const MPI_Comm com = MPI_COMM_WORLD;
  const int root = 0;
  const int myid = MPI::COMM_WORLD.Get_rank();
  const int num_proc = MPI::COMM_WORLD.Get_size();

  MPI_Datatype pixel_type;
  init_pixel_type(pixel_type);
  const int pixel_type_size = 3;

  int xsize, ysize, colmax;
  pixel *src = nullptr;
  if(myid == root){
    src = new pixel[MAX_PIXELS];
    /* Take care of the arguments */
    if (argc != 3) {
      cout << "Usage: " << argv[0] << " infile outfile" << endl;
      return 1;
    }

    /* read file */
    if(read_ppm (argv[1], &xsize, &ysize, &colmax, (char *) src) != 0)
      return 1;
  
    if (colmax > 255) {
      cout << "Too large maximum color-component value" << endl;
      return 1;
    }
    cout << "Has read the image, calling filter" << endl;
  }

  struct timespec stime, etime;
  if(myid == root)
    clock_gettime(CLOCK_REALTIME, &stime);

  //Broadcast the variables from root to the other processes
  MPI_Bcast(&xsize, 1, MPI_INT, root, com);
  MPI_Bcast(&ysize, 1, MPI_INT, root, com);

  int local_y = ysize/num_proc;
  //The last process gets the leftover rows
  if(myid == num_proc -1)
    local_y += ysize % num_proc;

  pixel *local_src = new pixel[local_y*xsize];

  int send_count[num_proc];
  int disp[num_proc];
  int disp_size = 0;
  for (int i=0; i<num_proc; i++){
    if(i == num_proc -1)
      send_count[i] = ((ysize/num_proc) + (ysize%num_proc))*(xsize*pixel_type_size);
    else
      send_count[i] = (ysize/num_proc)*xsize*pixel_type_size;
    
    disp[i] = disp_size;
    disp_size += send_count[i];
  }
    
  MPI_Scatterv(src, send_count, disp, pixel_type,
	       local_src, send_count[myid], pixel_type, root, com);
    
  //Calculate the sum
  const double local_sum = (double)get_sum(xsize, local_y, local_src);
  double sum = 0;
  MPI_Allreduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, com);
  sum = sum/(xsize*ysize);

  thresfilter(xsize, local_y, local_src, sum);

  MPI_Gatherv(local_src, send_count[myid], pixel_type, src,
	      send_count, disp, pixel_type, root, com);
  
  //print results
  if(myid == root){
    clock_gettime(CLOCK_REALTIME, &etime);
    
    const double delta_time = (etime.tv_sec  - stime.tv_sec) + 1e-9*(etime.tv_nsec  - stime.tv_nsec);
    cout << "Filtering took: " << delta_time <<  " secs" << endl  ;

    /* write result */
    cout << "Writing output file" << endl;
    
    if(write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
      return 1;
  }

  //cleanup
  delete[] src;
  delete[] local_src;

  MPI_Finalize();
  return 0;
}
