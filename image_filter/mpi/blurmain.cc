#include <mpi.h>
#include <iostream>
#include <math.h>
#include <time.h>

#include "ppmio.h"

#define MAX_RAD 1000
#define MAX_X 1.33
#define Pi 3.14159
#define MAX_PIXELS (1000*1000)

using namespace std;

struct pixel
{
  unsigned char r,g,b;
};

namespace std{
  ostream& operator<<(ostream& os, pixel& p){
    os << "r=" << p.r << ", g=" << p.g << ", b=" << p.b;
    return os;
  }
}

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

/* Requires that (num_proc%n) == 0. */
void get_gauss_weights(const int n, double* work){
  const MPI_Comm com = MPI_COMM_WORLD;
  const int root = 0; //root process
  const int myid = MPI::COMM_WORLD.Get_rank();
  const int num_proc = MPI::COMM_WORLD.Get_size();
  
  int work_size_proc = n / num_proc;
  //Last process gets the last indexes as well
  if(myid == num_proc -1)
    work_size_proc += n % num_proc;

  //Create work buffer for each process
  double work_proc[const_cast<const int&>(work_size_proc)];
  const int start_index = myid*work_size_proc;
  double x;

  for(int i=start_index; i<start_index+work_size_proc; i++){    
    x=(double)i * MAX_X/10;
    work_proc[i-start_index] = exp(-x*x * Pi);
  }

  int recv_count[num_proc];
  int displs[num_proc];

  //Specify how much data will be recv for each process
  for(int i=0; i<num_proc; i++){
    recv_count[i] = n/num_proc;
    if(i == num_proc-1)
      recv_count[i] += n%num_proc;
    displs[i] = i * (n/num_proc);
  }
  
  MPI_Allgatherv(&work_proc, work_size_proc, MPI_DOUBLE, work,
  		 recv_count, displs, MPI_DOUBLE, com);
  
}

pixel* pix(pixel* image, const int xx, const int yy, const int xsize)
{
  register int off = xsize*yy + xx;

#ifdef DBG
  if(off >= MAX_PIXELS) {
    fprintf(stderr, "\n Terribly wrong: %d %d %d\n",xx,yy,xsize);
  }
#endif
  return (image + off);
}


void blurfilter(const int xsize, const int ysize, 
		pixel* src, const int radius, const double *w)
{
  MPI_Status status;
  const MPI_Comm com = MPI_COMM_WORLD;
  const int root = 0; //root process
  const int myid = MPI::COMM_WORLD.Get_rank();  
  const int num_proc = MPI::COMM_WORLD.Get_size();

  //pixeltype
  MPI_Datatype pixel_type;
  init_pixel_type(pixel_type);
  const int pixel_type_size = 3;

  int x,y,x2,y2, wi;
  double r,g,b,n, wc;
  
  int row_interval = ysize/num_proc;
  //The last process gets the leftover rows
  if(myid == num_proc -1)
    row_interval += ysize % num_proc;
  
  const int row_start = row_interval*myid;
  const int local_size = row_interval*xsize;
  const int local_radius_size = radius*xsize;

  pixel local_src[local_size];

  //Allocate mem for the radius as well.
  pixel local_dst[local_size+(2*local_radius_size)];
  
  MPI_Scatter(src, local_size*pixel_type_size, pixel_type, local_src,
	      local_size*pixel_type_size, pixel_type, root, com);

  if( (myid == root) && (ysize % num_proc) != 0 ){
    const int leftover_size = (ysize % num_proc)*xsize;
    const int leftover_offset = ysize*xsize - leftover_size;
    
    MPI_Send(src+leftover_offset, leftover_size*pixel_type_size,
	     pixel_type, num_proc-1, 0, com);
  }

  else if( (myid == num_proc-1) && (ysize % num_proc) != 0){
    const int leftover_size = (ysize % num_proc)*xsize;

    MPI_Recv(local_src+local_size-leftover_size, leftover_size*pixel_type_size,
	     pixel_type, 0, 0, com, &status);
  }    
  
  for (y=0; y<row_interval; y++) {
    for (x=0; x<xsize; x++) {
      r = w[0] * pix(local_src, x, y, xsize)->r;
      g = w[0] * pix(local_src, x, y, xsize)->g;
      b = w[0] * pix(local_src, x, y, xsize)->b;
      n = w[0];
      for ( wi=1; wi <= radius; wi++) {
	wc = w[wi];
	x2 = x - wi;
	if(x2 >= 0) {
	  r += wc * pix(local_src, x2, y, xsize)->r;
	  g += wc * pix(local_src, x2, y, xsize)->g;
	  b += wc * pix(local_src, x2, y, xsize)->b;
	  n += wc;
	}
	x2 = x + wi;
	if(x2 < xsize) {
	  r += wc * pix(local_src, x2, y, xsize)->r;
	  g += wc * pix(local_src, x2, y, xsize)->g;
	  b += wc * pix(local_src, x2, y, xsize)->b;
	  n += wc;
	}
      }
      pix(local_dst+local_radius_size,x,y, xsize)->r = r/n;
      pix(local_dst+local_radius_size,x,y, xsize)->g = g/n;
      pix(local_dst+local_radius_size,x,y, xsize)->b = b/n;
    }
  }
			   

  /* Send the required pixels needed for blur to neighbour processes */
    const int data_size = pixel_type_size*local_radius_size;
  if(myid != root){

    //Send radius pixels needed to the lower process.
    MPI_Send(local_dst+local_radius_size, data_size,
	     pixel_type, myid-1, 1, com);
    
    //Recv low index rows needed for y-blur
    MPI_Recv(local_dst, data_size,
	     pixel_type, myid-1, 2, com, &status);
  }
  
  if(myid != num_proc-1){

    //Recv higher index rows needed for y-blur
    MPI_Recv(local_dst+local_radius_size+(row_interval*xsize),
	     data_size, pixel_type, myid+1, 1, com, &status);

    //Send radius pixels needed to the higher process.
    MPI_Send(local_dst+(row_interval*xsize), data_size,
	     pixel_type, myid+1, 2, com);
  }
  
  int min, max;
  if(myid == 0){
    min=radius;
    max=row_interval + (2*radius);
  }
  else if(myid == num_proc-1){
    min=0;
    max=row_interval+radius;
  }
  else{
    min=0;
    max=row_interval + (2*radius);
  }
  
  for (y=radius; y<row_interval+radius; y++) {
    for (x=0; x<xsize; x++) {
      r = w[0] * pix(local_dst, x, y, xsize)->r;
      g = w[0] * pix(local_dst, x, y, xsize)->g;
      b = w[0] * pix(local_dst, x, y, xsize)->b;
      n = w[0];
      for ( wi=1; wi <= radius; wi++) {
	wc = w[wi];
	y2 = y - wi;
	if(y2 >= min) {
	  r += wc * pix(local_dst, x, y2, xsize)->r;
	  g += wc * pix(local_dst, x, y2, xsize)->g;
	  b += wc * pix(local_dst, x, y2, xsize)->b;
	  n += wc;
	}
	y2 = y + wi;
	if(y2 < max) {
	  r += wc * pix(local_dst, x, y2, xsize)->r;
	  g += wc * pix(local_dst, x, y2, xsize)->g;
	  b += wc * pix(local_dst, x, y2, xsize)->b;
	  n += wc;
	}
      }
      pix(local_src,x,y-radius, xsize)->r = r/n;
      pix(local_src,x,y-radius, xsize)->g = g/n;
      pix(local_src,x,y-radius, xsize)->b = b/n;
    }
  }
  
  const int gather_size = (ysize%num_proc)*xsize*pixel_type_size; 
  
  int recv_count[num_proc];
  int displs[num_proc];

  //Specify how much data will be recv from each process
  for(int i=0; i<num_proc; i++){
    recv_count[i] = (ysize/num_proc)*xsize*pixel_type_size;
    if(i == num_proc-1)
      recv_count[i] += (ysize%num_proc)*xsize*pixel_type_size;
    displs[i] = i * (ysize/num_proc)*xsize*pixel_type_size;
  }
  
  MPI_Gatherv(local_src, local_size*pixel_type_size, pixel_type, src,
	     recv_count, displs, pixel_type, root, com);
	     
}

int main(int argc, char** argv)
{
  MPI::Status status;
  MPI_Init( &argc, &argv );
  const MPI_Comm com = MPI_COMM_WORLD;
  const int root = 0; //root process
  const int myid = MPI::COMM_WORLD.Get_rank();  
  const int num_proc = MPI::COMM_WORLD.Get_size();

  int radius, xsize, ysize, colmax;
  pixel *src = new pixel[MAX_PIXELS];

  if(myid == 0){
    /* Take care of the arguments */

    if (argc != 4) {
	fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
	exit(1);
    }
    radius = atoi(argv[1]);
    if((radius > MAX_RAD) || (radius < 1)) {
	fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
	exit(1);
    }

    /* read file */
    if(read_ppm ("../test_img/im1.ppm", &xsize, &ysize, &colmax, (char *) src) != 0)
      exit(1);
    
    if (colmax > 255) {
      cerr << "Too large maximum color-component value" << endl;
      exit(1);
    }

    cout << "Has read the image, generating coefficients" << endl;
  }
  
  //Notify all processes about radius
  MPI_Bcast(&radius, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* filter */
  double* weights= new double[radius];
  get_gauss_weights(radius, weights);

  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Bcast(&xsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ysize, 1, MPI_INT, 0, MPI_COMM_WORLD);

  struct timespec stime, etime;
  if(myid == 0){
    cout << "Calling filter\n" << endl;
    clock_gettime(CLOCK_REALTIME, &stime);
  }

  blurfilter(xsize, ysize, src, radius, weights);

  clock_gettime(CLOCK_REALTIME, &etime);

  if(myid == 0){
    const int delta_time = ((etime.tv_sec - stime.tv_sec) 
			    + 1e-9*(etime.tv_nsec - stime.tv_nsec));
    cout << "Filtering took: " <<  delta_time << "secs" << endl;
    
    cout << "Writing output file" << endl;
    if( write_ppm(argv[3], xsize, ysize, (char *)src) != 0)
      exit(1);
  }

  //cleanup
  delete[] src;
  delete[] weights;

  MPI_Finalize();

  return 0;
}
