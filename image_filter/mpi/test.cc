//#include <mpi.h>
#include <iostream>
#include <math.h>
//#include <cstdlib>

#define MAX_X 1.33
#define Pi 3.14159

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


void get_gauss_weights(int i, int n, double* weights_out)
{
  double x;
  for(; i<n; i++ ){
    x=(double)i * MAX_X/n;
    //weights_out[i] = exp(-x*x * Pi);
  }
}
 
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
}
*/
int main(int argc, char** argv)
{
  //MPI::Status status;
  //MPI::Init( argc, argv );
  int myid = 0;//MPI::COMM_WORLD.Get_rank();
  int size = 1;//MPI::COMM_WORLD.Get_size();
  
  //MPI_Datatype pixel_type;
  //init_pixel_type(pixel_type);

  pixel test[3];

  int tag = 1;
  int rc;
  if(myid==0){
    cout << "NUMBER OF THREADS = " << size << endl;
    int dst = 1;
    int src = 1;
    
    test[0] = pixel{'a','b','c'};
    test[1] = pixel{'d','e','f'};
    test[2] = pixel{'g','h','i'};
    //MPI::COMM_WORLD.Send( &test[0], 3, pixel_type, dst, tag );
    
  }
  else if(myid==1) {
    int src = 0;
    int dst = 0;
    
    // MPI::COMM_WORLD.Recv( &test[0], 3, pixel_type, src, tag);
    //cout << test[0] << endl;
  }

  //MPI::Finalize();
  return 0;
}
