#include <mpi.h>
#include <iostream>

struct pixel
{
  unsigned char r,g,b;
};

using namespace std;

namespace std{
  ostream& operator<<(ostream& os, pixel& p){
    os << "r=" << p.r << ", g=" << p.g << ", b=" << p.b;
    return os;
  }
} 

int main(int argc, char** argv)
{
  MPI::Status status;
  int myid;
  MPI::Init( argc, argv );
  myid = MPI::COMM_WORLD.Get_rank();
  cout << "My id is " << myid << endl;

  int size = 3;
  int blockcounts[] = {1, 1, 1};
  MPI_Aint disp[] = {0, 0, 0};
  MPI_Datatype old_types[3] = {MPI::UNSIGNED_CHAR, MPI::UNSIGNED_CHAR, MPI::UNSIGNED_CHAR};

  MPI_Datatype pixel_type;
  MPI_Type_struct(size, blockcounts, disp, old_types, &pixel_type);
  MPI_Type_commit(&pixel_type);

  char buf[5];

  pixel test;

  int tag = 1;
  int rc;
  if(myid==0){
    int dst = 1;
    int src = 1;
    
    test.r = 'a';
    test.g = 'b';
    test.b = 'c';
    MPI::COMM_WORLD.Send( &test, 3, pixel_type, dst, tag );
    //string msg = "Hello";
    //MPI::COMM_WORLD.Send( msg.c_str(), 5, MPI::CHAR, dst, tag );
    //MPI::COMM_WORLD.Recv( buf, 5, MPI::CHAR, src, tag, status);
    //cout << "Got " << buf << " from " << src << endl;
  }
  if(myid==1) {
    int src = 0;
    int dst = 0;
    
    MPI::COMM_WORLD.Recv( &test, 3, pixel_type, src, tag);
    cout << test << endl;
    /*
    MPI::COMM_WORLD.Recv( buf, 5, MPI::CHAR, src, tag, status);
    cout << "Got " << buf << " from " << src << endl;
    MPI::COMM_WORLD.Send( buf, 5, MPI::CHAR, dst, tag ); */
  }
  MPI::Finalize();
  return 0;
}
