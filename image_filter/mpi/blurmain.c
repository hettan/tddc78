#include <mpi.h>

typedef struct _pixel 
{
  usigned char r,g,b;
} pixel;
  
int main(int argc, char** argv)
{
  MPI_Status status;
  char *string = "xxxxx"; //rcv buff
  int myid;

  MPI_Init( NULL, NULL );
  MPI_Comm_rank( MPI_COMM_WORLD, &myid);
  
  if(myid==2)
    MPI_Send( "Hello", 5, MPI_CHAR, 7, 1234, MPI_COMM_WORLD );
  
  if(myid==7) {
    MPI_Recv( string, 5, MPI_CHAR, 2, MPI_ANY_TAG,
	      MPU_COMM_WORLD, &status);
    printf( "Got %s from P&d, tag %d\n",
	    string, status.MPI_SOURCE, status.MPI_TAG);
  }
  MPI_Finalize();
  }
  return 0;
}
