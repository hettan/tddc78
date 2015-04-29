#include <iostream>
#include <mpi.h>
#include <time.h>
#include <math.h>

#include "coordinate.h"
#include "physics.h"
#include "definitions.h"
#include "particle_util.h"

#define BOX_WIDTH 10000
#define BOX_HEIGHT 10000
#define MAX_TIME 1000.0

using namespace std;

//Returns the wall in the grid, if no wall the value will be negative and therefor not possible to reach.
void set_local_wall(int* neighbours, const int xsize,
		     const int ysize, cord_t& wall)
{
  const int no_wall_low = -10000;
  const int no_wall_high = 1000000;

  //Left wall
  if(neighbours[0] < 0)
    wall.x0 = 0;
  else
    wall.x0 = no_wall_low;
  

  //Right wall
  if(neighbours[1] < 0)
    wall.x1 = xsize;
  else
    wall.x1 = no_wall_high;

    
  //Upper wall
  if(neighbours[2] < 0)
    wall.y0 = 0;
  else
    wall.y0 = no_wall_low;
  

  //Lower wall
  if(neighbours[3] < 0)
    wall.y1 = ysize;
  else
    wall.y1 = no_wall_high;
  

}

//Returns neighbours, if neighbour rank is negative it means that it's not directly connected
void get_neighbours(const MPI_Comm& g_com, int myid, int* neighbours)
{
  //left
  MPI_Cart_shift(g_com, 1, -1, &myid, &(neighbours[0]));
  //right
  MPI_Cart_shift(g_com, 1, 1, &myid, &(neighbours[1]));
    //up
  MPI_Cart_shift(g_com, 0, -1, &myid, &(neighbours[2]));
  //down
  MPI_Cart_shift(g_com, 0, 1, &myid, &(neighbours[3]));

}

void print_results(const float pressure, const long area,
		   const int num_particles, const int temperature)
{
  const double R = (pressure*area) / (num_particles*temperature);
  cout << endl << "###RESULTS###" << endl
       << "Pressure: " << pressure << endl
       << "Area: " << area << endl
       << "Number of particles: " << num_particles << endl
       << "Magic constant R:  " << R << endl
       << "Temperature: " << temperature << endl;
}

int main(int argc, char** argv)
{
  //MPI init
  MPI::Status status;
  MPI_Init( &argc, &argv );
  const MPI_Comm com = MPI_COMM_WORLD;
  const int root = 0; //root process
  const int myid = MPI::COMM_WORLD.Get_rank();  
  const int num_proc = MPI::COMM_WORLD.Get_size();

  struct timespec stime, etime;
  if(myid == 0){
    clock_gettime(CLOCK_REALTIME, &stime);
  }
  
  //Initiate particles change numparticles to 10000*numproc
  //MPI_Bcast(&particles, num_particles*sizeof(particle), MPI_BYTE, root, com);
 
  int dims[2] = {0, 0};
  int periods[2] = {0, 0};

  //Specify the grid layout
  MPI_Dims_create(num_proc, 2, dims);

  MPI_Comm g_com;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &g_com);
  
  int grid_coord[2];
  MPI_Cart_coords(g_com, myid, 2, grid_coord);
  
  double local_xsize = BOX_WIDTH/(dims[0]+1);
  double local_ysize = BOX_HEIGHT/(dims[1]+1);
  
  //Generate particles on each grid 
  const int num_particles = 10000*num_proc;
  int local_num_particles = num_particles/num_proc;

  //Add leftover particles to the last process
  if(myid == num_proc - 1){
    local_num_particles += num_particles%num_proc;
  }

  particle_t particles[local_num_particles];
  
  init_particles(particles, local_num_particles, local_xsize, local_ysize);
  
  
  
  /*
  cout << "local_coord[0]=" << local_coord[0]
       << ", local_coord[1]=" << local_coord[1] << endl
       << "local_xsize=" << local_xsize
       << "  local_ysize=" << local_ysize << endl;
  */
  int neighbours[4];
  get_neighbours(g_com, myid, neighbours);
  cout << myid << ": grid_coord= (" << grid_coord[0] << ", " << grid_coord[1] << ")  -  ";
  for(int i=0; i<4; i++)
    cout << neighbours[i] << ", ";
  cout << endl;

  float pressure = 0.0;
  //cord_t wall{0, BOX_WIDTH, 0, BOX_HEIGHT};
  cord_t local_wall;
  set_local_wall(dims, local_xsize, local_ysize, local_wall);
  
  particle* send_left = nullptr;
  particle* send_right = nullptr;
  particle* send_up = nullptr;
  particle* send_down = nullptr;

  //Main loop: for each time-step do following
  double time_steps = STEP_SIZE; //start at time 1 
  while(time_steps < MAX_TIME){

    particle_t* collitions = nullptr;
    
    // - for all paricles do    
    int add_pressure;
    particle_t* p = particles;
    p = particles;
    while(p){
      /* Check for collition, if true then the pixels are moved to
	 collitions list instead. If false, the pixel will be moved. */
      add_pressure = check_collition(p, collitions, local_wall, time_steps);

      // - - Check for wall interaction and add the momentum
      if(add_pressure != 0)
	pressure += add_pressure;
      
      // - - Move paricles that has not collided with another.
      else
	feuler(&(p->pcord), time_steps);
      
      //Save next of p since it might get reassigned
      particle_t* next = p->next;
      if(p->pcord.x < 0){
	insert_particle(p, send_left);
	cout << myid << ": "<<  "Outside left " << p->pcord.x << endl;
      }
      else if(p->pcord.x > local_xsize){
	insert_particle(p, send_right);
	cout << myid << ": "<<  "Outside right " << p->pcord.x << endl;
      }
      else if(p->pcord.y < 0){
	insert_particle(p, send_up);
	cout << myid << ": "<<  "Outside up" << endl;
      }
      else if(p->pcord.y > local_ysize){
	insert_particle(p, send_down);
	cout << myid << ": "<< "Outside down" << endl;
      }
     
      p = next;
    }
    
    cout << myid << ": done! time " << time_steps << endl;
    
    // - Communicate if needed
    if((neighbours[0] >= 0) && (send_left == nullptr))
      cout << "should not be here!!!" << endl;

    //Re-add the collided particles
    particles_list_merge(particles, collitions);

    
    //make sure all is done before begining next timestep
    MPI_Barrier(com);
    time_steps += STEP_SIZE;
  }

  //Calculate pressure
  clock_gettime(CLOCK_REALTIME, &etime);
  
  if(myid == 0){
    const int delta_time = ((etime.tv_sec - stime.tv_sec) 
			    + 1e-9*(etime.tv_nsec - stime.tv_nsec));
    cout << "Calculation took: " <<  delta_time << "secs" << endl;
  
    const double temperature = 1.0;
    print_results(pressure, BOX_WIDTH*BOX_HEIGHT, num_particles, temperature);

  }

  MPI_Finalize();

  return 0;
}
