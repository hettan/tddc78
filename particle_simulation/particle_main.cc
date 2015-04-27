#include <iostream>
#include <mpi.h>
#include <time.h>

#include "coordinate.h"
#include "physics.h"
#include "definitions.h"
#include "particle_util.h"

#define BOX_WIDTH 10000
#define BOX_HEIGHT 10000
#define MAX_TIME 1000.0

using namespace std;

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
  const int num_particles = 10000*num_proc;
  
  particle_t particles[num_particles];
  if(myid == root)
    init_particles(particles, num_particles, BOX_WIDTH, BOX_HEIGHT); 
  
  //MPI_Bcast(&particles, num_particles*sizeof(particle), MPI_BYTE, root, com);
  
  float pressure = 0.0;
  cord_t wall{0, BOX_WIDTH, 0, BOX_HEIGHT};
  //Main loop: for each time-step do following
  double time_steps = STEP_SIZE; //start at time 1 
  while(time_steps < MAX_TIME){

    particle_t* collitions = nullptr;

    // - for all paricles do    
    particle_t *p = particles;
    p = particles;
    while(p){
      // - - check for collisions
      pressure += check_collition(p, collitions, wall, time_steps);
      p = p->next;
    }
    
    // - - Move paricles that has not collided with another.
    // - - Check for wall interaction and add the momentum
    time_steps += STEP_SIZE;
    
    //Re-add the collided particles
    particles_list_merge(particles, collitions);
  }
  
  // - Communicate if needed

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
