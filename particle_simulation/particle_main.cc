#include <iostream>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <sstream>

#include "coordinate.h"
#include "physics.h"
#include "definitions.h"
#include "util.hh"
#include "particle_list.hh"

#define BOX_WIDTH 10000
#define BOX_HEIGHT 10000
#define MAX_TIME 10.0

//Triolith doesn't have c++11 lib in mpi version
#define nullptr NULL

using namespace std;

int main(int argc, char** argv)
{
  //MPI init
  MPI::Status status;
  MPI_Init( &argc, &argv );
  const int root = 0; //root process
  const int myid = MPI::COMM_WORLD.Get_rank();  
  const int num_proc = MPI::COMM_WORLD.Get_size();

  struct timespec stime, etime;
  if(myid == root){
    clock_gettime(CLOCK_REALTIME, &stime);
  }
  
  int dims[2] = {0, 0};
  int periods[2] = {0, 0};

  //Specify the grid layout
  MPI_Dims_create(num_proc, 2, dims);

  MPI_Comm g_com;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &g_com);
  
  int grid_coord[2];
  MPI_Cart_coords(g_com, myid, 2, grid_coord);
  
  double local_xsize = BOX_WIDTH/(dims[0]);
  double local_ysize = BOX_HEIGHT/(dims[1]);

  //Generate particles on each grid 
  const int num_particles = 10000*num_proc;
  int local_num_particles = num_particles/num_proc;

  //Add leftover particles to the last process
  if(myid == num_proc - 1){
    local_num_particles += num_particles%num_proc;
  }
  
  ParticleList* particles = new ParticleList();
  particles->create_random_particles(local_num_particles, local_xsize, local_ysize);
  ParticleList* collitions = new ParticleList();

  //get neighbouring grids
  int neighbours[4];
  get_neighbours(g_com, myid, neighbours);
  
  ParticleList particles_send[4];
  MPI_Request requests[4];

  //Setup the wall
  float pressure = 0.0;
  float local_pressure = 0.0;
  cord_t local_wall;
  set_local_wall(neighbours, local_xsize, local_ysize, local_wall);

  //Main loop: for each time-step do following
  double time_steps = STEP_SIZE; //start at time 1 
  while(time_steps < MAX_TIME){
    
    // - for all paricles do    
    int add_pressure;
    particle_t* p = particles->get_first();

    int counter = 0;
    while(p){
      counter++;

      //Save and second-next of p since it might get reassigned
      particle_t* next = p->next;     
      particle_t* second_next = nullptr;
      if(next)
	second_next = next->next;
      
      /* Check for collition, if true then the pixels are moved to
	 collitions list instead. If false, the pixel will be moved. */
      particle* collided_particle;
      if( (collided_particle = get_collition(p)) ){
	
	//Skip next if it collided with p
	if(collided_particle == next){
	  next = second_next;
	}
	
	if(! prepare_send(collided_particle, particles_send, particles,
			  local_xsize, local_ysize, neighbours)){
	  particles->remove(collided_particle);
	  collitions->insert(collided_particle);
	}
      }   
      else {
	feuler(&(p->pcord), time_steps);
      }
      
      // - - Check for wall interaction and add the momentum
      add_pressure = wall_collide(&(p->pcord), local_wall);      
      local_pressure += add_pressure;
      if(! prepare_send(p, particles_send, particles,
			local_xsize, local_ysize, neighbours)){
	
	if(collided_particle){
	  particles->remove(p);
	  collitions->insert(p);
	}
      }
      p = next;
    }
    //Re-add the collided particles
    particles->insert_all(collitions);

    // - Communicate if needed
    for(int i=0; i<4; i++){
      if(neighbours[i] >= 0){
	send_to_neighbour(particles_send[i], neighbours[i], myid, requests, g_com); 
	recv_from_neighbour(particles, neighbours[i], myid, g_com);
      }
      particles_send[i].clear();
    }

    //make sure all is done before begining next timestep
    MPI_Barrier(g_com);
    time_steps += STEP_SIZE;
  }  

  //Calculate pressure
  MPI_Reduce(&local_pressure, &pressure, 1, MPI_DOUBLE, MPI_SUM, root, g_com);
  if(myid == root){
    clock_gettime(CLOCK_REALTIME, &etime);
    const int delta_time = ((etime.tv_sec - stime.tv_sec) 
			    + 1e-9*(etime.tv_nsec - stime.tv_nsec));
    cout << "Calculation took: " <<  delta_time << "secs" << endl;

    //Should temperature be updated in the algorithm instead?
    const double temperature = 1.0;
    print_results(pressure, BOX_WIDTH*BOX_HEIGHT, num_particles, temperature);
  }
  
  //cleanup
  delete particles;
  delete collitions;
  MPI_Finalize();

  return 0;
}
