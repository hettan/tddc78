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
  //Some high values that will never be reached
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

void send_to_neighbour(particle* data, const int count, const int dst,
		       MPI_Request* requests, const int myid){
  //The unique tag
  const int tag = stoi(to_string(dst) + to_string(myid));
  const int size = count * sizeof(particle_t);
  MPI_Isend(data, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, requests);
}

void prepare_send(particle_t*& p, particle_t** particles_send, int* particles_send_count, 
		  particle_t*& particles_first, particle_t*& collitions, int& particle_count,
		  const int local_xsize, const int local_ysize)
{
  if(p->pcord.x < 0){
    insert_particle_copy(p, particles_send[0], particles_send_count[0]++, particles_first, collitions);
    particle_count--;
  }
  else if(p->pcord.x > local_xsize){
    insert_particle_copy(p, particles_send[1], particles_send_count[1]++, particles_first, collitions);
    particle_count--;
  }
  else if(p->pcord.y < 0){
    insert_particle_copy(p,  particles_send[2], particles_send_count[2]++, particles_first, collitions);
    particle_count--;
  }
  else if(p->pcord.y > local_ysize){
    insert_particle_copy(p,  particles_send[3], particles_send_count[3]++, particles_first, collitions);
    particle_count--;
  }      
}

particle_t* recv_from_neighbour(int& total_count, const int from, const int myid){
  
  //The unique tag
  const int tag = stoi(to_string(myid) + to_string(from));

  MPI_Status status;
  MPI_Probe( from, tag, MPI_COMM_WORLD, &status);
  
  int size;
  MPI_Get_count(&status, MPI_BYTE, &size);
  const int count = size/sizeof(particle_t);
  total_count += count;

  particle_t* buff = static_cast<particle_t*>(malloc(size));  
  MPI_Recv(buff, size, MPI_BYTE, from, tag, MPI_COMM_WORLD, &status);
  return buff;
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
  
  //particle_t* particles_buff = new particle_t[num_particles];
  particle_t* particles_buff = init_particles(local_num_particles, local_xsize, local_ysize);
  particle_t* particles_first = particles_buff;
  
  //init_particles(particles_first, local_num_particles, local_xsize, local_ysize);
  int particle_count = local_num_particles;
  
  //get neighbouring grids
  int neighbours[4];
  get_neighbours(g_com, myid, neighbours);
  
  //Setup the wall
  float pressure = 0.0;
  float local_pressure = 0.0;
  cord_t local_wall;
  set_local_wall(neighbours, local_xsize, local_ysize, local_wall);
  
  //Communication variables
  particle_t* particles_send[4];
  init_particles_send(particles_send);
  int particles_send_count[4] = {0,0,0,0};
  MPI_Request requests[4];

  //Main loop: for each time-step do following
  double time_steps = STEP_SIZE; //start at time 1 
  while(time_steps < MAX_TIME){
    int col_counter=0;
    particle_t* collitions=nullptr;
    
    // - for all paricles do    
    int add_pressure;
    particle_t* p = particles_first;
    p = particles_first;
  
    while(p){
      //Save and second-next of p since it might get reassigned
      particle_t* next = p->next;     
      particle_t* second_next = nullptr;
      if(next)
	second_next = next->next;
      
      /* Check for collition, if true then the pixels are moved to
	 collitions list instead. If false, the pixel will be moved. */
      check_collition(p, collitions);
          
      // - - Move paricles that has not collided with another.
      if(collitions != particles_first)
	feuler(&(p->pcord), time_steps);
      
      // - - Check for wall interaction and add the momentum
      add_pressure = wall_collide(&(p->pcord), local_wall);      
      local_pressure += add_pressure;

      prepare_send(p, particles_send, particles_send_count, particles_first, 
		   collitions, particle_count, local_xsize, local_ysize);
      
      //If there was a collition check the next particle as well
      if(collitions == p){
	prepare_send(next, particles_send, particles_send_count, particles_first, 
		     collitions, particle_count, local_xsize, local_ysize);
	
	//Skip next if it collided with p
	if(collitions->next  == next)
	  next = second_next;
	
	col_counter += 2;
      }
     
      //Assign new first to particles_first
      if(collitions == particles_first) {
	particles_first = next;
      }
     
      p = next;
    }

    //Re-add the collided particles
    particles_list_merge(particles_first, collitions);
    
    // - Communicate if needed
    for(int i=0; i<4; i++){
      if(neighbours[i] >= 0){
	send_to_neighbour(particles_send[i], particles_send_count[i], 
			  neighbours[i], requests, myid);

	const int prev_particle_count = particle_count;
	particle_t* recv_buff = recv_from_neighbour(particle_count, neighbours[i], myid);
      
	const int recv_particle_count = (particle_count - prev_particle_count);
	insert_particles(recv_buff, recv_particle_count, particles_first);
      }
      particles_send_count[i] = 0;
    }

    //make sure all is done before begining next timestep
    MPI_Barrier(com);
    time_steps += STEP_SIZE;
  }  
  
  //Calculate pressure
  MPI_Reduce(&local_pressure, &pressure, 1, MPI_DOUBLE, MPI_SUM, root, g_com);

  if(myid == 0){
    clock_gettime(CLOCK_REALTIME, &etime);
    const int delta_time = ((etime.tv_sec - stime.tv_sec) 
			    + 1e-9*(etime.tv_nsec - stime.tv_nsec));
    cout << "Calculation took: " <<  delta_time << "secs" << endl;

    //Should temperature be updated in the algorithm instead?
    const double temperature = 1.0;
    print_results(pressure, BOX_WIDTH*BOX_HEIGHT, num_particles, temperature);
  }
  
  //cleanup
  cleanup_particles_send(particles_send);
  particles_cleanup(particles_first);

  MPI_Finalize();

  return 0;
}
