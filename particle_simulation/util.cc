#include <string>
#include <sstream>

#include "util.hh"
#include "physics.h"

using namespace std;

namespace std{
  ostream& operator<<(ostream& os, particle_t *p){
    os << "x=" << p->pcord.x << ", y=" << p->pcord.y 
       << ", vx=" << p->pcord.vx << ", vy=" << p->pcord.vy;
    return os;
  }
  ostream& operator<<(ostream& os, const particle_t& p){
    os << "x=" << p.pcord.x << ", y=" << p.pcord.y 
       << ", vx=" << p.pcord.vx << ", vy=" << p.pcord.vy;
    return os;
  }
  ostream& operator<<(ostream& os, ParticleList* l){
    particle_t* p = l->get_first();
    while(p){
      os << p << endl;
      p = p->next;
    }
    return os;
  }
}

string to_string(const int i){
  ostringstream oss;
  oss << i;
  return oss.str();
}

int to_int(const string str){
  int i;
  istringstream (str) >> i;
  return i;
}

//Check for collition of other particels and walls, returning the wall_collide pressure
particle_t* get_collition(particle_t* p1)
{
  if (p1 == nullptr)
    return nullptr;

  particle* p2 = p1->next;
  while(p2){
    const float c = collide(&(p1->pcord), &((p2)->pcord));
    if( c != -1 ){

      interact(&(p1->pcord), &(p2->pcord), c);
      return p2;
    }
    p2 = p2->next;
  }
  
  return nullptr;
} 

bool prepare_send(particle_t* p, ParticleList* particles_send, 
		   ParticleList* particles,
		   const int local_xsize, const int local_ysize,
		   const int* neighbours)
{
  //const int myid = MPI::COMM_WORLD.Get_rank();  
  if(p->pcord.x < 0 && neighbours[0] >= 0){
    particles->remove(p);
    particles_send[0].insert(p);
  }
  else if(p->pcord.x > local_xsize && neighbours[1] >= 0){
    particles->remove(p);
    particles_send[1].insert(p);
  }
  else if(p->pcord.y < 0 && neighbours[2] >= 0){
    particles->remove(p);
    particles_send[2].insert(p);
  }
  else if(p->pcord.y > local_ysize && neighbours[3] >= 0){
    particles->remove(p);
    particles_send[3].insert(p);
  }      
  else{
    return false;
  }

  return true;
}

void print_send_list(ParticleList* send_list){
  cout << "#send_list count#"<< endl;
  for(int i=0;i<4;i++)
    cout << i << "=>" << send_list[i].get_size() << endl;
}

void send_to_neighbour(ParticleList& send_list,
			const int dst, const int myid,
			MPI_Request* requests,
			MPI_Comm& com){
  particle_t* buff = send_list.get_arr_copy();
  const int count = send_list.get_size();

  //The unique tag
  const int tag = to_int(to_string(dst) + to_string(myid));
  const int size = count * sizeof(particle_t);
  //cout << myid << ": sending " << count << " to " << dst << endl;
  MPI_Isend(buff, size, MPI_BYTE, dst, tag, com, requests);
}

void recv_from_neighbour(ParticleList* particles,
			  const int from, const int myid,
			  MPI_Comm& com){
  
  //The unique tag
  const int tag = to_int(to_string(myid) + to_string(from));

  MPI_Status status;
  MPI_Probe( from, tag, com, &status);
  
  int size;
  MPI_Get_count(&status, MPI_BYTE, &size);
  const int count = size/sizeof(particle_t);
  
  particle_t* buff = static_cast<particle_t*>(malloc(size));  
  MPI_Recv(buff, size, MPI_BYTE, from, tag, com, &status);
  //cout << myid << "recved " << count << " from " << from << endl;
  
  //Add the particles list
  particle_t* p;
  for(int i=0; i<count; i++){
    p = static_cast<particle_t*>(malloc(sizeof(particle_t)));
    p[0] = *(buff+i);
    particles->insert(p);
  }
  free(buff);
}

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
