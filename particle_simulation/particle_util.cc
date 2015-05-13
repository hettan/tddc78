#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <mpi.h>
#include <algorithm>

#include "particle_util.h"
#include "physics.h"

using namespace std;

//Triolith doesn't have c++11 lib in mpi version
#define nullptr NULL

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
  ostream& operator<<(ostream& os, cord_t *c){
    os << "x0=" << c->x0 << ", x1=" << c->x1 
       << ", y0=" << c->y0 << ", y1=" << c->y1;
    return os;
  }
  ostream& operator<<(ostream& os, const cord_t& c){
    os << "x0=" << c.x0 << ", x1=" << c.x1 
       << ", y0=" << c.y0 << ", y1=" << c.y1;
    return os;
  }
}

double f_rand(double fMin, double fMax)
{
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}


//Generate random values for n particles
particle_t* init_particles(long n, int box_width, int box_height)
{
  srand(time(NULL));
  particle_t* first = nullptr;  
  particle_t* prev = nullptr;
  
  for(int i=0; i<n; i++){
    particle_t* p = static_cast<particle_t*>(malloc(sizeof(particle_t)));
    if(i == 0)
      first = p;

    p->pcord.x = f_rand(0, box_width);
    p->pcord.y = f_rand(0, box_height);
    
    const double r = f_rand(0, MAX_INITIAL_VELOCITY);
    const double v = f_rand(0, 2*PI);
    p->pcord.vx = r*cos(v);
    p->pcord.vy = r*sin(v);

    p->prev = prev;
    
    if(prev)
      prev->next = p;
    
    if(i == n-1)
      p->next = nullptr;
    
    prev = p;
  }
  return first;
}

void particles_cleanup(particle_t*& current)
{
  particle_t* next;

  int counter = 0;
  while(current){
    next = current->next;
    free(current);
    current = next;
  }
}

void insert_particle(particle_t* p , particle_t*& arr)
{
  if(arr){
    arr->prev = p;
    p->next = arr;
  }
  else
    p->next = nullptr;
    
  p->prev = nullptr;
  arr = p;
}

void insert_particles(particle_t* add_particles, const int count, particle_t*& arr)
{
  for(int i=0; i<count; i++) {
    insert_particle(add_particles+i, arr);
  }
}

void insert_particle_copy(particle_t*& p, particle_t* arr, int index, 
			  particle_t*& particles_first, particle_t*& collitions)
{ 
  
  unlink_particle(p);
  arr[index] = *p;

  remove_particle(p, particles_first);
  remove_particle(p, collitions);
  free(p);
}

void remove_particle(particle_t* p, particle_t*& arr)
{
  if(p == arr)
    arr = arr->next;
}

//Remove and fix prev/next pointers
void unlink_particle(particle_t *p)
{
  if(p->prev)
    (p->prev)->next = p->next;
  
  if(p->next)
    (p->next)->prev = p->prev;
}

//Add m to the front of the linking chain of dst, return new start particle
void particles_list_merge(particle_t*& dst, particle_t*& m)
{
  if (m == nullptr)
    return;

  if (dst == nullptr) {
    dst = m;
    return;
  }

  //move to last elem in m
  particle* first_m = m;
  while(m->next)
    m = m->next;
  
  m->next = dst;
  dst->prev = m;  
  
  dst = first_m;
  
  //set m to "empty"
  m = nullptr;
}

int get_size(particle_t* arr)
{
  int counter=0;
  while(arr){    
    counter++;
    arr=arr->next;
  }
  
  return counter;
}

void init_particles_send(particle_t** send)
{
  for(int i=0; i<4; i++)
    send[i] = new particle[SEND_BUFF_SIZE];
}

void cleanup_particles_send(particle_t** send)
{
  for(int i=0; i<4; i++)
    delete[] send[i];
}

//Check for collition of other particels and walls, returning the wall_collide pressure
float check_collition(particle_t* p1, particle_t*& collitions)
{
  if (p1 == nullptr)
    return 0.0;

  particle* p2 = p1->next;
  while(p2){
    const float c = collide(&(p1->pcord), &((p2)->pcord));
    if( c != -1 ){
      add_collitions(p1, p2, collitions);
      interact(&(p1->pcord), &(p2->pcord), c);
      break;
    }
    p2 = p2->next;
  }
  
  return 0.0;
} 

//insert particles in front of collitions, making p1 the new front followed by p2
void add_collitions(particle_t* p1, particle_t* p2, particle_t*& collitions)
{
  unlink_particle(p1);
  unlink_particle(p2);

  //prepare for collitions linking
  p1->prev = nullptr;
  p1->next = p2;
  p2->prev = p1;

  //If there's particles added already link together with p2
  if(collitions){
    p2->next = collitions;
    collitions->prev = p2;    
  }
  else
    p2->next = nullptr;
 
  collitions = p1;
}
