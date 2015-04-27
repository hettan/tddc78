#include <math.h>
#include <time.h>
#include <stdlib.h>

#include "particle_util.h"
#include "physics.h"

using namespace std;

namespace std{
  ostream& operator<<(ostream& os, particle_t *p){
    os << "x=" << p->pcord.x << ", y=" << p->pcord.y 
       << ", vx=" << p->pcord.vx << ", vy=" << p->pcord.vy;
    return os;
  }
}

double f_rand(double fMin, double fMax)
{
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

//Generate random values for n particles
void init_particles(particle_t *particles, long n, int box_width, int box_height)
{
  srand(time(NULL));
  for(int i=0; i<n; i++){
    particle_t *p = (particles+i);
    p->pcord.x = f_rand(0, box_width);
    p->pcord.y = f_rand(0, box_height);
    
    const double r = f_rand(0, MAX_INITIAL_VELOCITY);
    const double v = f_rand(0, 2*PI);
    p->pcord.vx = r*cos(v);
    p->pcord.vy = r*sin(v);

    //First elem doesn't have a prev
    if(i==0)
      p->prev = nullptr;
    else
      p->prev = (p-1);
    
    //Last elem doesn't have a next
    if(i == n-1)
      p->next = nullptr;
    else
    p->next = (p+1);
  }
}


//Check for collition of other particels and walls, returning the wall_collide pressure
float check_collition(particle_t *p1, particle_t *collitions, const cord_t wall, double time_steps)
{
  if (p1 == nullptr)
    return 0.0;

  particle *p2 = p1->next;
  while(p2){
    const float c = collide(&(p1->pcord), &((p2)->pcord));
    if( c != -1 ){
      cout << "collition! " << endl 
	   << "1: " << p1 << endl
	   << "2: " << p2 << endl;
      add_collitions(p1, p2, collitions);
      interact(&(p1->pcord), &(p2->pcord), c);
      break;
    }
    p2 = p2->next;
  }
  
  //If no collition move particle, p2 is used to check it last was reached or of break before.
  if(!p2)
    feuler(&(p1->pcord), time_steps);

  //Return the pressure even though a wall havn't' been it (then it will be 0)
  return wall_collide(&(p1->pcord), wall);
} 

//insert particles in front of collitions, making p1 the new front followed by p2
void add_collitions(particle_t *p1, particle_t *p2, particle* collitions)
{
  unlink_particle(p1);
  unlink_particle(p2);

  //prepare for collitions linking
  p1->prev = nullptr;
  p1->next = p2;
  p2->prev = p1;

  //If there's particles added already link together with p2
  if(collitions != nullptr){
    p2->next = collitions;
    collitions->prev = p2;    
  }
 
  collitions = p1;
}

//Remove and fix prev/next pointers
void unlink_particle(particle_t *p)
{
  if(p->prev != nullptr)
    (p->prev)->next = p->next;
  
  if(p->next != nullptr)
    (p->next)->prev = p->prev;
}

//Add m to the front of the linking chain of dst
void particles_list_merge(particle_t *dst, particle_t *m)
{
  if (m == nullptr)
    return;

  if (dst == nullptr) {
    dst = m;
    return;
  }

  //move to last elem in m
  while(m->next)
    m = m->next;
  
  m->next = dst;
  dst->prev = m;  
  
  //set m to "empty"
  m = nullptr;
}
