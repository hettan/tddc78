#ifndef _PARTICLE_LIST_HH
#define _PARTICLE_LIST_HH

#include <iostream>
#include <mpi.h>

#include "definitions.h"

using namespace std;

//#define nullptr NULL
#define MAX_SEND_BUFF_SIZE 100000

class ParticleList
{
public:
  ParticleList();
  ~ParticleList();

  void create_random_particles(const long n, const int box_width, 
			       const int box_height);
  void insert(particle_t* p);
  particle_t* remove(particle_t* p);
  void erase(particle_t* p);
  particle_t* pop();
  void insert_all(ParticleList* l);
  void clear(){ first=nullptr; size=0; }

  particle_t* get_first(){ return first; }
  int get_size() const{ return size; }
  
  particle_t* get_arr_copy();
      
private:
  double f_rand(const double f_min, const double f_max);

  particle_t* first;
  particle_t* arr;
  int size;
};

int get_size(particle_t*);



#endif // _PARTICLE_LIST_HH
