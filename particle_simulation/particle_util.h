#ifndef _particle_util_h
#define _particle_util_h

#include <iostream>

#include "definitions.h"

#define SEND_BUFF_SIZE 1000

namespace std{
  ostream& operator<<(ostream& os, particle_t* p);
  ostream& operator<<(ostream& os, const particle_t& p);

  ostream& operator<<(ostream& os, cord_t* c);
  ostream& operator<<(ostream& os, const cord_t& c);
}

double f_rand(double fMin, double fMax);

//Generate random values for n particles
particle_t* init_particles(long n, int box_width, int box_height);
void particles_cleanup(particle_t*& first);

void insert_particle(particle_t* p , particle_t*& arr);
void insert_particles(particle_t* particles, const int count, particle_t*& arr);
void insert_particle_copy(particle_t*& p, particle_t* arr, int index, particle_t*& particles_first, particle_t*& collitions);
void remove_particle(particle_t* p, particle_t*& arr);
void unlink_particle(particle_t* p);
void particles_list_merge(particle_t*& dst, particle_t*& m);

int get_size(particle_t* arr);

void init_particles_send(particle** send);
void cleanup_particles_send(particle** send);

float check_collition(particle_t* p1, particle_t*& collitions);
void add_collitions(particle_t* p1, particle_t* p2, particle*& collitions);




#endif
