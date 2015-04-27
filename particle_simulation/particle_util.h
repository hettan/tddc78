#ifndef _particle_util_h
#define _particle_util_h

#include <iostream>

#include "definitions.h"

double f_rand(double fMin, double fMax);

//Generate random values for n particles
void init_particles(particle_t *particles, long n, int box_width, int box_height);

float check_collition(particle_t *p1, particle_t* collitions, const cord_t wall, const double time_steps);

void init_particles(particle_t *particles, long n, int box_width, int box_height);

void add_collitions(particle_t *p1, particle_t *p2, particle* collitions);

void unlink_particle(particle_t *p);

void particles_list_merge(particle_t *dst, particle_t *m);
#endif
