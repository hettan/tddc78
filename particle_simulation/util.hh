#ifndef _UTIL_HH
#define _UTIL_HH

#include <mpi.h>

#include "definitions.h"
#include "particle_list.hh"

#define nullptr NULL

namespace std{
  ostream& operator<<(ostream& os, particle_t *p);
  ostream& operator<<(ostream& os, const particle_t& p);
  ostream& operator<<(ostream& os, ParticleList* l);
}

particle_t* get_collition(particle_t* p1);

bool prepare_send(particle_t* p, ParticleList* particles_send,
		   ParticleList* particles,
		   const int local_xsize, const int local_ysize,
		  const int* neighbours);

void print_send_list(ParticleList* send_list);

void send_to_neighbour(ParticleList& send_list, 
		       const int dst, const int myid,
			MPI_Request* requests, 
			MPI_Comm& com);

void recv_from_neighbour(ParticleList* particles,
			  const int from, const int myid,
			  MPI_Comm& com);

void set_local_wall(int* neighbours, const int xsize,
		    const int ysize, cord_t& wall);

void get_neighbours(const MPI_Comm& g_com, int myid, int* neighbours);


void print_results(const float pressure, const long area,
		   const int num_particles, const int temperature);

#endif
