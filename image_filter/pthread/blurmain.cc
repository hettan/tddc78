#include <stdlib.h>
//#include <stdio.h>
#include <iostream>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <semaphore.h>

#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"

using namespace std;

#define MAX_THREADS 50

int get_y_start(const int thread_id, const int num_threads, const int ysize) {
  return thread_id * (ysize/num_threads);
}

//add the leftoverpixels to the last thread
int get_y_end(const int thread_id, const int num_threads, const int ysize) {
  int y_end = (thread_id+1) * (ysize/num_threads);
  if(thread_id == num_threads-1)
    y_end += ysize%num_threads;
  return y_end;
}

int main (int argc, char ** argv) {
  int radius;
  int xsize, ysize, colmax;
  pixel* src = new pixel[MAX_PIXELS];
  pixel* dst = new pixel[MAX_PIXELS];
  struct timespec stime, etime;
#define MAX_RAD 1000

  double w[MAX_RAD];

  /* Take care of the arguments */

  if (argc != 5) {
    cerr << "Usage: " << argv[0] <<" radius threads infile outfile" << endl;
    return 1;
  }
  radius = atoi(argv[1]);
  if((radius > MAX_RAD) || (radius < 1)) {
    cerr << "Radius (" << radius << ") must be greater than zero and less then "
	 << MAX_RAD << endl;
    return 2;
  }

  /* read file */
  if(read_ppm (argv[3], &xsize, &ysize, &colmax, (char *) src) != 0)
    return 3;

  if (colmax > 255) {
    cerr << "Too large maximum color-component value" << endl;
    return 4;
  }
  
  const int num_threads = atoi(argv[2]);
  if(num_threads < 1){
    cerr << "Number of threads need to be higher than zero" << endl;
    return 5;
  }

  cout << "Has read the image, generating coefficients" << endl;

  /* filter */
  get_gauss_weights(radius, w);

  cout << "Calling filter" << endl;

  clock_gettime(CLOCK_REALTIME, &stime);

  //Start the threads and blurfilter_x
  pthread_t threads[MAX_THREADS];
  struct thread_data_blurfilter thread_data[MAX_THREADS];
  for(int i=0; i<num_threads; i++){
    thread_data[i].threadId = i;
    thread_data[i].xsize = xsize;
    thread_data[i].ysize = ysize;

    thread_data[i].src = src;
    thread_data[i].dst = dst;
    thread_data[i].radius = radius;
    thread_data[i].w = w;
    
    thread_data[i].y_start = get_y_start(i, num_threads, ysize);
    thread_data[i].y_end = get_y_end(i, num_threads, ysize);
   
    int rc = pthread_create( &threads[i], NULL, blurfilter_x, (void*)&thread_data[i] );
    if(rc){
      cout << "Error creating thread, " << rc << endl;
      return 6;
    } 
  }

  //Wait for all threads to finish
  void* status;
  for(int i=0; i<num_threads; i++){
    int rc = pthread_join(threads[i], &status);
    if(rc){
      cout << "Error in thread join: " << rc << endl;
      return 7;
    }
    
    //cout << "Thread " << i << " done, status=" << status << endl;
  }

  
  //Start the threads again but this time with blurfilter_y
  for(int i=0; i<num_threads; i++){
    thread_data[i].threadId = i;
    thread_data[i].xsize = xsize;
    thread_data[i].ysize = ysize;

    thread_data[i].src = src;
    thread_data[i].dst = dst;
    thread_data[i].radius = radius;
    thread_data[i].w = w;

    thread_data[i].y_start = get_y_start(i, num_threads, ysize);
    thread_data[i].y_end = get_y_end(i, num_threads, ysize);
   
    int rc = pthread_create( &threads[i], NULL, blurfilter_y, (void*)&thread_data[i] );
    if(rc){
      cout << "Error creating thread, " << rc << endl;
      return 6;
    } 
  }

  //Wait for all threads to finish
  for(int i=0; i<num_threads; i++){
    int rc = pthread_join(threads[i], &status);
    if(rc){
      cout << "Error in thread join: " << rc << endl;
      return 7;
    }
    
    //cout << "Thread " << i << " done, status=" << status << endl;
  }
  
  clock_gettime(CLOCK_REALTIME, &etime);
  
  const int delta_time = (etime.tv_sec  - stime.tv_sec) + 1e-9*(etime.tv_nsec  - stime.tv_nsec);
  cout << "Filtering took: "<< delta_time <<" secs" << endl;

  /* write result */
  cout << "Writing output file" << endl;
    
  if(write_ppm (argv[4], xsize, ysize, (char *)src) != 0)
    return 8;

  //Cleanup
  delete[] src;
  delete[] dst;

  pthread_exit(NULL);

  return 0;
}
