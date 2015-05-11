#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <semaphore.h>

#include "ppmio.h"
#include "thresfilter.h"

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
  int xsize, ysize, colmax;
  pixel* src = new pixel[MAX_PIXELS];
  struct timespec stime, etime;
    
  /* Take care of the arguments */

  if (argc != 4) {
    cerr << "Usage: " << argv[0] <<" threads infile outfile" << endl;
    return 1;
  }

  /* read file */
  if(read_ppm (argv[2], &xsize, &ysize, &colmax, (char *) src) != 0)
    return 2;

  if (colmax > 255) {
    cerr << "Too large maximum color-component value" << endl;
    return 3;
  }
     
  const int num_threads = atoi(argv[1]);
  if(num_threads < 1){
    cerr << "Number of threads need to be higher than zero" << endl;
    return 4;
  }


  cout << "Has read the image, calling filter" << endl;

  clock_gettime(CLOCK_REALTIME, &stime);
  
  pthread_t threads[MAX_THREADS];
  struct thread_data_thresfilter thread_data[MAX_THREADS];
  
  //Start the threads
  for(int i=0; i<num_threads; i++){
    thread_data[i].threadId = i;
    thread_data[i].xsize = xsize;

    //Set the y_interval to be processed by the thread
    const int y_start = get_y_start(i, num_threads, ysize);
    thread_data[i].src = src + y_start;
    thread_data[i].ysize = get_y_end(i, num_threads, ysize);
   
    int rc = pthread_create( &threads[i], NULL, thresfilter, (void*)&thread_data[i] );
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
    
    //    cout << "Thread " << i << " done, status=" << status << endl;
  }

  clock_gettime(CLOCK_REALTIME, &etime);

  const int delta_time = (etime.tv_sec  - stime.tv_sec) + 1e-9*(etime.tv_nsec  - stime.tv_nsec);
  cout << "Filtering took: "<< delta_time <<" secs" << endl;
  
  /* write result */
  cout << "Writing output file" << endl;
    
  if(write_ppm (argv[3], xsize, ysize, (char *)src) != 0)
    return 5;

  //cleanup
  delete[] src;

  return(0);
}
