#include "thresfilter.h"

void* thresfilter(void* t_param){
#define uint unsigned int 
  
  struct thread_data_thresfilter* thread_data = (thread_data_thresfilter*)t_param;
  const int xsize = thread_data->xsize;
  const int ysize = thread_data->ysize;
  pixel* src = thread_data->src;

  uint sum, i, psum, nump;

  nump = xsize * ysize;

  for(i = 0, sum = 0; i < nump; i++) {
    sum += (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
  }

  sum /= nump;

  for(i = 0; i < nump; i++) {
    psum = (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
    if(sum > psum) {
      src[i].r = src[i].g = src[i].b = 0;
    }
    else {
      src[i].r = src[i].g = src[i].b = 255;
    }
  }

  return nullptr;
}
