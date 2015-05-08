/*
  File: blurfilter.h

  Declaration of pixel structure and blurfilter function.
    
 */

#ifndef _BLURFILTER_H_
#define _BLURFILTER_H_

/* NOTE: This structure must not be padded! */
typedef struct _pixel {
    unsigned char r,g,b;
} pixel;


struct thread_data_blurfilter{
  int threadId;
  int xsize;
  int ysize;
  pixel* src;
  pixel* dst;
  int radius;
  double* w;
  int y_start;
  int y_end;
};


void* blurfilter_x(void* t_param);
void* blurfilter_y(void* t_param);

#endif
