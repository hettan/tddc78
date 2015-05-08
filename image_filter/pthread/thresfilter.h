/*
  File: thresfilter.h

  Declaration of pixel structure and thresfilter function.
    
 */
#ifndef _THRESFILTER_H_
#define _THRESFILTER_H_
/* NOTE: This structure must not be padded! */
typedef struct _pixel {
    unsigned char r,g,b;
} pixel;

struct thread_data_thresfilter{
  int threadId;
  int xsize;
  int ysize;
  pixel* src;
};


void* thresfilter(void* t_param);

#endif
