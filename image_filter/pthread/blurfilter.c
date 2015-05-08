/*
  File: blurfilter.c

  Implementation of blurfilter function.
    
 */
#include <stdio.h>
#include <pthread.h>
#include <semaphore.h>

#include "blurfilter.h"
#include "ppmio.h"

pixel* pix(pixel* image, const int xx, const int yy, const int xsize)
{
  register int off = xsize*yy + xx;

#ifdef DBG
  if(off >= MAX_PIXELS) {
    fprintf(stderr, "\n Terribly wrong: %d %d %d\n",xx,yy,xsize);
  }
#endif
  return (image + off);
}

void* blurfilter_x(void* t_param){
  struct thread_data_blurfilter* thread_data = (thread_data_blurfilter*)t_param;
  const int xsize = thread_data->xsize;

  pixel* src = thread_data->src;
  pixel* dst = thread_data->dst;
  const int radius = thread_data->radius;
  const double* w = thread_data->w;

  const int y_start = thread_data->y_start;
  const int y_end = thread_data->y_end;

  int x,y,x2, wi;
  double r,g,b,n, wc;
  
  for (y=y_start; y<y_end; y++) {
    for (x=0; x<xsize; x++) {
      r = w[0] * pix(src, x, y, xsize)->r;
      g = w[0] * pix(src, x, y, xsize)->g;
      b = w[0] * pix(src, x, y, xsize)->b;
      n = w[0];
      for ( wi=1; wi <= radius; wi++) {
	wc = w[wi];
	x2 = x - wi;
	if(x2 >= 0) {
	  r += wc * pix(src, x2, y, xsize)->r;
	  g += wc * pix(src, x2, y, xsize)->g;
	  b += wc * pix(src, x2, y, xsize)->b;
	  n += wc;
	}
	x2 = x + wi;
	if(x2 < xsize) {
	  r += wc * pix(src, x2, y, xsize)->r;
	  g += wc * pix(src, x2, y, xsize)->g;
	  b += wc * pix(src, x2, y, xsize)->b;
	  n += wc;
	}
      }
      pix(dst,x,y, xsize)->r = r/n;
      pix(dst,x,y, xsize)->g = g/n;
      pix(dst,x,y, xsize)->b = b/n;
    }
  }
  
  return nullptr;
}

void* blurfilter_y(void* t_param){
  struct thread_data_blurfilter* thread_data = (thread_data_blurfilter*)t_param;
  const int xsize = thread_data->xsize;
  const int ysize = thread_data->ysize;

  pixel* src = thread_data->src;
  pixel* dst = thread_data->dst;  
  const int radius = thread_data->radius;
  const double* w = thread_data->w;
  
  const int y_start = thread_data->y_start;
  const int y_end = thread_data->y_end;

  int x,y,y2, wi;
  double r,g,b,n, wc;
  
  for (y=y_start; y<y_end; y++) {
    for (x=0; x<xsize; x++) {
      r = w[0] * pix(dst, x, y, xsize)->r;
      g = w[0] * pix(dst, x, y, xsize)->g;
      b = w[0] * pix(dst, x, y, xsize)->b;
      n = w[0];
      for ( wi=1; wi <= radius; wi++) {
	wc = w[wi];
	y2 = y - wi;
	if(y2 >= 0) {
	  r += wc * pix(dst, x, y2, xsize)->r;
	  g += wc * pix(dst, x, y2, xsize)->g;
	  b += wc * pix(dst, x, y2, xsize)->b;
	  n += wc;
	}
	y2 = y + wi;
	if(y2 < ysize) {
	  r += wc * pix(dst, x, y2, xsize)->r;
	  g += wc * pix(dst, x, y2, xsize)->g;
	  b += wc * pix(dst, x, y2, xsize)->b;
	  n += wc;
	}
      }
      pix(src,x,y, xsize)->r = r/n;
      pix(src,x,y, xsize)->g = g/n;
      pix(src,x,y, xsize)->b = b/n;
    }
  }
  
  for(y=y_start; y<y_end; y++)
    printf("%d", pix(src, xsize-1, y, xsize)->r);
  
  return nullptr;
  
}



