#include <stdlib.h>
#include <string.h>

#include "darray.h"

struct darray{
  int *A;
  int a; /* Space alloced for array */
  int n; /* Number or elements in array */
  int i; /* Current iterator position */
};

darrays
allocDA(void
	)
{
  darrays d = (darrays) malloc(sizeof(struct darray));

  d->i = 0;
  d->n = 0;
  d->a = 4;
  d->A = (int *)malloc(d->a*sizeof(int));

  return d;
}

darrays
dupDA(darrays d
	)
{
  darrays dup = (darrays) malloc(sizeof(struct darray));
  d->A = NULL;
  return cpyDA(dup, d);
}

darrays
cpyDA(darrays dst,
      darrays src
	   )
{
  if (!dst) return dupDA(src);
  
  dst->i = src->i;
  dst->n = src->n;
  dst->a = src->a;
  dst->A = (int*)realloc(dst->A, dst->a*sizeof(int));
  memcpy(dst->A, src->A, src->a*sizeof(int));

  return dst;
}

void
freeDA(darrays d
       )
{
  free(d->A);
  free(d);
}

void
pushDA(darrays d,
       int i
       )
{
  if(d->a == d->n){
    d->a *= 2;
    d->A = (int*)realloc(d->A, d->a*sizeof(int));
  }

  d->A[d->n] = i;
  d->n++;
}

/* Dump the array into a buffer */
void
dumpDA(darrays d,
       int *B
       )
{
  memcpy(B, d->A, d->n*sizeof(int));
}

void
expandDA(darrays d,
	 int m
	 )
{
  if(m > d->a){
    d->a *= 2;
    if(m > d->a)
      d->a = m;

    d->A = (int*)realloc(d->A, d->a*sizeof(int));
  }
}

void
resetDA(darrays d
	)
{
  d->n = 0;
}

void
resetIterator(darrays d
	      )
{
  d->i = 0;
}

int
hasNextDA(darrays d
	  )
{
  return d->i < d->n;
}

int
getNextDA(darrays d
	  )
{
  int r = -1;

  if(d->i < d->n){
    r = d->A[d->i];
    d->i++;
  }

  return r;
}
