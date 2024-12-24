#include <stdlib.h>
#include <string.h>

#include "queue.h"

struct queue{
  int *A; /* Array of elements */
  int front;
  int back;
  int n; /* size of the queue */
};

/* Create a queue that supports n elements */
queues
allocQ(int n
       )
{
  queues q = (queues)malloc(sizeof(struct queue));

  q->n = n+1;
  q->A = (int*) malloc(q->n*sizeof(int));
  q->front = 0;
  q->back = 0;

  return q;
}

queues
dupQ(queues q
    )
{
  queues dup = (queues)calloc(1, sizeof(struct queue));
  return cpyQ(dup, q);
}

queues
cpyQ(queues dst,
     queues src
    )
{
  if (!dst) return dupQ(src);
  dst->n = src->n;
  dst->back = src->back;
  dst->front = src->front;
  dst->A = (int*)realloc(dst->A, src->n*sizeof(int));
  memcpy(dst->A, src->A, src->n*sizeof(int));
  return dst;
}

void
freeQ(queues q
      )
{
  free(q->A);
  free(q);
}

/* True iff queue is empty */
int
emptyQ(queues q
       )
{
  return q->front == q->back;
}

/* Current element in front */
int
frontQ(queues q
       )
{
  return q->A[q->front];
}

/* Put this element at the end of the queue */
void
pushQ(queues q,
      int e /* element to push to queue */
      )
{
  q->A[q->back] = e;
  q->back++;
  q->back %= q->n;
}

/* Remove the element in front of the queue */
void
popQ(queues q
     )
{
  q->front++;
  q->front %= q->n;
}

