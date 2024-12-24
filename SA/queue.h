#ifndef _QUEUE_H
#define _QUEUE_H

typedef struct queue *queues;

/* Create a queue that supports n elements */
queues
allocQ(int n
       );

queues
dupQ(queues q
    );

queues
cpyQ(queues dst,
     queues src
    );

void
freeQ(queues q
      );

/* True iff queue is empty */
int
emptyQ(queues q
       );


/* Current element in front */
int
frontQ(queues q
       );

/* Put this element at the end of the queue */
void
pushQ(queues q,
      int e /* element to push to queue */
      );

/* Remove the element in front of the queue */
void
popQ(queues q
     );

#endif /* _QUEUE_H */
