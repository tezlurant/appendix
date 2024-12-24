#ifndef _DARRAY_H
#define _DARRAY_H

typedef struct darray *darrays;

darrays
allocDA(void
	);

darrays
dupDA(darrays d
	);

darrays
cpyDA(darrays dst,
      darrays src
	   );

void
freeDA(darrays d
       );

void
pushDA(darrays d,
       int i
       );

/* Dump the array into a buffer */
void
dumpDA(darrays d,
       int *B
       );

void
expandDA(darrays d,
	 int m
	 );

void
resetDA(darrays d
	);

void
resetIterator(darrays d
	      );

int
hasNextDA(darrays d
	  );

int
getNextDA(darrays d
	  );

#endif /* _DARRAY_H */
