#ifndef _REORGANIZE_H
#define _REORGANIZE_H

#include "graph.h"
#include "sa_state.h"

typedef struct organizer *organizers;

/* Creates a new organizing structure */
organizers
allocO(graph G,
       sTrees t
       );

organizers
dupO(organizers o
      );

organizers
cpyO(organizers dst,
     organizers src
    );

void
freeO(organizers o
      );

/* Checks to see if the heuristic amortizes.
   Otherwise increment internal counter. */
int
checkO(organizers o,
       int dE
       );

/* Inserts topological order into t */
void
reorganize(organizers o,
	   int *L /* List of break nodes */
	   );

#endif /* _REORGANIZE_H */
