#ifndef _SPLAY_TREE_H
#define _SPLAY_TREE_H

typedef struct sTree *sTrees;
typedef struct node *nodes;

sTrees
allocSTree(int n
           );

sTrees
dupSTree(sTrees t
        );

sTrees
cpySTree(sTrees dst,
         sTrees src
        );

void
freeSTree(sTrees t
          );

/* Turns into an empty tree */
void
clearTree(sTrees t
          );

int
size(sTrees t
     );

enum childI {left, right};

/* Return the key associated with this node */
int
key(sTrees t,
    nodes n
    );

/* Return the value stored in this node. */
int
value(sTrees t,
      nodes n,
      int *cv
      );

/* Return node pointer by index */
nodes
getNode(sTrees t,
        int i
        );

/* Return node pointer in order */
nodes
getNodeInOrder(sTrees t,
	       int i
	       );

/* Get the inoder keys */
int *
getInorder(sTrees t,
           int *L,
	   int *map,
	   enum childI d
           );

/* Insert by key */
void
insertKey(sTrees t,
	  int k
	  );

/* Insert inorder */
void
insertInorderKey(sTrees t,
		 int p, /* position */
		 int k
		 );

/* Insert by relative position */
void
insertN(sTrees t,
        nodes n, /* Reference node */
        int i, /* Node index on t */
        enum childI d /* Direction of reference */
        );

/* Make the node corresponding to i the new root */
void
reRoot(sTrees t,
       int i, /* New root index */
       enum childI d /* Direction of reference */
       );

void
removeN(sTrees t,
        nodes n /* Node to remove */
        );

void
splitSt(sTrees t,
	int k, /* The key */
	enum childI d /* Direction to keep */
	);

/* Next or prev element, depending on direction. */
nodes
diressor(sTrees t,
	 nodes n, /* Reference node */
	 enum childI d
	 );

/* Maximum or minimum, depending on direction. */
nodes
dirum(sTrees t,
      enum childI d
      );

/* Get the floor on key space */
void
roundSt(sTrees t,
        int k, /* The key */
	nodes *floor,
	nodes *ceil
        );

#endif /* _SPLAY_TREE_H */
