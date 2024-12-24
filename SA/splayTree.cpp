#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "splayTree.h"

struct node{
  int s; /* Number of nodes in sub-tree */
  nodes c[2]; /* child array */
  nodes *hook; /* How is this node connected? */
};

struct sTree{
  int n; /* sizeof N */
  nodes N; /* Node array */
  nodes root; /* Root of the Splay tree */
};

sTrees
allocSTree(int n
           )
{
  sTrees t = (sTrees)calloc(1, sizeof(struct sTree));
  t->N = (nodes)calloc(n, sizeof(struct node));
  t->n = n;

  return t;
}

sTrees
dupSTree(sTrees t
        )
{
  sTrees dup = (sTrees)calloc(1, sizeof(struct sTree));
  return cpySTree(dup, t);
}

sTrees
cpySTree(sTrees dst,
         sTrees src
        )
{
  if (!dst) return dupSTree(src);

  if (dst->n != src->n) {
    dst->n = src->n;
    // printf("realloc(%p, %d*%d): ", (void*)dst->N, src->n, sizeof(struct node));
    dst->N = (nodes)realloc(dst->N, src->n*sizeof(struct node));
    // printf("%p\n", (void*)dst->N);
  }

  /* fix pointers */
  if (src->root) {
    dst->root = dst->N + (src->root - src->N);
  } else {
    dst->root = NULL;
  }
  nodes srcN = src->N, dstN = dst->N; 
  for (int i = 0; srcN != src->N + src->n/*i < src->n*/; ++i) {
    if (srcN->hook) {
      if (srcN->hook == &(src->root)) { // hook points the root
        dstN->hook = &(dst->root);
        // printf("root hook store  (src=%p dst=%p)\n", (void*)srcN->hook, (void*)dstN->hook);
      } else { // hook points some child node
        dstN->hook = (nodes*) ((long unsigned)dst->N + ((long unsigned)(srcN->hook) - (long unsigned)src->N));
        // printf(" srcN[%i]->c (src=%p dst=%p) [src:%p-%p, dst:%p-%p]\n", i,
        //   (void*)srcN->hook, (void*)dstN->hook,
        //   (void*)&src->N[0], (void*)&src->N[src->n-1],
        //   (void*)&dst->N[0], (void*)&dst->N[dst->n-1]);
      }
      if (*srcN->hook) {
        *dstN->hook = dst->N + (*srcN->hook - src->N);
        // printf("  - hook (src=%p dst=%p) [src:%p-%p, dst:%p-%p]\n",
        //   (void*)*srcN->hook, (void*)*dstN->hook,
        //   (void*)&src->N[0], (void*)&src->N[src->n-1],
        //   (void*)&dst->N[0], (void*)&dst->N[dst->n-1]);
      }
    }
    else dstN->hook = NULL;
    if (srcN->c[left]) dstN->c[left] = dst->N + (srcN->c[left] - src->N);
    else               dstN->c[left] = NULL;
    if (srcN->c[right]) dstN->c[right] = dst->N + (srcN->c[right] - src->N);
    else                dstN->c[right] = NULL;

    dstN->s = srcN->s;
    srcN++;
    dstN++;
  }

  return dst;
}

void
freeSTree(sTrees t
          )
{
  free(t->N);
  free(t);
}

static void
clearTreeR(nodes n
	   )
{
  if(NULL != n){
    if(n->hook) *(n->hook) = NULL;
    n->hook = NULL;
    clearTreeR(n->c[left]);
    clearTreeR(n->c[right]);
  }
}

/* Turns into an empty tree */
void
clearTree(sTrees t
          )
{
  clearTreeR(t->root);
}

static void
clearNode(sTrees t,
	  int i
          )
{
  bzero(&t->N[i], sizeof(struct node));
}

int
size(sTrees t
     )
{
  int r = 0;
  if(NULL != t->root)
    r = t->root->s;

  return r;
}

static nodes
father(nodes n
       )
{
  assert(NULL != n->hook && "Father call on root node.");

  long int d = (long int) n->hook - (long int) n;

  if(0 > d){
    d -= sizeof(struct node);
    d++;
  }

  nodes f = &n[d/(int)sizeof(struct node)];

  return f;
}

static enum childI
childIdx(nodes n, /* The child */
	 nodes f /* The father */
	 )
{
  enum childI i = left;

  if(f->c[right] == n)
    i = right;

  return i;
}

static enum childI
other(enum childI dir
      )
{
  enum childI r = left;
  if(r == dir)
    r = right;

  return r;
}

int
key(sTrees t,
    nodes n
    )
{
  return ((long unsigned int)n - (long unsigned int)t->N)/sizeof(struct node);
}

static int *
getInorderR(sTrees t,
            nodes n,
            int *L,
            int *map,
            enum childI d
            )
{
  if(NULL != n->c[d])
    L = getInorderR(t, n->c[d], L, map, d);

  if(NULL == map)
    *L = key(t, n);
  else
    *L = map[key(t, n)];
  L++;

  if(NULL != n->c[other(d)])
    L = getInorderR(t, n->c[other(d)], L, map, d);

  return L;
}

int *
getInorder(sTrees t,
           int *L,
	         int *map,
	         enum childI d
          )
{
  int *F = L;

  if(NULL != t->root)
    F = getInorderR(t, t->root, L, map, d);

  return F;
}

static void
resetNode(nodes n
          )
{
  *(n->hook) = n;

  n->s = 1;
  if(NULL != n->c[left])
    n->s += n->c[left]->s;
  if(NULL != n->c[right])
    n->s += n->c[right]->s;
}

#ifndef NDEBUG
static void
checkTree(sTrees t
          )
{
  for(int i = 0; i < t->n/*(int)(sizeof(t->N)/sizeof(struct node))*/; i++){
    if(NULL != t->N[i].hook) {
      // printf("check %p father (%p)\n", &(t->N[i]), father(&(t->N[i])));
      assert(&(t->N[i]) != father(&(t->N[i])) && "Found loop in tree");
    }
  }
}
#endif /* NDEBUG */

static void
splay(sTrees t,
      nodes n
      )
{ /* WARNING: Apply only to nodes in the tree */
  assert(NULL != n->hook && "Splay call on out of tree node");
#ifndef NDEBUG
  checkTree(t);
#endif /* NDEBUG */

  while(&(t->root) != n->hook){
    nodes f = father(n); /* Father node */
    assert(f != n && "Strange loop in tree");
    enum childI fi = childIdx(n, f);
    if(&(t->root) == f->hook){
      /* Simple Zig transformation
         |
         |     f    =>    n
         |    / \   =>   / \
         |   n   C  =>  A   f
         |  / \     =>     / \
         | A   B    =>    B   C
         |
         | fi = /   ;  other(fi) = \   */

      n->hook = f->hook; /* Copy hook */
      f->c[fi] = n->c[other(fi)];
      if(NULL != f->c[fi])
        f->c[fi]->hook = &f->c[fi];
      f->hook = &n->c[other(fi)];
    } else {
      nodes gf = father(f); /* Grand Father node */
      assert(gf != n && gf != f && "Weird loop in tree");
      enum childI gfi = childIdx(f, gf);

      if(fi == gfi){ /* Zig Zig */
        /* Zig Zig transformation
          |
          |       gf   =>    n
          |      /  \  =>   / \
          |     f    D =>  A   f
          |    / \     =>     / \
          |   n   C    =>    B   gf
          |  / \       =>       /  \
          | A   B      =>      C    D
          |
          | fi = /   ;  other(fi) = \    */

        n->hook = gf->hook;
        f->hook = &n->c[other(fi)];
        gf->hook = &f->c[other(fi)];

        f->c[fi] = n->c[other(fi)];
        if(NULL != f->c[fi])
          f->c[fi]->hook = &f->c[fi];
        gf->c[fi] = f->c[other(fi)];
        if(NULL != gf->c[fi])
          gf->c[fi]->hook = &gf->c[fi];
      } else { /* Zig Zag */
        /* Zig Zag transformation
          |
          |     gf    =>       n
          |    /  \   =>      / \
          |   f    D  =>     /   \
          |  / \      =>    /     \
          | A   n     =>   f       gf
          |    / \    =>  / \     /  \
          |   B   C   => A   B   C    D
          |
          | fi = \   ;  gfi = /     */

        n->hook = gf->hook;
        f->hook = &n->c[gfi];
        gf->hook = &n->c[fi];

        f->c[fi] = n->c[gfi];
        if(NULL != f->c[fi])
          f->c[fi]->hook = &f->c[fi];
        gf->c[gfi] = n->c[fi];
        if(NULL != gf->c[gfi])
          gf->c[gfi]->hook = &gf->c[gfi];
      }
      resetNode(gf);
    }
    resetNode(f);
    resetNode(n);
  }
  resetNode(n);

#ifndef NDEBUG
  checkTree(t);
#endif /* NDEBUG */
}

/* Return value stored in this node. */
int
value(sTrees t,
      nodes n,
      int *cv /* complementar value */
      )
{ /* WARNING: Apply only to nodes in the tree */
  assert(NULL != n->hook && "Splay call on out of tree node");

  int r = 0;
  splay(t, n);
  if(NULL != n->c[left])
    r = n->c[left]->s;

  if(NULL != cv){
    *cv = 0;
    if(NULL != n->c[right])
      *cv = n->c[right]->s;
  }

  return r;
}

/* Return node pointer by index */
nodes
getNode(sTrees t,
        int i
        )
{
  nodes r = NULL;
  if(NULL != t->N[i].hook)
    r = &(t->N[i]);

  return r;
}

nodes
getNodeInOrder(sTrees t,
	       int i
	       )
{
  nodes r = t->root;

  assert(NULL != r && "Request node on empty tree.");

  int l = 0; /* Number of elements on the left */

  while(NULL != r && l != i){
    if(NULL != r->c[left])
      l += r->c[left]->s;

    if(l != i){
      enum childI d = left;
      if(l < i)
	      d = right;

      r = r->c[d];
    }
  }

  return r;
}

void
insertKey(sTrees t,
	  int k
	  )
{
  assert(NULL == t->N[k].hook && "Inserting existing node");
  // if (NULL != t->N[k].hook) { printf("duplicate!\n"); return; }

  if(NULL == t->root){ /* Insert on empty tree */
    reRoot(t, k, right);
  } else {
    nodes floor;
    nodes ceil;
    roundSt(t, k, &floor, &ceil);

    if(NULL == floor)
      reRoot(t, k, right);
    else if(NULL == ceil)
      reRoot(t, k, left);
    else
      insertN(t, floor, k, right);
  }
}

void
insertInorderKey(sTrees t,
		 int p, /* position */
		 int k
		 )
{
  assert(NULL == t->N[k].hook && "Inserting existing node");

  nodes *h = &t->root;
  nodes u = *h;

  while(NULL != u){
    int lsize = 0;
    if(NULL != u->c[left])
      lsize = u->c[left]->s;

    if(lsize < p){
      p -= 1 + lsize;
      h = &u->c[right];
    } else {
      h = &u->c[left];
    }
    u = *h;
  }

  clearNode(t, k);
  nodes v = &t->N[k];
  v->hook = h;
  *h = v;
  splay(t, v);
}


/* Insert this node into the tree */
void
insertN(sTrees t,
        nodes n, /* Reference node */
        int i, /* Node to insert */
        enum childI d /* Direction of reference */
        )
{ /* WARNING: Apply only to nodes in the tree */
  assert(NULL != n->hook && "Splay call on out of tree node");
  assert(NULL == t->N[i].hook && "Inserting existing node");

  if(NULL != n->c[d]){
    n = n->c[d];
    d = other(d);
    while(NULL != n->c[d])
      n = n->c[d];
  }

  clearNode(t, i);
  nodes v = &t->N[i];
  v->hook = &n->c[d];
  *(v->hook) = v;
  splay(t, v);
}

/* Make the node corresponding to i the new root */
void
reRoot(sTrees t,
       int i, /* New root index */
       enum childI d /* Direction of reference */
       )
{
  assert(NULL == t->N[i].hook && "Rooting to active node");

  clearNode(t, i);

  t->N[i].c[d] = t->root;
  if(NULL != t->N[i].c[d])
    t->N[i].c[d]->hook = &(t->N[i].c[d]);

  t->N[i].hook = &t->root;
  resetNode(&t->N[i]);
}

void
removeN(sTrees t,
        nodes n /* Node to remove */
        )
{ /* WARNING: Apply only to nodes in the tree */
  assert(NULL != n->hook && "Splay call on out of tree node");

  if(1 == t->root->s) /* Tree has only 1 node */
    t->root = NULL;
  else {  /* Tree has several nodes */
    splay(t, n);

    nodes o = n->c[left];
    if(NULL == o)
      o = n->c[right];
    else if(NULL != n->c[right]) { /* General case */
      while(NULL != o->c[right])
        o = o->c[right];

      splay(t, o); /* Splay the other node */
      splay(t, n); /* Splay n again */

      /* Removal transform
        |
        |     n    =>
        |    / \   =>
        |   o   B  =>    o
        |  /       =>   / \
        | A        =>  A   B
        |
        | d = \     ; other(d) = /     */

      o->c[right] = n->c[right];
      o->c[right]->hook = &o->c[right];
    }

    o->hook = &t->root;
    resetNode(o);
  }
  n->hook = NULL;
}

void
splitSt(sTrees t,
	int k, /* The key */
	enum childI d /* Direction */
	)
{
  nodes floor;
  nodes ceil;
  roundSt(t, k, &floor, &ceil);

  switch(d){
  case left:
    if(NULL == floor)
      clearTree(t);
    else{
      splay(t, floor);
      clearTreeR(floor->c[right]);
    }
    break;

  case right:
    if(NULL == ceil)
      clearTree(t);
    else{
      splay(t, ceil);
      clearTreeR(ceil->c[left]);
    }
    break;

  default:
    assert(0 && "Invalid split direction ");
    break;
  }
}

/* Extreme movement */
static nodes
extreme(sTrees t,
	nodes n,
	enum childI d
	)
{
  assert(NULL != n && "extreme call on NULL node.");

  while(NULL != n->c[d])
    n = n->c[d];

  splay(t, n);

  return n;
}

nodes
diressor(sTrees t,
	 nodes n, /* Reference node */
	 enum childI d
	 )
{
  assert(NULL != n && "diressor call on NULL node.");

  splay(t, n); /* Make it into root */
  if(NULL != n->c[d]){
    n = n->c[d];

    n = extreme(t, n, other(d));
  }

  return n;
}

nodes
dirum(sTrees t,
      enum childI d
      )
{
  nodes r = t->root;

  if(NULL != r)
    r = extreme(t, r, d);

  return r;
}

void
roundSt(sTrees t,
        int k, /* The key */
	nodes *floor,
	nodes *ceil
        )
{
  nodes n = t->root;
  *floor = NULL;
  *ceil = NULL;
  nodes last = NULL;

  while(NULL != n && key(t, n) != k){
    last = n;

    if(key(t, n) < k){
      *floor = n;
      n = n->c[right];
    } else {
      *ceil = n;
      n = n->c[left];
    }
  }

  if(key(t, n) == k){
    last = n;
    *floor = n;
    *ceil = n;
  }

  if(NULL != last)
    splay(t, last);
}
