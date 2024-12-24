#include <stdlib.h>
#include <bsd/stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#include "graph.h"
#include "splayTree.h"
#include "splayTree.cpp"
#include "darray.h"
#include "darray.cpp"
#include "reorganize.h"
#include "reorganize.cpp"
#include "sa_state.h"

#define NDEBUG

#ifdef POSTSCRIPT
FILE *LOG = NULL;
int POST_NB_VERT = 0;
#endif /* POSTSCRIPT */

#ifndef NDEBUG
static void
assertState(sa_state s);
#endif /* NDEBUG */

//#define ARC4RAND_UNIF(_num) arc4random_uniform(_num)
#define ARC4RAND_UNIF(_num) random() % _num

unsigned long nbProcVert = 0;



struct sa_state_ {
  graph G; /* The underlying graph */
  sTrees t; /* K = vertex index. V = Topological Order index. */
  sTrees O[2]; /* K = Topological index. V = sub counter. */
  int *Top2V; /* Maps topological order to vertex. */
  int *condemned; /* List of nodes to remove. Normaly small */

  /* Candidate List stuff */
  // char *cBool; /* Candidate booleans */
  int cLs; /* Size of candidate list */
  int *cList; /* List of candidates */
  int reorderCount;
  /* The current candidate is cList[0] */

  /* Excluded adjacency in check */
  darrays DegEx[2]; /* Excluded vertexes */

  /* Position lives in the topological order */
  int p; /* The position selected by check */

  enum adjacencies dir; /* Direction selected by check */
  /* Delta E for last check. Used to trim execute. */
  int checkDE;

  /* Supporting the reorganizing heuristic */
  organizers o;
};

sa_state
allocS(graph G)
{
  sa_state s = (sa_state) calloc(1, sizeof(struct sa_state_));

  s->G = G;
  s->t = allocSTree(1+G->v); /* Add a sentinel */
  s->O[in] = allocSTree(G->v);
  s->O[out] = allocSTree(G->v);
  s->Top2V = (int *)calloc(G->v, sizeof(int));
  s->condemned = (int *)calloc(G->v, sizeof(int));
  // s->cBool = (char *)calloc(G->v, sizeof(char));
  s->reorderCount = 0;
  s->cList = (int *)malloc(G->v*sizeof(int));

  s->DegEx[in] = allocDA();
  s->DegEx[out] = allocDA();

  s->o = allocO(G, s->t);

  return s;
}

sa_state // duplicates the sa_state
dupS(sa_state s)
{
  sa_state dup = (sa_state) calloc(1, sizeof(struct sa_state_));
  return cpyS(dup, s);
}

sa_state
cpyS(sa_state dst,
     sa_state src
    )
{
  // printf("DST %p sa_state\n", dst);
  // printS(dst);
  // printf("SRC %p sa_state\n", src);
  // printS(src);
  if (!dst) return dupS(src);

#ifndef NDEBUG
  assertState(src);
#endif /* NDEBUG */

  dst->G       = src->G;
  dst->cLs     = src->cLs;
  dst->checkDE = src->checkDE;
  dst->p       = src->p;
  dst->dir     = src->dir;

  dst->t = cpySTree(dst->t, src->t); /* Add a sentinel */
  dst->O[in] = cpySTree(dst->O[in], src->O[in]);
  dst->O[out] = cpySTree(dst->O[out], src->O[out]);

  dst->Top2V = (int *)realloc(dst->Top2V, src->G->v*sizeof(int));
  memcpy(dst->Top2V, src->Top2V, src->G->v*sizeof(int));

  dst->condemned = (int *)realloc(dst->condemned, src->G->v*sizeof(int));
  memcpy(dst->condemned, src->condemned, src->G->v*sizeof(int));

  // dst->cBool = (char *)realloc(dst->cBool, src->G->v*sizeof(char));
  // memcpy(dst->cBool, src->cBool, src->G->v*sizeof(char));

  dst->cList = (int *)realloc(dst->cList, src->G->v*sizeof(int));
  memcpy(dst->cList, src->cList, src->G->v*sizeof(int));

  dst->DegEx[in] = cpyDA(dst->DegEx[in], src->DegEx[in]);
  dst->DegEx[out] = cpyDA(dst->DegEx[out], src->DegEx[out]);

  dst->o = cpyO(dst->o, src->o);

#ifndef NDEBUG
  assertState(dst);
#endif /* NDEBUG */

  // printf("dst: %p\n", (void*)dst);
  // printS(dst);
  // printf("src: %p\n", (void*)src);
  // printS(src);

  return dst;
}

void
prepareS(sa_state s,
         int*  A, /* Array with vertexes */
         int   l	/* length of the array */
) {
  s->cLs = 0;
  clearTree(s->t);

#ifndef NUSE_INIT_APPROX
  /* Dump array into the tree */
  for (int i = 0; i < l; i++) {
    reRoot(s->t, A[i], left);
    // s->cBool[A[i]] = 1;

#ifdef POSTSCRIPT
    fprintf(LOG, "%d %d i\n",
	    1+(A[i]%POST_NB_VERT),
	    1+(A[i]/POST_NB_VERT));
    fflush(LOG);
#endif /* POSTSCRIPT */
  }

  /* Run approximation discarding */
  for(int i = 0; i < l; i++){
    int v = A[i];
    /* Identify position in topological order */
    nodes nv = getNode(s->t, v);
    if(NULL != nv){
      int top = value(s->t, nv, NULL);
      int *outA = s->G->E[v][out];

      nodes u;
      int valid = 1;
      while (-1 != *outA && valid) {
        // s->cBool[*outA] = 1;
        u = getNode(s->t, *outA);
        if (NULL != u &&
            value(s->t, u, NULL) <= top
        ) valid = 0;
        outA++;
        nbProcVert++;
      }

      if (!valid) { /* Means early temination */
        removeN(s->t, u);
        s->cList[s->cLs] = outA[-1];
        s->cLs++;
#ifdef POSTSCRIPT
	fprintf(LOG, "%d %d c\n",
		1+(outA[-1]%POST_NB_VERT),
		1+(outA[-1]/POST_NB_VERT)); // TODO
	fflush(LOG);
#endif /* POSTSCRIPT */
        removeN(s->t, nv);
        s->cList[s->cLs] = v;
        s->cLs++;
#ifdef POSTSCRIPT
	fprintf(LOG, "%d %d c\n",
		1+(v%POST_NB_VERT),
		1+(v/POST_NB_VERT));
	fflush(LOG);
#endif /* POSTSCRIPT */
      }
    }
  }
#else
  for (int i = 0; i < l; i++) {
    s->cList[s->cLs] = A[i];
    s->cLs++;
  }
#endif /* NUSE_INIT_APPROX */
}

void
prepareS2(sa_state s,
          int* A, /* Array with vertexes */
          int l,	/* length of the array */
          int *B, /* Array with vertexes */
          int lB /* length of the array */
)
{
  s->cLs = 0;
  clearTree(s->t);
  for (int i = 0; i < lB; i++) 
  {
    reRoot(s->t, B[i], left);
  }
  int n = SA_get_nbVerts();
  for(int i = lB; i < n; i++)
  {
    s->cList[s->cLs] = B[i];
    s->cLs++;
  }
}

void
prepareSwithKnowSolution(sa_state s,
         int*  A, /* Array with vertexes */
         int   l,	/* length of the array */
         int*  C, /* Array all vertexes in SCC */
         int  sC /* size of the SCC */
) {
  s->cLs = 0;
  clearTree(s->t);
  /* Dump array into the tree */
  for (int i = 0; i < l; i++) {
    reRoot(s->t, A[i], left);
  }
  for (int i = 0; i < sC; i++) {
    if (NULL == getNode(s->t, C[i])) {
      s->cList[s->cLs] = C[i];
      s->cLs++;
    }
  }
}

void
freeS(sa_state s)
{
  freeSTree(s->t);
  freeSTree(s->O[in]);
  freeSTree(s->O[out]);
  free(s->Top2V);
  free(s->condemned);
  // free(s->cBool);
  free(s->cList);

  freeDA(s->DegEx[in]);
  freeDA(s->DegEx[out]);
  freeO(s->o);

  free(s);
}

int
getE(sa_state s)
{
  return size(s->t);
}

#ifndef NDEBUG
static void
gdbBreak(void)
{}

static void
assertState(sa_state s)
{
  for(int i = 0; i < s->G->v; i++){
    nodes u = getNode(s->t, i);
    if(NULL != u){
      int posU = value(s->t, u, NULL);
      int *outA = s->G->E[i][out];
      while(-1 != *outA){
        nodes v = getNode(s->t, *outA);
        if(NULL != v){
          int posV = value(s->t, v, NULL);
          // printf("test edge %i(pos=%i)->%i(pos=%i) \n", i, posU, *outA, posV);
          if(posU >= posV)
            gdbBreak();
          assert(posU < posV && "Invalid sa_state configuration");
        }
        outA++;
      }
    }
  }
}
#endif /* NDEBUG */

/* Modified sa_state print for torus */
void
printS(sa_state s)
{
  int side = ceil(sqrt(s->G->v));
  int A[side*side];

  for(int i = 0; i < side*side; i++)
    A[i] = -1;

  if(NULL != s->t){
    int L[1+size(s->t)];
    *getInorder(s->t, &L[0], NULL, left) = -1;

    int j = 0;
    int *T = L;
    while(-1 != *T){
      A[*T] = j++;
      T++;
    }
  }

  int v = s->cList[0];
  for(int i = 0; i < side*side; i++){
    if(v == i) {
      fprintf(stderr, " * ");
    } else if(-1 != A[i]) {
      fprintf(stderr, "%.2d ", A[i]);
    } else {
      fprintf(stderr, "   ");
    }
    if(0 == ((1+i)%side)){
      fprintf(stderr, "\n");
    }
  }

  /* for(int j = 0; j<20; j++) */
    /* fprintf(stderr, "\n"); */
}

static void
swapI(int *A,
      int *B
      )
{
  static int T;

  T = *A;
  *A = *B;
  *B = T;
}

int 
swapSelect(sa_state s,
        int p
        )
{
  //assert(0 < s->cLs && "Empty candidate list");
  if(0 > p) return -1;
  int v = s->cList[0];
  int u = s->cList[p];
  int iterV = iterCList[v];
  iterCList[v] = iterCList[u];
  iterCList[u] = iterV;
  swapI(&(s->cList[0]), &(s->cList[p]));
  return 0;
}

int
selectV(sa_state s)
{
  return s->cList[0];
}

int
swapSelectSE(sa_state s,
        int p,
        int q
        )
{
  //assert(0 < s->cLs && "Empty candidate list");
  if(0 > p || 0 > q) return -1;
  swapI(&(s->cList[p]), &(s->cList[q]));
  return 0;
}

// ここで次に非巡回グラフに追加する候補を選ぶ
// int swapIdxの部分で一様に選ぶ
// この確率を調整する

int
choose(sa_state s
       )
{
  assert(0 < s->cLs && "Empty candidate list");
  if (0 >= s->cLs) return -1;
  int swapIdx = ARC4RAND_UNIF(s->cLs);
  // printf("choose(s) cList[0](%i)<-->cList[%i](%i)\n",
  //   s->cList[0], swapIdx, s->cList[swapIdx]);
  int v = s->cList[0];
  int u = s->cList[swapIdx];
  int iterV = iterCList[v];
  iterCList[v] = iterCList[u];
  iterCList[u] = iterV;
  swapI(&(s->cList[0]), &(s->cList[swapIdx]));

  // printf("-cList (size t = %i): ", size(s->t));
  // for (int i = 0; i < s->cLs; ++i) printf("%i ", s->cList[i]);
  // printf("\n");
  return 0;
}

int
chooseSB(sa_state s,
        int selectBorder
       )
{
  assert(0 < selectBorder && "Empty candidate list");
  if (0 >= selectBorder) return -1;
  int swapIdx = ARC4RAND_UNIF(selectBorder);
  // printf("choose(s) cList[0](%i)<-->cList[%i](%i)\n",
  //   s->cList[0], swapIdx, s->cList[swapIdx]);
  swapI(&(s->cList[0]), &(s->cList[swapIdx]));
  // printf("-cList (size t = %i): ", size(s->t));
  // for (int i = 0; i < s->cLs; ++i) printf("%i ", s->cList[i]);
  // printf("\n");
  return 0;
}

int
chooseLogdV(sa_state s
      )
{
  assert(0 < s->cLs && "Empty candidate list");
  if(cLsSize(s) == 0) return 0;
  if(cLsSize(s) == 1) return choose(s);
  int logCLss = (int)log2(cLsSize(s) + 1); //候補の数(+1は0の場合の対策)
  if(logCLss == 0) 
  {
    logCLss = 1;
  }
  choose(s);
  double lambda = 0.3;
  int v = s->cList[0];
  double inDegV = s->G->inDeg[v];
  double outDegV = s->G->outDeg[v];
  double p = inDegV + outDegV - lambda*abs(inDegV - outDegV);
  //double p = sqrt((inDegV-1)*(outDegV-1));
  logCLss--;
  while(logCLss--)
  {
    int swapIdxC = ARC4RAND_UNIF(s->cLs - 1);
    swapIdxC++;
    int u = s->cList[swapIdxC];
    double inDegU = s->G->inDeg[u];
    double outDegU = s->G->outDeg[u];
    double q = inDegU + outDegU - lambda*abs(inDegU - outDegU);
    //double q = sqrt((inDegU-1)*(outDegU-1));
    if(q < p)
    {
      swapI(&(s->cList[0]), &(s->cList[swapIdxC]));
      inDegV = inDegU;
      outDegV = outDegU;
      v = u;
      p = q;
    }
  }
  return 0;
}

int
chooseAssign(sa_state s,
     int cLss          
     )
{
  assert(0 < s->cLs && "Empty candidate list");
  if(cLss == 0) return 0;
  int swapIdx = ARC4RAND_UNIF(cLss);
  swapI(&(s->cList[0]), &(s->cList[swapIdx]));
  return 0;
}

int
cLsSize(sa_state s
     )
{
  return s->cLs;
}

void
cLsZero(sa_state s
     )
{
  s->cLs = 0;
}

int
cLsAdd(sa_state s,
     int v
     )
{
  if (0 > s->cLs) return -1;
  iterCList[v] = s->cLs;
  s->cList[s->cLs] = v;
  s->cLs++;
  if(s->cLs > s->G->v) return -1;
  return 0;
}

int
cLsGet(sa_state s,
     int i
     )
{
  return s->cList[i];
}

/* Returns how many elements will be removed if position
   p is used. */
static int
incompatibleSize(sa_state s,
		 int p
		 )
{
  int r = 0;
  nodes floor;
  nodes ceil;

  roundSt(s->O[out], p, &floor, &ceil);
  if(NULL != ceil){
    r += value(s->O[out], ceil, NULL);
  } else
    r += size(s->O[out]);
  /* Counting strictly less than p */

  roundSt(s->O[in], p, &floor, &ceil);
  if(NULL != floor){
    int c;
    value(s->O[in], floor, &c);
    r += c;

    if(floor == ceil){ /* Also equal to p */
      r++;
    }
  } else
    r += size(s->O[in]);
  /* Counting less than or equal to p */

  return r;
}

static int
next(sa_state s,
     int pos,
     int bound
     )
{
  int p = INT_MAX;
  /* int pout = INT_MAX; */

  nodes floor;
  nodes ceil;

  nbProcVert++;

  roundSt(s->O[in], pos, &floor, &ceil);
  if(NULL != ceil)
    p = 1 + key(s->O[in], ceil);

  /* roundSt(s->O[out], pos, &floor, &ceil); */
  /* if(NULL != ceil) */
  /*   pout = 1 + key(s->O[out], ceil); */

  /* if(pout < p) */
  /*   p = pout; */

  if(pos < bound && p > bound)
    p = bound;

  return p;
}

static int
prev(sa_state s,
     int pos,
     int bound
     )
{
  int p = -1;
  int pout = -1;

  nodes floor;
  nodes ceil;

  nbProcVert++;

  /* roundSt(s->O[in], pos-1, &floor, &ceil); */
  /* if(NULL != floor) */
  /*   p = key(s->O[in], floor); */

  roundSt(s->O[out], pos-1, &floor, &ceil);
  if(NULL != floor)
    pout = key(s->O[out], floor);

  if(pout > p)
    p = pout;

  if(pos > bound && p < bound)
    p = bound;

  return p;
}

void
neighborInitPdELog(sa_state s,
                  int v
                  )
{
  const int INF = 0x5fffffff;
  repushL = 0;
  for(int i=0;i<s->G->outDeg[v];i++){
    int to = s->G->E[v][out][i];
    pdELog[to] = INF;
    if(iterCList2[to] == -1)
    {
      repushList[repushL] = to;
    }
  }
  for(int i=0;i<s->G->inDeg[v];i++){
    int to = s->G->E[v][in][i];
    pdELog[to] = INF;
    if(iterCList2[to] == -1)
    {
      repushList[repushL] = to;
    }
  }
  pdELog[v] = INF;
}

double H[MAX_V][2];

void mergeSort2(int left, int right)
{
	if (left + 1 < right)
	{
		int mid;
		mid = (left + right) / 2;
		// 分割
		mergeSort2(left, mid);
		mergeSort2(mid, right);
		// 結合
		auto merge = [&](auto &self, int left, int mid, int right) -> void
		{
			int i, j, k;
			int n1, n2; // 部分列L, Rの要素数決定に利用

			n1 = mid - left;
			n2 = right - mid;

			double **L = new double *[n1 + 1];
			for (int i = 0; i < n1 + 1; i++)
			{
				L[i] = new double[2];
			}
			double **R = new double *[n2 + 1];
			for (int i = 0; i < n2 + 1; i++)
			{
				R[i] = new double[2];
			}
			for (i = 0; i < n1; i++)
			{
				L[i][0] = H[left + i][0];
        L[i][1] = H[left + i][1];
			}
			for (i = 0; i < n2; i++)
			{
				R[i][0] = H[mid + i][0];
        R[i][1] = H[mid + i][1];
			}
			L[n1][0] = 1e9;
			R[n2][0] = 1e9;

			j = 0;
			k = 0;
			for (i = left; i < right; i++)
			{
				if (L[j][0] <= R[k][0])
				{
					H[i][0] = L[j][0];
          H[i][1] = L[j][1];
					j++;
				}
				else
				{
					H[i][0] = R[k][0];
          H[i][1] = R[k][1];
					k++;
				}
			}
		};
		merge(merge, left, mid, right);
	}
}

// 初期解を求める
void
initSolution(sa_state s)
{
  double lambda = 0.3;
  int clss = cLsSize(s);
  cLsZero(s);
  for(int i = 0; i < clss; i++)
  {
    int v = s->cList[i];
    H[i][1] = (double)v;
    H[i][0] = s->G->outDeg[v] + s->G->inDeg[v] - lambda*abs(s->G->outDeg[v] - s->G->inDeg[v]);
  }
  mergeSort2(0, clss);
  for(int i = 0; i < clss; i++)
  {
    int v = (int)H[i][1];
    //追加判定
    int dE = 1;
    if(dE <= checkB(s,dE,v))
    {
      insertV(s, v);
    }
    else
    {
      cLsAdd(s, v);
    }
  }
}

/* Returns the proposed energy variation */
// 先行研究では乱数のところが動くようにする
int
check(sa_state s,
      int dE
      )
{
  // printf("cList (size t = %i): ", size(s->t));
  // for (int i = 0; i < s->cLs; ++i) printf("%i ", s->cList[i]);
  // printf("\n");
  assert(2 > dE && "Check with impossible dE");
  
  if(checkO(s->o, 0)){
    /*s->reorderCount = 1;*/ // this fix does not work
    s->cList[s->cLs] = -1;
    reorganize(s->o, s->cList);
    const int INF = 0x5fffffff;
    for(int i=0;i<MAX_V;i++) pdELog[i] = INF;
  }
  
  /* Apply re-organization heuristic */
  #ifdef EXESTING_TECHNIQUE
  if(checkO(s->o, 0)){
    /*s->reorderCount = 1;*/ // this fix does not work
    s->cList[s->cLs] = -1;
    reorganize(s->o, s->cList);
  #ifndef NDEBUG
    assertState(s);
  #endif /* NDEBUG */
  }
  #endif

  /* Candidate */
  int v = s->cList[0];

  /* printf("State before check\n"); */
  /* printS(s); */
  /* printf(">>>> Candidate is : %d\n", v); */

  /* Load trees O[in] and O[out] */
  clearTree(s->O[in]);
  clearTree(s->O[out]);
  resetDA(s->DegEx[in]);
  resetDA(s->DegEx[out]);

  int *inA = s->G->E[v][in];
  int *outA = s->G->E[v][out];
  /* Positions over Topological order */
  int pa = 0;
  int pb = size(s->t);

  while(pa <= pb &&
        !(-1 == *inA && -1 == *outA)
	){
    /* Increment the heuristic */
    checkO(s->o, 2);

    /* Add one incoming */
    if(pa <= pb && -1 != *inA){
      if(NULL == getNode(s->t, *inA))
	      pushDA(s->DegEx[in], *inA);
      else { /* Vertex exists in DAG */
        int top = value(s->t, getNode(s->t, *inA), NULL);
        s->Top2V[top] = *inA;
        insertKey(s->O[in], top);

        /* Move pa forward */
        while(
          pa <= pb &&
	        incompatibleSize(s, pa) > 1 - dE
          //incompatibleSize(s, pa) > MAX_V - dE
        ) {
          pa = next(s, pa, 1+pb);
        }
      }
      inA++; /* Cycle step */
    }

    /* Add one outgoing */
    if(pa <= pb && -1 != *outA){
      if(NULL == getNode(s->t, *outA))
	      pushDA(s->DegEx[out], *outA);
      else { /* Vertex exists in DAG */
        int top = value(s->t, getNode(s->t, *outA), NULL);
        s->Top2V[top] = *outA;
        insertKey(s->O[out], top);

        /* Move pb backward */
        while(
          pa <= pb &&
          incompatibleSize(s, pb) > 1 - dE
          //incompatibleSize(s, pb) > MAX_V - dE
        ) {
          pb = prev(s, pb, pa); // count this as processed node
        }
      }
      outA++; /* Cycle step */
    }
  }

  /* Selecting the actual split point */
  s->p = pa;
  if(pa <= pb){ /* Otherwise check failed */
    int count = 0;
    int pmin = INT_MAX;
    int p = pa;
    while(p <= pb){
      if(incompatibleSize(s, p) < pmin){
        pmin = incompatibleSize(s, p);
        count = 0;
      }
      if(incompatibleSize(s, p) == pmin) count++;
      p = next(s, p, 1+pb);
    }

    if(1 < count){
      count = ARC4RAND_UNIF(count);
      count++;
    }

    p = pa;
    while(0 < count){
      if(incompatibleSize(s, p) == pmin)
	      count--;
      if(0 < count){
	      p = next(s, p, 1+pb);
      }
    }

    /* Choose uniformly in interval */
    // if (p < pb && p+1 < next(s, p, 1+pb)){
    //   p += ARC4RAND_UNIF(next(s, p, 1+pb) - p);
    // }

    s->p = p;
  }

  /* Reload graph Adj */
  int *T = getInorder(s->O[in], s->G->E[v][in], s->Top2V, right);
  dumpDA(s->DegEx[in], T);
  T = getInorder(s->O[out], s->G->E[v][out], s->Top2V, left);
  dumpDA(s->DegEx[out], T);

  /* printf("State after check\n"); */
  /* printS(s); */

#ifndef NDEBUG
  assertState(s);
#endif /* NDEBUG */

  s->checkDE = 1 - incompatibleSize(s, s->p);

  /* fprintf(stderr, "pa = %d ", pa); */
  /* fprintf(stderr, "p = %d ", s->p); */
  /* fprintf(stderr, "pb = %d ", pb); */
  /* fprintf(stderr, "dE = %d ", dE); */
  /* fprintf(stderr, "\n"); */

  // printf("--- checkDE = %i, p = %i\n", s->checkDE, s->p);
  return s->checkDE;
}


int
checkDelete(sa_state s,
      int deteriorationLimit,
      int dE//vを追加時に削除してもよい頂点数
      )
{
  assert(deteriorationLimit + 1 > dE && "Check with impossible dE");
  
  if(checkO(s->o, 0)){
    /*s->reorderCount = 1;*/ // this fix does not work
    s->cList[s->cLs] = -1;
    reorganize(s->o, s->cList);
    const int INF = 0x5fffffff;
    for(int i = 0; i < MAX_V; i++) pdELog[i] = INF;
  }

  /* Candidate */
  int v = s->cList[0];

  /* Load trees O[in] and O[out] */
  clearTree(s->O[in]);
  clearTree(s->O[out]);
  resetDA(s->DegEx[in]);
  resetDA(s->DegEx[out]);

  int *inA = s->G->E[v][in];
  int *outA = s->G->E[v][out];
  /* Positions over Topological order */
  int pa = 0;
  int pb = size(s->t);

  while(pa <= pb &&
        !(-1 == *inA && -1 == *outA)
	){
    /* Increment the heuristic */
    checkO(s->o, 2);

    /* Add one incoming */
    if(pa <= pb && -1 != *inA){
      if(NULL == getNode(s->t, *inA))
	      pushDA(s->DegEx[in], *inA);
      else { /* Vertex exists in DAG */
        int top = value(s->t, getNode(s->t, *inA), NULL);
        s->Top2V[top] = *inA;
        insertKey(s->O[in], top);
        /* Move pa forward */
        while(
          pa <= pb &&
	        incompatibleSize(s, pa) > 1 - dE // 元が1 - dE
          //incompatibleSize(s, pa) > dE // 元が1 - dE
          //incompatibleSize(s, pa) > 1 - dE
        ) {
          pa = next(s, pa, 1 + pb);     
        }
      }
      inA++; /* Cycle step */
    }

    /* Add one outgoing */
    if(pa <= pb && -1 != *outA){
      if(NULL == getNode(s->t, *outA))
	      pushDA(s->DegEx[out], *outA);
      else { /* Vertex exists in DAG */
        int top = value(s->t, getNode(s->t, *outA), NULL);
        s->Top2V[top] = *outA;
        insertKey(s->O[out], top);
        /* Move pb backward */
        while(
          pa <= pb &&
          incompatibleSize(s, pb) > 1 - dE // 元が1 - dE
          //incompatibleSize(s, pb) > deteriorationLimit - dE // 元が1 - dE
          //incompatibleSize(s, pb) > dE
        ) {
          pb = prev(s, pb, pa); // count this as processed node
        }
      }
      outA++; /* Cycle step */
    }
  }

  /* Selecting the actual split point */
  s->p = pa;
  if(pa <= pb){ /* Otherwise check failed */
    int count = 0;
    int pmin = INT_MAX;
    int p = pa;
    while(p <= pb){
      if(incompatibleSize(s, p) < pmin){
        pmin = incompatibleSize(s, p);
        count = 0;
      }
      if(incompatibleSize(s, p) == pmin) count++;
      p = next(s, p, 1+pb);
    }

    if(1 < count){
      count = ARC4RAND_UNIF(count);
      count++;
    }

    p = pa;
    while(0 < count){
      if(incompatibleSize(s, p) == pmin)
	      count--;
      if(0 < count){
	      p = next(s, p, 1+pb);
      }
    }
    s->p = p;
  }


  /* Reload graph Adj */
  int *T = getInorder(s->O[in], s->G->E[v][in], s->Top2V, right);
  dumpDA(s->DegEx[in], T);
  T = getInorder(s->O[out], s->G->E[v][out], s->Top2V, left);
  dumpDA(s->DegEx[out], T);

  s->checkDE = 1 - incompatibleSize(s, s->p);
  //s->checkDE = -incompatibleSize(s, s->p);
  return s->checkDE;
}


int
checkNoTopo(sa_state s,
      int dE
      )
{
  // printf("cList (size t = %i): ", size(s->t));
  // for (int i = 0; i < s->cLs; ++i) printf("%i ", s->cList[i]);
  // printf("\n");
  assert(2 > dE && "Check with impossible dE");

  /* Candidate */
  int v = s->cList[0];

  /* printf("State before check\n"); */
  /* printS(s); */
  /* printf(">>>> Candidate is : %d\n", v); */

  /* Load trees O[in] and O[out] */
  clearTree(s->O[in]);
  clearTree(s->O[out]);
  resetDA(s->DegEx[in]);
  resetDA(s->DegEx[out]);

  int *inA = s->G->E[v][in];
  int *outA = s->G->E[v][out];
  /* Positions over Topological order */
  int pa = 0;
  int pb = size(s->t);

  while(pa <= pb &&
        !(-1 == *inA && -1 == *outA)
	){
    /* Increment the heuristic */
    checkO(s->o, 2);

    /* Add one incoming */
    if(pa <= pb && -1 != *inA){
      if(NULL == getNode(s->t, *inA))
	      pushDA(s->DegEx[in], *inA);
      else { /* Vertex exists in DAG */
        int top = value(s->t, getNode(s->t, *inA), NULL);
        s->Top2V[top] = *inA;
        insertKey(s->O[in], top);

        /* Move pa forward */
        while(
          pa <= pb &&
	        incompatibleSize(s, pa) > 1 - dE
        ) {
          pa = next(s, pa, 1+pb);
        }
      }
      inA++; /* Cycle step */
    }

    /* Add one outgoing */
    if(pa <= pb && -1 != *outA){
      if(NULL == getNode(s->t, *outA))
	      pushDA(s->DegEx[out], *outA);
      else { /* Vertex exists in DAG */
        int top = value(s->t, getNode(s->t, *outA), NULL);
        s->Top2V[top] = *outA;
        insertKey(s->O[out], top);

        /* Move pb backward */
        while(
          pa <= pb &&
          incompatibleSize(s, pb) > 1 - dE
        ) {
          pb = prev(s, pb, pa); // count this as processed node
        }
      }
      outA++; /* Cycle step */
    }
  }

  /* Selecting the actual split point */
  s->p = pa;
  if(pa <= pb){ /* Otherwise check failed */
    int count = 0;
    int pmin = INT_MAX;
    int p = pa;
    while(p <= pb){
      if(incompatibleSize(s, p) < pmin){
        pmin = incompatibleSize(s, p);
        count = 0;
      }
      if(incompatibleSize(s, p) == pmin) count++;
      p = next(s, p, 1+pb);
    }

    if(1 < count){
      count = ARC4RAND_UNIF(count);
      count++;
    }

    p = pa;
    while(0 < count){
      if(incompatibleSize(s, p) == pmin)
	      count--;
      if(0 < count){
	      p = next(s, p, 1+pb);
      }
    }

    /* Choose uniformly in interval */
    // if (p < pb && p+1 < next(s, p, 1+pb)){
    //   p += ARC4RAND_UNIF(next(s, p, 1+pb) - p);
    // }

    s->p = p;
  }

  /* Reload graph Adj */
  int *T = getInorder(s->O[in], s->G->E[v][in], s->Top2V, right);
  dumpDA(s->DegEx[in], T);
  T = getInorder(s->O[out], s->G->E[v][out], s->Top2V, left);
  dumpDA(s->DegEx[out], T);

  /* printf("State after check\n"); */
  /* printS(s); */

#ifndef NDEBUG
  assertState(s);
#endif /* NDEBUG */

  s->checkDE = 1 - incompatibleSize(s, s->p);

  /* fprintf(stderr, "pa = %d ", pa); */
  /* fprintf(stderr, "p = %d ", s->p); */
  /* fprintf(stderr, "pb = %d ", pb); */
  /* fprintf(stderr, "dE = %d ", dE); */
  /* fprintf(stderr, "\n"); */

  // printf("--- checkDE = %i, p = %i\n", s->checkDE, s->p);
  return s->checkDE;
}

int
checkRepush(sa_state s,
      int v
      )
{
  int dE = 0;
  assert(2 > dE && "Check with impossible dE");
  neighborInitPdELog(s, v);

  /* Candidate */

  /* Load trees O[in] and O[out] */
  clearTree(s->O[in]);
  clearTree(s->O[out]);
  resetDA(s->DegEx[in]);
  resetDA(s->DegEx[out]);

  int *inA = s->G->E[v][in];
  int *outA = s->G->E[v][out];
  /* Positions over Topological order */
  int pa = 0;
  int pb = size(s->t);

  while(pa <= pb &&
        !(-1 == *inA && -1 == *outA)
	){
    /* Increment the heuristic */
    checkO(s->o, 2);
    /* Add one incoming */
    if(pa <= pb && -1 != *inA){
      if(NULL == getNode(s->t, *inA))
	      pushDA(s->DegEx[in], *inA);
      else { /* Vertex exists in DAG */
        int top = value(s->t, getNode(s->t, *inA), NULL);
        s->Top2V[top] = *inA;
        insertKey(s->O[in], top);

        /* Move pa forward */
        while(
          pa <= pb &&
	        incompatibleSize(s, pa) > 1 - dE
        ) {
          pa = next(s, pa, 1+pb);
        }
      }
      inA++; /* Cycle step */
    }

    /* Add one outgoing */
    if(pa <= pb && -1 != *outA){
      if(NULL == getNode(s->t, *outA))
	      pushDA(s->DegEx[out], *outA);
      else { /* Vertex exists in DAG */
        int top = value(s->t, getNode(s->t, *outA), NULL);
        s->Top2V[top] = *outA;
        insertKey(s->O[out], top);

        /* Move pb backward */
        while(
          pa <= pb &&
          incompatibleSize(s, pb) > 1 - dE
        ) {
          pb = prev(s, pb, pa); // count this as processed node
        }
      }
      outA++; /* Cycle step */
    }
  }

  /* Selecting the actual split point */
  s->p = pa;
  if(pa <= pb){ /* Otherwise check failed */
    int count = 0;
    int pmin = INT_MAX;
    int p = pa;
    while(p <= pb){
      if(incompatibleSize(s, p) < pmin){
        pmin = incompatibleSize(s, p);
        count = 0;
      }
      if(incompatibleSize(s, p) == pmin) count++;
      p = next(s, p, 1+pb);
    }

    if(1 < count){
      count = ARC4RAND_UNIF(count);
      count++;
    }

    p = pa;
    while(0 < count){
      if(incompatibleSize(s, p) == pmin)
	      count--;
      if(0 < count){
	      p = next(s, p, 1+pb);
      }
    }

    /* Choose uniformly in interval */
    // if (p < pb && p+1 < next(s, p, 1+pb)){
    //   p += ARC4RAND_UNIF(next(s, p, 1+pb) - p);
    // }

    s->p = p;
  }

  /* Reload graph Adj */
  int *T = getInorder(s->O[in], s->G->E[v][in], s->Top2V, right);
  dumpDA(s->DegEx[in], T);
  T = getInorder(s->O[out], s->G->E[v][out], s->Top2V, left);
  dumpDA(s->DegEx[out], T);

  s->checkDE = 1 - incompatibleSize(s, s->p);

  /* fprintf(stderr, "pa = %d ", pa); */
  /* fprintf(stderr, "p = %d ", s->p); */
  /* fprintf(stderr, "pb = %d ", pb); */
  /* fprintf(stderr, "dE = %d ", dE); */
  /* fprintf(stderr, "\n"); */

  // printf("--- checkDE = %i, p = %i\n", s->checkDE, s->p);
  return s->checkDE;
}

void
topoSort(sa_state s
        )
{
    s->cList[s->cLs] = -1;
    reorganize(s->o, s->cList);
}

int
checkB(sa_state s,
      int dE,
      int v  /* Candidate */
      )
{
  /* Apply re-organization heuristic */
  // if(checkO(s->o, 0)){
  //   /*s->reorderCount = 1;*/ // this fix does not work
  //   s->cList[s->cLs] = -1;
  //   reorganize(s->o, s->cList);
  // }
  /* Load trees O[in] and O[out] */
  clearTree(s->O[in]);
  clearTree(s->O[out]);
  resetDA(s->DegEx[in]);
  resetDA(s->DegEx[out]);
  int *inA = s->G->E[v][in];
  int *outA = s->G->E[v][out];
  /* Positions over Topological order */
  int pa = 0;
  int pb = size(s->t);
  while(pa <= pb &&
        !(-1 == *inA && -1 == *outA)
	){
    /* Increment the heuristic */
    //checkO(s->o, 2);
    /* Add one incoming */
    if(pa <= pb && -1 != *inA){
      if(NULL == getNode(s->t, *inA))
	      pushDA(s->DegEx[in], *inA);
      else { /* Vertex exists in DAG */
        int top = value(s->t, getNode(s->t, *inA), NULL);
        s->Top2V[top] = *inA;
        insertKey(s->O[in], top);
        /* Move pa forward */
        while(
          pa <= pb &&
	        incompatibleSize(s, pa) > 1 - dE
        ) {
          pa = next(s, pa, 1+pb);
        }
      }
      inA++; /* Cycle step */
    }
    /* Add one outgoing */
    if(pa <= pb && -1 != *outA){
      if(NULL == getNode(s->t, *outA))
	      pushDA(s->DegEx[out], *outA);
      else { /* Vertex exists in DAG */
        int top = value(s->t, getNode(s->t, *outA), NULL);
        s->Top2V[top] = *outA;
        insertKey(s->O[out], top);
        /* Move pb backward */
        while(
          pa <= pb &&
          incompatibleSize(s, pb) > 1 - dE
        ) {
          pb = prev(s, pb, pa); // count this as processed node
        }
      }
      outA++; /* Cycle step */
    }
  }
  /* Selecting the actual split point */
  s->p = pa;
  if(pa <= pb){ /* Otherwise check failed */
    int count = 0;
    int pmin = INT_MAX;
    int p = pa;
    while(p <= pb){
      if(incompatibleSize(s, p) < pmin){
        pmin = incompatibleSize(s, p);
        count = 0;
      }
      if(incompatibleSize(s, p) == pmin) count++;
      p = next(s, p, 1+pb);
    }
    if(1 < count){
      count = ARC4RAND_UNIF(count);
      count++;
    }
    p = pa;
    while(0 < count){
      if(incompatibleSize(s, p) == pmin)
	      count--;
      if(0 < count){
	      p = next(s, p, 1+pb);
      }
    }
    /* Choose uniformly in interval */
    // if (p < pb && p+1 < next(s, p, 1+pb)){
    //   p += ARC4RAND_UNIF(next(s, p, 1+pb) - p);
    // }
    s->p = p;
  }
  /* Reload graph Adj */
  int *T = getInorder(s->O[in], s->G->E[v][in], s->Top2V, right);
  dumpDA(s->DegEx[in], T);
  T = getInorder(s->O[out], s->G->E[v][out], s->Top2V, left);
  dumpDA(s->DegEx[out], T);
  s->checkDE = 1 - incompatibleSize(s, s->p);
  return s->checkDE;
}

/* Executes validated transition */
/* Returns number of new candidates */
void
execute(sa_state s
        )
{
  int v = s->cList[0]; /* The candidate */

  s->p = pLog[v];
  neighborInitPdELog(s, v);

  cList2[cLs2] = v;
  iterCList2[v] = cLs2;
  cLs2++;

  s->cLs--;
  int u = s->cList[s->cLs];
  int iterV = iterCList[v];
  iterCList[v] = iterCList[u];
  iterCList[u] = iterV;
  swapI(&(s->cList[0]), &(s->cList[s->cLs]));
  /* printf("Adopting candidate %d\n", v); */
#ifdef POSTSCRIPT
  fprintf(LOG, "%d %d i\n", 1+(v%POST_NB_VERT), 1+(v/POST_NB_VERT));
#endif /* POSTSCRIPT */

  /* 1 - add v */
  insertInorderKey(s->t, s->p, v);

  assert(size(s->t) > 0);

  /* 2 - remove incompatible nodes */
  int *condemned = s->condemned; /* Relabel and re-use */
  splitSt(s->O[in], s->p, right);
  condemned = getInorder(s->O[in], condemned, s->Top2V, right);

  splitSt(s->O[out], s->p-1, left);
  condemned = getInorder(s->O[out], condemned, s->Top2V, left);
  *condemned = -1;

  condemned = s->condemned;
  while(-1 != *condemned && 0 >= s->checkDE){
    int u = *condemned;
    removeN(s->t, getNode(s->t, u));
    neighborInitPdELog(s, u);
    /* Return back to candidate list, as current candidate. */
#ifdef POSTSCRIPT
    fprintf(LOG, "%d %d c\n",
	    1+(*condemned%POST_NB_VERT),
	    1+(*condemned/POST_NB_VERT));
#endif /* POSTSCRIPT */

    cLs2--;
    int iterU = iterCList2[u];
    cList2[iterU] = cList2[cLs2];
    iterCList2[cList2[iterU]] = iterU;

    iterCList[u] = s->cLs;
    s->cList[s->cLs] = u;
    s->cLs++;
    condemned++;
    s->checkDE++;
  }
}

void
executeNature(sa_state s
        )
{
  int v = s->cList[0]; /* The candidate */
  s->cLs--;
  swapI(&(s->cList[0]), &(s->cList[s->cLs]));
  /* printf("Adopting candidate %d\n", v); */
#ifdef POSTSCRIPT
  fprintf(LOG, "%d %d i\n", 1+(v%POST_NB_VERT), 1+(v/POST_NB_VERT));
#endif /* POSTSCRIPT */

  /* 1 - add v */
  insertInorderKey(s->t, s->p, v);

  assert(size(s->t) > 0);

  /* 2 - remove incompatible nodes */
  int *condemned = s->condemned; /* Relabel and re-use */
  splitSt(s->O[in], s->p, right);
  condemned = getInorder(s->O[in], condemned, s->Top2V, right);

  splitSt(s->O[out], s->p-1, left);
  condemned = getInorder(s->O[out], condemned, s->Top2V, left);
  *condemned = -1;

  condemned = s->condemned;
  while(-1 != *condemned && 0 >= s->checkDE){
    int u = *condemned;
    removeN(s->t, getNode(s->t, u));
    neighborInitPdELog(s, u);
    /* Return back to candidate list, as current candidate. */
#ifdef POSTSCRIPT
    fprintf(LOG, "%d %d c\n",
	    1+(*condemned%POST_NB_VERT),
	    1+(*condemned/POST_NB_VERT));
#endif /* POSTSCRIPT */

    s->cList[s->cLs] = u;
    s->cLs++;
    condemned++;
    s->checkDE++;
  }
}


void
insertV(sa_state s,
      int v
      )
{
  insertInorderKey(s->t, s->p, v);
  splitSt(s->O[in], s->p, right);
  splitSt(s->O[out], s->p-1, left);
}

void 
insertPV(sa_state s,
      int v,
      int p
      )
{
  neighborInitPdELog(s, v);
  insertInorderKey(s->t, p, v);
  splitSt(s->O[in], p, right);
  splitSt(s->O[out], p-1, left);
}

int
getSP(sa_state s
      )
{
  return s->p;
}

int
swapV(sa_state s,
      int v
      )
{
  insertV(s, v);
  int u;
  int *condemned = s->condemned; /* Relabel and re-use */
  splitSt(s->O[in], s->p, right);
  condemned = getInorder(s->O[in], condemned, s->Top2V, right);
  splitSt(s->O[out], s->p-1, left);
  condemned = getInorder(s->O[out], condemned, s->Top2V, left);
  *condemned = -1;
  condemned = s->condemned;
  while(-1 != *condemned)
  {
    u = *condemned;
    //swapVers[*condemned] = true;
    removeN(s->t, getNode(s->t, *condemned));
    condemned++;
  }
  if(u < 0) {
    fprintf(stderr, "Warning: Swap non existing vertex %d\n", v);
    return -1;
  }
  return u;
}

int
swapPV(sa_state s,
      int v,
      int p
      )
{
  insertPV(s, v, p);
  int u;
  int *condemned = s->condemned; /* Relabel and re-use */
  splitSt(s->O[in], p, right);
  condemned = getInorder(s->O[in], condemned, s->Top2V, right);
  splitSt(s->O[out], p-1, left);
  condemned = getInorder(s->O[out], condemned, s->Top2V, left);
  *condemned = -1;
  condemned = s->condemned;
  const int INF = 0x5fffffff;
  while(-1 != *condemned)
  {
    u = *condemned;
    for(int i=0;i<s->G->outDeg[u];i++){
      int to = s->G->E[u][out][i];
      pdELog[to] = INF;
      pLog[to] = -1;
    }
    for(int i=0;i<s->G->inDeg[u];i++){
      int to = s->G->E[u][in][i];
      pdELog[to] = INF;
      pLog[to] = -1;
    }
    removeN(s->t, getNode(s->t, u));
    condemned++;
  }
  if(u < 0) {
    fprintf(stderr, "Warning: Swap non existing vertex %d\n", v);
    return -1;
  }
  return u;
}

int
eraseV(sa_state s,
      int v
      )
{
  nodes n = getNode(s->t, v);
  if(NULL != n){
    removeN(s->t, n);
    return 0;
  }
  else
  {
    fprintf(stderr, "Warning: Erase non existing vertex %d\n", v);
    return -1;
  }
}

int *
getSketch(sa_state s,
	        int *L
          )
{
  int *R = getInorder(s->t, L, NULL, left);

  *R = -1;
  return R;
}
