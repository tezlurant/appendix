#include "sa_lib.h"
#include "aux.h"
#include "simulatedAnnealing.h"
#include "simulatedAnnealing.cpp"
#include "sa_state.cpp"
#include "graph.cpp"
#include <math.h>
#include <assert.h>
#include <time.h>

//#define EXESTING_TECHNIQUE

#define LIMIT_TIME
#define MAXTIME 1000

__thread SA_s internal_lib_sa_;
__thread FILE* internal_prof_file_;

extern unsigned long nbProcVert;
extern unsigned long nbDfsCalls;
unsigned long notUsedBudget = 0;

// call this after the graph is known
static void
allocSAAlgStructs()
{
  SA_s sa = internal_lib_sa_;
  int const minNumberOfBatches = 10; // tentative
  int const minSizeOfBatches = 100; 

  sa->accScore = 0;
  sa->ops = 0;
  sa->G->maxDAG = 0;

  // int randval = arc4random();
  // srandom(randval); // for testing
  // printf(" >>> seed is %i <<< \n", randval);
  // srandom(-32777779); // for testing

  // printf("Found %i SSCs in graph\n", sa->order[0]); // TODO: check BB for this

#ifdef POSTSCRIPT
  /* Use this flag only to create animations for the 100x100 torus */
  POST_NB_VERT = ceil(sqrt(sa->G->v));
  if(NULL == LOG)
    LOG = fopen("log", "w");

  fprintf(LOG, "/init { /N exch def 340 N div dup scale -1 -1 translate \
             /Palatino-Roman 0.8 selectfont 0 setlinewidth \
             1 1 N { 1 1 N { 1 index c } for pop } for } bind def \
     /r { translate 0 0 1 1 4 copy rectfill 0 setgray rectstroke } bind def \
     /i { gsave 0.5 setgray r grestore } bind def \
     /d { gsave 0.0 setgray r grestore } bind def \
     /q { gsave r 0.5 0.28 translate (Q) dup stringwidth pop -2 div 0 moveto \
          1 setgray show grestore } bind def \
     /c { gsave 1 setgray r grestore } bind def\n");
  fprintf(LOG, "%d init\n", POST_NB_VERT);
  for (int d = sa->G->v; d < POST_NB_VERT*POST_NB_VERT; ++d) {
    fprintf(LOG, "%d %d d\n",
      1+(d%POST_NB_VERT),
      1+(d/POST_NB_VERT)
    );
  }
#endif /* POSTSCRIPT */

  /* printG(rG); */
  int numberOfBatches = minNumberOfBatches;
  int sizeForTheBatch = 1 + (sa->n / numberOfBatches);

  while (sizeForTheBatch > 10 * minSizeOfBatches) {
    numberOfBatches++;
    sizeForTheBatch = 1 + (sa->n / numberOfBatches);
  }


  for (int i = 0; i < sa->order[0]; ++i) {
    if (sa->SCCsize[i] < 3) {
      notUsedBudget += 1+((sa->n * sa->SCCsize[i]) / sa->rG->v);
    }
  }
  // sa->n += notUsedBudget;
  int *A = sa->order + 1;
  // int B[1000];
  // //int *A = sa->order;
  // int lengB = 0;
  // //ここAの中身を更新
  // int n = sa->G->v;
  // int cntCandGreedyV = n;
  // for(int i = 0; i < n; i++)
  // {
  //   forestV[i] = false;
  //   candGreedyV[i] = i;
  // }
  // bool lpW = false;
  // auto dfs = [&](auto &&self, int v, int oriV) ->void
  // {
  //   if(lpW) return;
  //   for(int i=0;i<sa->G->outDeg[v];i++)
  //   {
  //     int nv = sa->G->E[v][out][i];
  //     if(nv == oriV)
  //     {
  //       lpW = true;
  //       return;
  //     }
  //     if(forestV[nv])
  //     {
  //       self(self, nv, oriV);
  //     }
  //   }
  //   return;
  // };

  // for(int i = 0; i < n; i++){
  //   int r = rand() % cntCandGreedyV;
  //   int v = candGreedyV[r];
  //   candGreedyV[r] = candGreedyV[0];
  //   candGreedyV[0] = v;
  //   lpW = false;
  //   for(int j=0;j<sa->G->outDeg[v];j++){
  //     int nv = sa->G->E[v][out][j];
  //     if(forestV[nv]){
  //       dfs(dfs, nv, v);
  //     }
  //     if(lpW) break;
  //   }
  //   if(!lpW){
  //     forestV[v] = true;
  //     B[lengB] = v;
  //     lengB++;
  //   }
  //   cntCandGreedyV--;
  //   candGreedyV[0] = candGreedyV[cntCandGreedyV];
  //   candGreedyV[cntCandGreedyV] = v;
  // }

  //ここまで
  //printf("n=%d, lengA=%d\n", n, lengB);
  int k = 0;
  for (int i = 0; i < sa->order[0]; ++i) {
      sa->G->maxDAG = getE(sa->s[0]);
    //printf("maxDAG1=%d\n", sa->G->maxDAG);
    prepareS(sa->s[i], A, sa->SCCsize[i]);
    //prepareS2(sa->s[i], A, sa->SCCsize[i], B, lengB);
    //sa->G->maxDAG = getE(sa->s[0]);
    //printf("maxDAG2=%d\n", sa->G->maxDAG);
#if BATCH_SIZE > 0
    sa->bestState[i] = cpyS(sa->bestState[i], sa->s[i]);
#endif
    A += sa->SCCsize[i];
#if BATCH_SIZE > 0
    if (sa->bestState[i])
      freeS(sa->bestState[i]);
    sa->bestState[i] = allocS(sa->rG);
#endif
    for (int j = 0; j < sa->SCCsize[i]; ++j)
      *(sa->vertsPerSCC[i]+j) = k++;
    *(sa->vertsPerSCC[i] + sa->SCCsize[i]) = -1;
    // sa->budgetPerSCC[i] = ((sa->n * sa->SCCsize[i]) / sa->rG->v);
    sa->budgetPerSCC[i] = (((sa->n+notUsedBudget) * sa->SCCsize[i]) / sa->rG->v) >> BATCH_FACTOR; // allows for a first run over all SCCs
    sa->budgetPerSCC[i] /= numberOfBatches;
    if (0 >= sa->budgetPerSCC[i]) sa->budgetPerSCC[i] = 2;
    // sa->budgetPerSCC[i] += sa->rG->v; // add a little extra to compensate small SCCs
    sa->initBudgetPerSCC[i] = sa->budgetPerSCC[i];
    sa->accSolPerSCC[i] = 0;
  }

  sa->numberOfBatch = numberOfBatches;
  // sa->hotArray = (double*)malloc(sizeof(double)*sa->numberOfBatch);
  // sa->coldArray = (double*)malloc(sizeof(double)*sa->numberOfBatch);

  for(int i = 0; i < sa->G->v; i++){
    sa->typeVerts[i] = 1;
    sa->randomVerts[i] = i;
  }
}

int
SA_get_nbVerts()
{
  SA_s sa = internal_lib_sa_;
  return sa->G->v;
}

int
SA_get_nbEdges()
{
  SA_s sa = internal_lib_sa_;
  return sa->G->e;
}

int
SA_get_nbLoop()
{
  SA_s sa = internal_lib_sa_;
  return sa->G->s1_count;
}

int
SA_reset()
{
  SA_s sa = internal_lib_sa_;
  allocSAAlgStructs();
  clock_gettime(CLOCK_MONOTONIC_RAW, &sa->startTime);
  return 0;
}

int
SA_init_G(SA_parameters_s params, graph G)
{
  SA_s sa = (SA_s)malloc(sizeof(struct SA_));
  internal_lib_sa_ = sa;
  internal_prof_file_ = stdout;

  sa->repOps = params.repOps;
  sa->n = params.n;
  sa->hot = params.hot;
  sa->hotD = params.hotD;
  sa->cold = params.cold;
  sa->coldD = params.coldD;

  G->t = NULL;
  G->d = NULL;
  sa->G = G;
  sa->order = (int*)malloc((1+sa->G->v)*sizeof(int));
  sa->SCCsize = (int*)malloc((1+sa->G->v)*sizeof(int));
  sa->SCCnbEdges = (int*)malloc((1+sa->G->v)*sizeof(int));
  sa->rG = sccRestrict(sa->G, sa->order, sa->SCCsize, sa->SCCnbEdges);

  sa->maxE = (int*)calloc(sa->order[0], sizeof(int));
  sa->s = (sa_state*)calloc(sa->order[0], sizeof(sa_state));
#if BATCH_SIZE > 0
  sa->bestState = (sa_state*)calloc(sa->G->v, sizeof(sa_state));
#endif

  sa->budgetPerSCC = (long*)malloc(sa->order[0]*sizeof(long));
  sa->initBudgetPerSCC = (long*)malloc(sa->order[0]*sizeof(long));
  sa->accSolPerSCC = (int*)malloc((sa->order[0]+1)*sizeof(int));
  sa->vertsPerSCC = (int**)calloc(sa->order[0]+1, sizeof(int*));
  sa->vertsPerSCC[0] = (int*)malloc((sa->G->v*2)*sizeof(int));
  sa->typeVerts = (int*)malloc(Node*sizeof(int));
  sa->randomVerts = (int*)malloc(Node*sizeof(int));

  for (int i = 0; i < sa->order[0]; ++i) {
    sa->s[i] = allocS(sa->rG);
#if BATCH_SIZE > 0
    sa->bestState[i] = allocS(sa->rG);
#endif
    if (i > 0)
      sa->vertsPerSCC[i] = sa->vertsPerSCC[i-1] + sa->SCCsize[i-1] + 1;
  }
  SA_reset();
  return 0;
}

int
SA_init_F(SA_parameters_s params, const char *filename)
{
  FILE *stream = fopen(filename, "r");
  graph G = loadG(stream);
  fclose(stream);
  return SA_init_G(params, G);
}

struct timespec
SA_getStartTime()
{
  SA_s sa = internal_lib_sa_;
  return sa->startTime;
}

void
SA_destroy()
{
  SA_s sa = internal_lib_sa_;
  for (int i = 0; i < sa->order[0]; ++i) {
    freeS(sa->s[i]);
  }
  free(sa->maxE);
  free(sa->vertsPerSCC[0]);
  free(sa->vertsPerSCC);
  free(sa->s);
  freeG(sa->G);
  free(sa->order);
  free(sa->SCCsize);
  free(sa->budgetPerSCC);
  free(sa->accSolPerSCC);
  free(sa->typeVerts);
  free(sa->randomVerts);
  free(sa);
}

int*  // loop with while(-1 != *bestSolution)
SA_getBestSolution()
{
  SA_s sa = internal_lib_sa_;
  return getBestSolutionBuffer(sa->G);
}

int
SA_maxDAG()
{
  SA_s sa = internal_lib_sa_;
  return sa->G->maxDAG;
}

void
SA_updateMaxDAG(int maxDAG)
{
  SA_s sa = internal_lib_sa_;
  sa->G->maxDAG = maxDAG;
}

void
SA_run()
{

  SA_printHeader();
  
  SA_s sa = internal_lib_sa_;
  /* With 10% prob. a -5 change */
  //double hotStart = getTemperature(sa->hot, sa->hotD);
  /* With 10% prob. a -1 change */
  //double coldStart = getTemperature(sa->cold, sa->coldD);

  // add small epsilon
  //double dT = (hotStart - coldStart - 0.00001)/((double)sa->numberOfBatch * (double)(1<<BATCH_FACTOR)); /* Update per batch */

  //double hot = hotStart;
  //double cold = hotStart - dT;

  int maxE = 0;

  // printf("hotStart=%f, coldStart=%f, dT=%f, hot=%f, cold=%f\n", hotStart, coldStart, dT, hot, cold);
  int *L = getBestSolutionBuffer(sa->G); /* Current best config */
  int *P;
  
  // clock_gettime(CLOCK_MONOTONIC_RAW, &sa->startTime);
  /* Running simulated Annealing */
  int *A = &(sa->order[1+sa->rG->v]);
  volatile struct timespec curTime, startTime = SA_getStartTime();

  sa->ops = 0;
  int noMoreBudget = 0;
  int isFirstRun = 1;

  SA_printLine(0);
  int countLoop = 0;
  sa->G->maxDAG = getE(sa->s[0]);
  #ifndef EXESTING_TECHNIQUE
  const int INF = 0x5fffffff;
  for(int i=0;i<MAX_V;i++){
    pdELog[i] = INF;
  }
  #endif

  while (
#ifdef LIMIT_TIME
  1
#else
  noMoreBudget < sa->order[0]
#endif
  ) {
    countLoop++;
    noMoreBudget = 0;
    for (int i = sa->order[0] - 1; 0 <= i; i--) {
      
      int *vertsPerSCC = sa->vertsPerSCC[i];
      long curOps = sa->ops;

      sa->currSCC = i;
      sa->currSCCsize = sa->SCCsize[i];

      if (0 >= sa->budgetPerSCC[i]) {
        noMoreBudget++;
        continue;
      }

      sa->accScore -= sa->accSolPerSCC[i];

      A -= sa->SCCsize[i];
      if (isFirstRun) {
        *vertsPerSCC = *A;
      }
      SA_printLine(i);
      if(sa->SCCnbEdges[i] == sa->SCCsize[i]*sa->SCCsize[i]-sa->SCCsize[i]) {
        vertsPerSCC++;
        *vertsPerSCC = -1;
        sa->budgetPerSCC[i] = 0;
        noMoreBudget++;
        sa->maxE[i] = 1;
        
        // always print initial solution (approx. alg)
        SA_printLine(i);
        
        continue;
      }
      initSolution(sa->s[i]);    
      vertsPerSCC = executeSA5(sa->s[i],
        sa->budgetPerSCC[i],
        vertsPerSCC, sa->repOps);

      assert(-1 == *vertsPerSCC && "Wrong list ending");
      assert(sa->currSCCsize > vertsPerSCC - sa->vertsPerSCC[i] && "Too many vertexes in list");

      sa->accSolPerSCC[i] = 0;
      vertsPerSCC = sa->vertsPerSCC[i];
      while(-1 != *vertsPerSCC && vertsPerSCC - sa->vertsPerSCC[i] < sa->currSCCsize) {
        vertsPerSCC++;
        sa->accSolPerSCC[i]++;
      }
      sa->accScore += sa->accSolPerSCC[i];

      if(0 != sa->repOps && 0 == ((++sa->ops) % sa->repOps) || sa->ops < INIT_OPS) {
        SA_printLine(i);
      }
    } // end loop per SCC
    SA_printLine(0);
    A += sa->rG->v;
    // if (cold > coldStart) {
    //   hot = cold;
    //   cold -= dT;
    // }
    isFirstRun = 0;   
    if (noMoreBudget == sa->order[0]) {
      break;
    }
    
#ifdef LIMIT_TIME
		TIMER_T curTime;
		TIMER_READ(curTime);
    //printf("%f\n",CLOCK_DIFF_MS(sa->startTime, curTime));
		if (CLOCK_DIFF_MS(sa->startTime, curTime) > MAXTIME)
      {
        //printf("TIMEOUT\n");
			  goto EXITPOS_SA;
      }
#endif
  } // end loop budget

EXITPOS_SA:
  P = L;
  for (int i = 0; i < sa->order[0]; i++) {
    int *vertsPerSCC = sa->vertsPerSCC[i];
    while (*vertsPerSCC != -1) {
      assert(sa->G->v > L - P);
      *L = *vertsPerSCC;
      assert(-1 != *L);
      L++;
      vertsPerSCC++;
    }
  }
  
  *L = -1;
  L = P;
  //printf("countLoop=%d\n", countLoop);
  //printf("ops=%lu\n", sa->ops);
}

void
SA_set_prof_file(FILE *fp)
{
  internal_prof_file_ = fp;
}

void SA_printHeader()
{
  return;
  fprintf(internal_prof_file_, "# iter\ttime_ms\tscore\tmaxScr\taccScr\tSCCid\tSCCsize\n"); // \tvertPerIt\tcurTemp\tnbDfs
}

void SA_printLine(int scc_id)
{
  return;
  SA_s sa = internal_lib_sa_;

  struct timespec curTime;
  clock_gettime(CLOCK_MONOTONIC_RAW, (struct timespec*)&curTime);
  int maxE = 0;
  for (int j = 0; j < sa->order[0]; ++j) maxE += sa->maxE[j];
  fprintf(internal_prof_file_, "%ld\t%.3f\t%d\t%d\t%d\t%d\t%d\n",
    sa->ops,
    CLOCK_DIFF_MS(sa->startTime, curTime),
    getE(sa->s[scc_id]),
    maxE,
    sa->accScore,
    scc_id,
    sa->SCCsize[scc_id]
  );
}

void SA_printLine2(long double T)
{
  SA_s sa = internal_lib_sa_;
  
  struct timespec curTime;
  clock_gettime(CLOCK_MONOTONIC_RAW, (struct timespec*)&curTime);
  int maxE = SA_maxDAG();
  int currentE = getE(sa->s[0]);
  int tt = CLOCK_DIFF_MS(sa->startTime, curTime);
  fprintf(internal_prof_file_, "%d %d %d %Lf\n",
    tt,
    SA_get_nbVerts() - maxE,
    SA_get_nbVerts() - currentE,
    T
  );

  // struct timespec curTime;
  // clock_gettime(CLOCK_MONOTONIC_RAW, (struct timespec*)&curTime);
  // int maxE = SA_maxDAG();
  //     int tt = CLOCK_DIFF_MS(sa->startTime, curTime);
  // fprintf(internal_prof_file_, "%d %d\n",
  //   tt,
  //   SA_get_nbVerts() - maxE
  // );
}
