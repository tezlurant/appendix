#include <stdlib.h>
#include <bsd/stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fenv.h>
#include <assert.h>
#include <limits.h>
#include <time.h>
#include <set>
#include <random> // 乱数ライブラリ

#include "sa_state.h"
#include "sa_lib.h"
#include "aux.h"

#define LIMIT_TIME
#define MAXTIME 1000
#define BATCH_SIZE 0
extern unsigned long nbProcVert;
extern unsigned long nbDfsCalls;

extern __thread SA_s internal_lib_sa_;

/* Used for defining calibration values */
double
getTemperature(double p, /* probability */
	       int dE /* Energy delta */
	       )
{
  assert(0.5 > p && "Invalid probability to define temp.");
  double res = dE / log2((1/p)-1);
  // printf("temp: %f\n", res);
  return res;
}

int *
executeSA5(
    sa_state s, /* State struct to use */
	  unsigned long int n, /* Number of iterations to run */
	  int *L, /* Buffer for storing solution */
	  unsigned long int repOps /* Report batch size */
	  )
{
  //int cnt = 0; /* Counter */
  const int INF = 0x5fffffff;
  SA_s sa = internal_lib_sa_;
  int *R = L; /* Return value */
  for(int i=0;i<MAX_V;i++){
    pdELog[i] = INF;
  }
  long double T0 = 1.0;
  long double T = T0;
  long double alpha = 0.99; //0.9とかもあり
  int sgv = SA_get_nbVerts(); //グラフ頂点数
  int pdE; //pの位置に挿入するために削除する必要のある要素数
  fesetround(FE_TONEAREST); /* For dE computation */
  //printf("%lu\n",n);
  // function re-enters multiple times (compute temperature) 
  //char profFileName[1<<10];
  //sprintf(profFileName, "SA_n=%lu_%d週目.txt", n, cnt);
  //FILE *fp = fopen(profFileName, "w");
  n = INF;
  int *typeVerts = sa->typeVerts; // 0: 削除頂点(追加候補), 1: 非巡回頂点
  int *randomVerts = sa->randomVerts; // ランダムに頂点を選択するための配列
  //ここでCList2関係の配列を書き換える
  cLs2 = 0;
  for(int i = 0; i < sgv; i++){
    cList2[i] = -1;
  }
  for(int i=0;i<cLsSize(s);i++){
    int v = cLsGet(s, i);
    cList2[v] = -2;
  }
  for(int i=0;i<sgv;i++){
    if(cList2[i] == -1){
      cList2[cLs2] = i;
      iterCList2[i] = cLs2;
      cLs2++;
    }
  }
  updateStagnantBitList(s);

  int stagnantCount = 0; // 停滞回数

  const long double eps = 1e-9; // 誤差
  const int initDeteriorationLimit = 1; // 改悪の初期値
  const int deteriorationLimitMax = 50; // 改悪の上限の最大値
  int deteriorationLimit = initDeteriorationLimit; // 改悪の上限(defaultは1,maxは5?)
  //std::mt19937 engine(0); // シードに現在時刻を使用
  //std::uniform_real_distribution<double> dist(0.0, 1.0); // 0.0以上1.0未満の一様乱数
  for (unsigned long int i = 0; i < n; i++, sa->ops++) {
    int dE = 0; // 削除してもいい要素数 (dE + 1を許容(交換まで))
    if(0 != repOps && 0 == (sa->ops % repOps) || sa->ops < INIT_OPS) {
      SA_printLine(sa->currSCC);
    }
    long double B = (long double)rand()/RAND_MAX; // for testing
    B+=eps;
    float t; /* temporary variable */
    t = logf(1/(double)B);
    t *= T;
    dE = lrintf(t);
    if(dE > deteriorationLimit) dE = deteriorationLimit;
    dE = -dE;
    if (choose(s)) {
      break; // invalid sa_state
    }
    int sp;
    int v = cLsGet(s, 0);
    pdE = check(s, dE);
    
    //SA_printLine2(T);
    if (
       dE <= pdE
       ) 
    {
      if(getE(s) > SA_maxDAG()){
        SA_updateMaxDAG(getE(s));
        deteriorationLimit = initDeteriorationLimit;
      }
      if (getE(s) > sa->maxE[sa->currSCC]) {
        sa->maxE[sa->currSCC] = getE(s);
        if (0 > pdE) {
          if(0 != repOps && 0 == (sa->ops % repOps) || sa->ops < INIT_OPS) {
            SA_printLine(sa->currSCC);
          }
          R = getSketch(s, L);
        }
      }
      int prevE = getE(s);
      executeNature(s);
      //execute(s);
      int currentE = getE(s);
      if(currentE > prevE){
        stagnantCount = 0;
      }
    }
    stagnantCount++;
    long long currentF =  SA_get_nbVerts()-getE(s) + SA_get_nbLoop();
    long double fp = pow(2,(long double)sgv/(long double)currentF);
    //long double fp = (long double)sgv/(long double)currentF;
    long long sgvFp = 2*sgv*fp;
    if(stagnantCount >= sgvFp){ //5,10,20
      //解のサイズで判定(3%で判定)
      // if(SA_maxDAG()-getE(s)>(3*sgv)/100){
      //   deteriorationLimit = initDeteriorationLimit - 1;
      // }
      //解の距離で判定 (何%?離れていればいい？)
      if(ansDistance(s) > sgv/5){
        deteriorationLimit = initDeteriorationLimit - 1;
        updateStagnantBitList(s);
      }
      //Tに何か加算or乗算
      // 改悪の上限を1加算する(ただし上限は定める1000なら5まで？)
      deteriorationLimit++;
      if(deteriorationLimit > deteriorationLimitMax){
        deteriorationLimit = deteriorationLimitMax;
      }
      // 解が最適解の更新が行われた場合は上限をリセット
      
      // Tをリセット
      TIMER_T curTime;
		  TIMER_READ(curTime);
      double currentTime = CLOCK_DIFF_MS(sa->startTime, curTime); 
      //T = T0*(1-currentTime/(double)MAXTIME); //deterLimによって倍率変更する？
      double beta = 0.5;
      T = T0*(1 - beta/((double)deteriorationLimit));
      //T = T0;
      //T = currentF*0.001*initDeteriorationLimit;
      for(int i = 0; i < sgv; i++){
        pdELog[i] = INF;
      }
      stagnantCount = 0;
    }
    
    /* End update */
    T*=alpha;

#ifdef LIMIT_TIME
		TIMER_T curTime;
		TIMER_READ(curTime);
		if (CLOCK_DIFF_MS(sa->startTime, curTime) > MAXTIME)
			goto EXITPOS_SA;
#endif
  }

EXITPOS_SA:
  // test1. 1回の実行でSA or greedyで局所最適解を求める
  //printf("dE = 1: %d, dE = 0: %d\n", dECount[1], dECount[0]);
  // ここまで
  //fclose(fp);
  if(getE(s) > SA_maxDAG()){
    SA_updateMaxDAG(getE(s));
  }
  if (getE(s) >= sa->maxE[sa->currSCC]) {
    sa->maxE[sa->currSCC] = getE(s);
    R = getSketch(s, L);
  }
  
  if (R == L) {
    R = getSketch(s, L);
  }
  assert(R != L && "No solution returned");
  assert(sa->currSCCsize > R - L && "Wrong sized solution");

  return R;
}