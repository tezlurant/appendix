#ifndef _SA_STATE_H
#define _SA_STATE_H

#include "graph.h"
const int MAX_V = 5000;
int pdELog[MAX_V]; //初期値は2で未計算を示す
int pLog[MAX_V]; //pdELogとセット
int iterCList[MAX_V]; // 解候補に含まれているかどうか(含まれてるなら何番目か)
int cList2[MAX_V];
int iterCList2[MAX_V]; // 解集合に含まれているかどうか(含まれてるなら何番目か)
int repushL = 0; // repushに入ってるサイズ
int repushList[MAX_V]; //repushで隣接点を評価するためのリスト
bool stagnantBitList[MAX_V];
bool currentStagnantBitList[MAX_V];
int cLs2 = 0;
#ifdef POSTSCRIPT
extern FILE *LOG;
extern int POST_NB_VERT;
#endif /* POSTSCRIPT */

typedef struct sa_state_ *sa_state;

sa_state
allocS(graph G
       );

sa_state
dupS(sa_state s
     );

sa_state
cpyS(sa_state dst,
     sa_state src
    );

void
freeS(sa_state s
      );

/* Set sa_state vertex to s */
void
prepareS(sa_state s,
	 int* A, /* Array with vertexes */
	 int l	 /* length of the array */
	 );

void
prepareS2(sa_state s,
          int* A, /* Array with vertexes */
          int l,	 /* length of the array */
          int B,
          int lB
          );

void
prepareSwithKnowSolution(sa_state s,
         int*  A, /* Array with vertexes */
         int   l, /* length of the array */
         int*  C, /* Array all vertexes in SCC */
         int   sC  /* size of the SCC */
        );

int 
swapSelect(sa_state s,
           int p
           );

int
selectV(sa_state s
        );

int
swapSelectSE(sa_state s,
       int p,
       int q
       );

/* Chooses a potential transition */
int
choose(sa_state s
       );

/* Chooses a potential transition 無駄を省く*/
int
chooseSB(sa_state s,
        int selectBorder
       );

/* Chooses a potential transition log(dv)のサイズから優先度選択.ver */
int
chooseLogdV(sa_state s
      );

/* Chooses a potential transition Assing.ver */
int
chooseAssign(sa_state s,
             int swapIdx 
             );

/* Returns the size of the candidate list */
int
cLsSize(sa_state s
        );

void
cLsZero(sa_state s
        );

int
cLsAdd(sa_state s,
        int v
        );
int
cLsGet(sa_state s,
        int i
        );

void
updateStagnantBitList(sa_state s
                      );

int
ansDistance(sa_state s
            );

void
neighborInitPdELog(sa_state s,
                   int v
                   );

void mergeSort2(int left, 
                int right
                );

void
initSolution(sa_state s
                );

/* Returns the proposed energy variation */
int
check(sa_state s,
      int dE
      );

int
checkDelete(sa_state s,
            int deteriorationLimit,
            int E
            );

int
checkNoTopo(sa_state s,
            int dE
            );

int
checkRepush(sa_state s,
            int v
            );

void 
topoSort(sa_state s
         );

int
checkB(sa_state s,
       int dE,
       int v
       );

/* Executes validated transition */
/* Returns number of new candidates */
void
execute(sa_state s
        );

void
executeNature(sa_state s
              );

void
insertV(sa_state s,
        int v
        );

void
insertPV(sa_state s,
        int v,
        int p
        );

//現在のs->pを返す
int
getSP(sa_state s
      );

int
swapV(sa_state s,
       int v
       );

int
swapPV(sa_state s,
        int v,
        int p
        );

// ノードを削除する関数, vがノード
int
eraseV(sa_state s,
        int v
        );

int
greedy(sa_state s
       );

/* Return sa_state energy */
int
getE(sa_state s
     );

void
printS(sa_state s
       );

/* Load sa_state into list L */
int *
getSketch(sa_state s,
	  int *L
          );

#endif /* _SA_STATE_H */
