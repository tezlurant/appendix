#include<iostream>
#include "graph.h"

using namespace std;

int init_vCount, now_vCount;//頂点数
int init_eCount, now_eCount;//辺数
int loopCount = 0, validVertexCount = 0; //ループ数、有効頂点数
int inDeg[MAX_V], outDeg[MAX_V];//入次数、出次数
int loop[MAX_V], validVertex[MAX_V]; //ループ、有効頂点
int Edge[MAX_V][MAX_V]; //辺の集合

void updateDegree()//現在のグラフの操作性を更新
{
	for (int i = 0; i < init_vCount; i++)
	{
		inDeg[i] = 0;
		outDeg[i] = 0;
	}
	for (int i = 0; i < init_vCount; i++)
		for (int j = 0; j < init_vCount; j++)
		{
			if (Edge[i][j])
			{
				inDeg[j]++;
				outDeg[i]++;
			}
		}
}

void DeleteVertex(int i)//iは削除する頂点
{
	int x;
	for (x = 0; x < init_vCount; x++)
	{
		Edge[x][i] = 0;
		Edge[i][x] = 0;
	}
	updateDegree();
}

void initReduction()
{
	init_vCount = 0;
	now_vCount = 0;
	init_eCount = 0;
	now_eCount = 0;
	loopCount = 0;
	validVertexCount = 0;
	for (int i = 0; i < MAX_V; i++)
	{
		inDeg[i] = 0;
		outDeg[i] = 0;
		loop[i] = 0;
		validVertex[i] = 0;
		for (int j = 0; j < MAX_V; j++)
		{
			Edge[i][j] = 0;
		}
	}
}

void reduction()//5つの基本的な削減処理を使用してダイアグラムを縮小
{
	bool flag = true;
	while (flag)
	{
		flag = false;
		for (int v = 0; v < init_vCount; v++)
		{
			if(inDeg[v] > 0 || outDeg[v] > 0)//点が存在すると判断
			{
				if (Edge[v][v] == 1) {//自己ループ
					loop[loopCount] = v;
					DeleteVertex(v);
					loopCount++;
					flag = true;
				}				
				else if (inDeg[v] == 0) {//入度0
					DeleteVertex(v);
					flag = true;
				}
				else if (outDeg[v] == 0) {//出度0
					DeleteVertex(v);
					flag = true;
				}
				else if (inDeg[v] == 1) {//入度1
					int u;
					for (u = 0; u < init_vCount; u++) {
						if (Edge[u][v] == 1) {
							Edge[u][v] = 0;
							for (int y = 0; y < init_vCount; y++)
							if (Edge[v][y] == 1)
							{
								Edge[v][y] = 0;
								Edge[u][y] = 1;
							}
							DeleteVertex(v);
							updateDegree();
							break;
						}
					}
					flag = true;
				}
				else if (outDeg[v] == 1) {//出度1
					int u;
					for (u = 0; u < init_vCount; u++) {
						if (Edge[v][u] == 1) {
							Edge[v][u] = 0;
							for (int x = 0; x < init_vCount; x++)
							if (Edge[x][v] == 1)
							{
								Edge[x][v] = 0;
								Edge[x][u] = 1;
							}
							DeleteVertex(v);
							updateDegree();
							break;
						}
					}
					flag = true;
				}
			}
		}
	}
	for (int i = 0; i < init_vCount; i++) {
		for (int j = 0; j < init_vCount; j++) {
			if (Edge[i][j]) {
				now_eCount++; //エッジ数の削減
			}
		}
	}
	for (int i = 0; i < init_vCount; i++) {
		if (inDeg[i] != 0 || outDeg[i] != 0) {
			now_vCount++; //縮小頂点の数
			validVertex[validVertexCount] = i; //縮小後の残りの頂点を記録
			validVertexCount++;
		}
	}
}
