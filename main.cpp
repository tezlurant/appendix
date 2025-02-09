#define MAIN_TEST_ALL_DATA
#define TEST_SA_COLD
// 追加データの実験ではMAX_Vの変更を忘れずに
//#define ADDITIONAL_DATA

#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>

#include "SA/sa_lib.h"
#include "SA/sa_lib.cpp"

//int Node = 1000;
// 追加データの実験ではMAX_Vの変更を忘れずに


static void run_bench(int seed, int selected_vert_num, std::vector<std::string> *dataFile);

#ifdef MAIN_TEST_ALL_DATA
#ifdef ADDITIONAL_DATA
	#if Node == 500
	std::vector<std::string> dataFile500 = {
		std::string("P500-5000-0.dat"), std::string("P500-5000-1.dat"), std::string("P500-5000-2.dat"),
		std::string("P500-5000-3.dat"), std::string("P500-5000-4.dat"), std::string("P500-5000-5.dat"),
		std::string("P500-5000-6.dat"), std::string("P500-5000-7.dat"), std::string("P500-5000-8.dat"),
		std::string("P500-5000-9.dat"),
		std::string("P500-5500-0.dat"), std::string("P500-5500-1.dat"), std::string("P500-5500-2.dat"),
		std::string("P500-5500-3.dat"), std::string("P500-5500-4.dat"), std::string("P500-5500-5.dat"),
		std::string("P500-5500-6.dat"), std::string("P500-5500-7.dat"), std::string("P500-5500-8.dat"),
		std::string("P500-5500-9.dat"),
		std::string("P500-6000-0.dat"), std::string("P500-6000-1.dat"), std::string("P500-6000-2.dat"),
		std::string("P500-6000-3.dat"), std::string("P500-6000-4.dat"), std::string("P500-6000-5.dat"),
		std::string("P500-6000-6.dat"), std::string("P500-6000-7.dat"), std::string("P500-6000-8.dat"),
		std::string("P500-6000-9.dat"),
		std::string("P500-6500-0.dat"), std::string("P500-6500-1.dat"), std::string("P500-6500-2.dat"),
		std::string("P500-6500-3.dat"), std::string("P500-6500-4.dat"), std::string("P500-6500-5.dat"),
		std::string("P500-6500-6.dat"), std::string("P500-6500-7.dat"), std::string("P500-6500-8.dat"),
		std::string("P500-6500-9.dat"),
		std::string("P500-7000-0.dat"), std::string("P500-7000-1.dat"), std::string("P500-7000-2.dat"),
		std::string("P500-7000-3.dat"), std::string("P500-7000-4.dat"), std::string("P500-7000-5.dat"),
		std::string("P500-7000-6.dat"), std::string("P500-7000-7.dat"), std::string("P500-7000-8.dat"),
		std::string("P500-7000-9.dat")
	};
	#endif
	#if Node == 1000
	std::vector<std::string> dataFile1000 = {
		std::string("P1000-10000-0.dat"), std::string("P1000-10000-1.dat"), std::string("P1000-10000-2.dat"),
		std::string("P1000-10000-3.dat"), std::string("P1000-10000-4.dat"), std::string("P1000-10000-5.dat"),
		std::string("P1000-10000-6.dat"), std::string("P1000-10000-7.dat"), std::string("P1000-10000-8.dat"),
		std::string("P1000-10000-9.dat"),
		std::string("P1000-15000-0.dat"), std::string("P1000-15000-1.dat"), std::string("P1000-15000-2.dat"),
		std::string("P1000-15000-3.dat"), std::string("P1000-15000-4.dat"), std::string("P1000-15000-5.dat"),
		std::string("P1000-15000-6.dat"), std::string("P1000-15000-7.dat"), std::string("P1000-15000-8.dat"),
		std::string("P1000-15000-9.dat"),
		std::string("P1000-20000-0.dat"), std::string("P1000-20000-1.dat"), std::string("P1000-20000-2.dat"),
		std::string("P1000-20000-3.dat"), std::string("P1000-20000-4.dat"), std::string("P1000-20000-5.dat"),
		std::string("P1000-20000-6.dat"), std::string("P1000-20000-7.dat"), std::string("P1000-20000-8.dat"),
		std::string("P1000-20000-9.dat"),
		std::string("P1000-25000-0.dat"), std::string("P1000-25000-1.dat"), std::string("P1000-25000-2.dat"),
		std::string("P1000-25000-3.dat"), std::string("P1000-25000-4.dat"), std::string("P1000-25000-5.dat"),
		std::string("P1000-25000-6.dat"), std::string("P1000-25000-7.dat"), std::string("P1000-25000-8.dat"),
		std::string("P1000-25000-9.dat"),
		std::string("P1000-30000-0.dat"), std::string("P1000-30000-1.dat"), std::string("P1000-30000-2.dat"),
		std::string("P1000-30000-3.dat"), std::string("P1000-30000-4.dat"), std::string("P1000-30000-5.dat"),
		std::string("P1000-30000-6.dat"), std::string("P1000-30000-7.dat"), std::string("P1000-30000-8.dat"),
		std::string("P1000-30000-9.dat")
	};
	#endif
	#if Node == 2000
	std::vector<std::string> dataFile2000 = {
		std::string("P2000-40000-0.dat"), std::string("P2000-40000-1.dat"), std::string("P2000-40000-2.dat"),
		//std::string("P2000-40000-3.dat"), std::string("P2000-40000-4.dat"), std::string("P2000-40000-5.dat"),
		//std::string("P2000-40000-6.dat"), std::string("P2000-40000-7.dat"), std::string("P2000-40000-8.dat"),
		//std::string("P2000-40000-9.dat"),
		std::string("P2000-60000-0.dat"), std::string("P2000-60000-1.dat"), std::string("P2000-60000-2.dat"),
		//std::string("P2000-60000-3.dat"), std::string("P2000-60000-4.dat"), std::string("P2000-60000-5.dat"),
		//std::string("P2000-60000-6.dat"), std::string("P2000-60000-7.dat"), std::string("P2000-60000-8.dat"),
		//std::string("P2000-60000-9.dat"),
		std::string("P2000-80000-0.dat"), std::string("P2000-80000-1.dat"), std::string("P2000-80000-2.dat"),
		//std::string("P2000-80000-3.dat"), std::string("P2000-80000-4.dat"), std::string("P2000-80000-5.dat"),
		//std::string("P2000-80000-6.dat"), std::string("P2000-80000-7.dat"), std::string("P2000-80000-8.dat"),
		//std::string("P2000-80000-9.dat")
	};
	#endif
	#if Node == 3000
	std::vector<std::string> dataFile3000 = {
		std::string("P3000-60000-0.dat"), std::string("P3000-60000-1.dat"), std::string("P3000-60000-2.dat"),
		//std::string("P3000-60000-3.dat"), std::string("P3000-60000-4.dat"), std::string("P3000-60000-5.dat"),
		//std::string("P3000-60000-6.dat"), std::string("P3000-60000-7.dat"), std::string("P3000-60000-8.dat"),
		//std::string("P3000-60000-9.dat"),
		std::string("P3000-90000-0.dat"), std::string("P3000-90000-1.dat"), std::string("P3000-90000-2.dat"),
		//std::string("P3000-90000-3.dat"), std::string("P3000-90000-4.dat"), std::string("P3000-90000-5.dat"),
		//std::string("P3000-90000-6.dat"), std::string("P3000-90000-7.dat"), std::string("P3000-90000-8.dat"),
		//std::string("P3000-90000-9.dat"),
		std::string("P3000-120000-0.dat"), std::string("P3000-120000-1.dat"), std::string("P3000-120000-2.dat"),
		//std::string("P3000-120000-3.dat"), std::string("P3000-120000-4.dat"), std::string("P3000-120000-5.dat"),
		//std::string("P3000-120000-6.dat"), std::string("P3000-120000-7.dat"), std::string("P3000-120000-8.dat"),
		//std::string("P3000-120000-9.dat")
	};
	#endif
	#if Node == 4000
	std::vector<std::string> dataFile4000 = {
		std::string("P4000-80000-0.dat"), std::string("P4000-80000-1.dat"), std::string("P4000-80000-2.dat"),
		//std::string("P4000-80000-3.dat"), std::string("P4000-80000-4.dat"), std::string("P4000-80000-5.dat"),
		//std::string("P4000-80000-6.dat"), std::string("P4000-80000-7.dat"), std::string("P4000-80000-8.dat"),
		//std::string("P4000-80000-9.dat"),
		std::string("P4000-120000-0.dat"), std::string("P4000-120000-1.dat"), std::string("P4000-120000-2.dat"),
		//std::string("P4000-120000-3.dat"), std::string("P4000-120000-4.dat"), std::string("P4000-120000-5.dat"),
		//std::string("P4000-120000-6.dat"), std::string("P4000-120000-7.dat"), std::string("P4000-120000-8.dat"),
		//std::string("P4000-120000-9.dat"),
		std::string("P4000-160000-0.dat"), std::string("P4000-160000-1.dat"), std::string("P4000-160000-2.dat"),
		//std::string("P4000-160000-3.dat"), std::string("P4000-160000-4.dat"), std::string("P4000-160000-5.dat"),
		//std::string("P4000-160000-6.dat"), std::string("P4000-160000-7.dat"), std::string("P4000-160000-8.dat"),
		//std::string("P4000-160000-9.dat")
	};
	#endif
	#if Node == 5000
	std::vector<std::string> dataFile5000 = {
		std::string("P5000-100000-0.dat"), std::string("P5000-100000-1.dat"), std::string("P5000-100000-2.dat"),
		//std::string("P5000-100000-3.dat"), std::string("P5000-100000-4.dat"), std::string("P5000-100000-5.dat"),
		//std::string("P5000-100000-6.dat"), std::string("P5000-100000-7.dat"), std::string("P5000-100000-8.dat"),
		//std::string("P5000-100000-9.dat"),
		std::string("P5000-150000-0.dat"), std::string("P5000-150000-1.dat"), std::string("P5000-150000-2.dat"),
		//std::string("P5000-150000-3.dat"), std::string("P5000-150000-4.dat"), std::string("P5000-150000-5.dat"),
		//std::string("P5000-150000-6.dat"), std::string("P5000-150000-7.dat"), std::string("P5000-150000-8.dat"),
		//std::string("P5000-150000-9.dat"),
		std::string("P5000-200000-0.dat"), std::string("P5000-200000-1.dat"), std::string("P5000-200000-2.dat"),
		//std::string("P5000-200000-3.dat"), std::string("P5000-200000-4.dat"), std::string("P5000-200000-5.dat"),
		//std::string("P5000-200000-6.dat"), std::string("P5000-200000-7.dat"), std::string("P5000-200000-8.dat"),
		//std::string("P5000-200000-9.dat")
	};
	#endif
#else
	#if Node==50
	std::vector<std::string> dataFile50 = {
		std::string("P50-100.dat"), std::string("P50-150.dat"), std::string("P50-200.dat"),
		std::string("P50-250.dat"), std::string("P50-300.dat"), std::string("P50-500.dat"),
		std::string("P50-600.dat"), std::string("P50-700.dat"), std::string("P50-800.dat"),
		std::string("P50-900.dat")
	};
	#endif
	#if Node==100
	std::vector<std::string> dataFile100 = {
		std::string("P100-200.dat"), std::string("P100-300.dat"), std::string("P100-400.dat"),
		std::string("P100-500.dat"), std::string("P100-600.dat"), std::string("P100-1000.dat"),
		std::string("P100-1100.dat"), std::string("P100-1200.dat"), std::string("P100-1300.dat"),
		std::string("P100-1400.dat")
	};
	#endif
	#if Node==500
	std::vector<std::string> dataFile500 = {
		std::string("P500-1000.dat"), std::string("P500-1500.dat"), std::string("P500-2000.dat"),
		std::string("P500-2500.dat"), std::string("P500-3000.dat"), std::string("P500-5000.dat"),
		std::string("P500-5500.dat"), std::string("P500-6000.dat"), std::string("P500-6500.dat"),
		std::string("P500-7000.dat")
	};
	#endif
	#if Node==1000
	std::vector<std::string> dataFile1000 = {
		std::string("P1000-3000.dat"), std::string("P1000-3500.dat"), std::string("P1000-4000.dat"),
		std::string("P1000-4500.dat"), std::string("P1000-5000.dat"), std::string("P1000-10000.dat"),
		std::string("P1000-15000.dat"), std::string("P1000-20000.dat"), std::string("P1000-25000.dat"),
		std::string("P1000-30000.dat")
	};
	#endif
#endif

// TODO: need to change the graph reading library
std::vector<std::string> ourGraphs = {
	std::string("genome1.in"), std::string("genome2.in"), std::string("tpcc.in"), std::string("tpcc30.in")
};


// %.warm_MT15000    : HotD ::= 1
// %.warm_MT15000    : ColdD ::= 1
// %.warm_MT15000    : ProbHot ::= 0.00585
// %.warm_MT15000    : ProbCold ::= 0.00006

#if defined(TEST_SA_COLD)
	const float p_sa_cold_p1 = 0.00001;
	const int   p_sa_cold_d1 = 1;
	const float p_sa_cold_p2 = 0.0000001;
	const int   p_sa_cold_d2 = 1;
	const int   p_sa_cold_iter = 1<<25;
	const int   p_sa_cold_reps = 1<<25;
#endif

#ifdef TEST_SA_COLD
	FILE *fout2, *fout2_prof;
#endif

int main()
{

	// TODO: our runs are not deterministic
	//int seed = static_cast<unsigned int>(time(0));
	//std::cout << "seed = " << seed << std::endl;
#ifdef TEST_SA_COLD
	char profFileName[1<<10];
	#ifdef ADDITIONAL_DATA
		sprintf(profFileName, "./SA_Result/SAColdRunTimeAdditionalP%d.txt", Node);
	#else
		sprintf(profFileName, "./SA_Result/SAColdRunTimeP%d.txt", Node);
	#endif
	fout2 = fopen(profFileName, "w");
	fprintf(fout2, "seed = %d\n", 0);
#endif
	std::vector<std::string>* dataFile = nullptr;
	#ifdef ADDITIONAL_DATA
		#if Node==500
			dataFile = &dataFile500;
		#endif
		#if Node==1000
			dataFile = &dataFile1000;
		#endif
		#if Node==2000
			dataFile = &dataFile2000;
		#endif
		#if Node==3000
			dataFile = &dataFile3000;
		#endif
		#if Node==4000
			dataFile = &dataFile4000;
		#endif
		#if Node==5000
			dataFile = &dataFile5000;
		#endif
	#else
		#if Node==50
			dataFile = &dataFile50;
		#endif
		#if Node==100
			dataFile = &dataFile100;
		#endif
		#if Node==500
			dataFile = &dataFile500;
		#endif
		#if Node==1000
			dataFile = &dataFile1000;
		#endif
	#endif

	//その他
	//dataFile = &ourGraphs;
	run_bench(0, Node, dataFile);

#ifdef TEST_SA_COLD
	fclose(fout2);
#endif
	return 0;
}
#endif

static void run_bench(int seed, int selected_vert_num, std::vector<std::string> *dataFile)
{

	const int runTimes = 30;
	TIMER_T startT, endT;

	for (int i = 0; i < dataFile->size(); ++i)
	{
		char profFileName[1<<10];
		#ifdef ADDITIONAL_DATA
			#if Node==500
				std::string folder("./fsp-data/P500/");
			#endif
			#if Node==1000
				std::string folder("./fsp-data/P1000/");
			#endif
			#if Node==2000
				std::string folder("./fsp-data/P2000/");
			#endif
			#if Node==3000
				std::string folder("./fsp-data/P3000/");
			#endif
			#if Node==4000
				std::string folder("./fsp-data/P4000/");
			#endif
			#if Node==5000
				std::string folder("./fsp-data/P5000/");
			#endif
		#else
			std::string folder("./fsp-data/");
		#endif
		//std::string folder("./our-graphs/");
		std::cout << (*dataFile)[i] << std::endl;

#ifdef TEST_SA_COLD
		{
			SA_parameters_s p {
				.n      = p_sa_cold_iter,
  				.repOps = p_sa_cold_reps,
  				.hot    = p_sa_cold_p1,
				.hotD   = p_sa_cold_d1,
				.cold   = p_sa_cold_p2,
				.coldD  = p_sa_cold_d2
			};

			sprintf(profFileName, "./SA_Result/SA_%s_prof.txt", (*dataFile)[i].c_str());
			fout2_prof = fopen(profFileName, "w");

			//std::cout << "-SA COLD" << std::endl;

			SA_init_F(p, (folder+(*dataFile)[i]).c_str());
			SA_set_prof_file(fout2_prof);

			double totalTime = 0;
			long totalFvsNum = 0;
			long minFvs = SA_get_nbVerts();
			long maxFvs = -1;
			
			std::vector<int> allFvs;
			fprintf(fout2, "%s,", (*dataFile)[i].c_str());
			for (int j = 0; j < runTimes; ++j)
			{
				//std::cout << "#" << j << " ";

				SA_reset(); // TODO: not working

				TIMER_READ(startT);
				SA_run();
				TIMER_READ(endT);

				totalTime += CLOCK_DIFF_MS(startT, endT);

				//int *P = SA_getBestSolution();
				//while (-1 != *P) P++;
				//long fvsNum = SA_get_nbVerts() - (P - SA_getBestSolution()) + SA_get_nbLoop();
				long fvsNum = SA_get_nbVerts() - SA_maxDAG() + SA_get_nbLoop();
				//std::cout<<"fvs="<<fvsNum<<std::endl;

				totalFvsNum += fvsNum;
				if (fvsNum > maxFvs)
					maxFvs = fvsNum;
				if (fvsNum < minFvs)
					minFvs = fvsNum;

				allFvs.push_back(fvsNum);
				//SA_destroy();

				fprintf(fout2, "%lu,", fvsNum);

			}
			fprintf(fout2, "\n");
			fflush(fout2);
			//SA_init_F(p, (folder+(*dataFile)[i]).c_str());
			double deviation = 0.0;
			double avg = (double)totalFvsNum / runTimes;
			printf("avg=%lf\n\n", avg);
			for (int k : allFvs)
			{
				deviation += ((k - (double)totalFvsNum / runTimes) * (k - (double)totalFvsNum / runTimes));
			}
			deviation /= allFvs.size();
			deviation = sqrt(deviation);

			// fprintf(fout2, "%d;%d;%lf;%lf;%ld;%ld;%lf\n",
			// 	SA_get_nbVerts(), SA_get_nbEdges(), (double)totalTime / runTimes,
			// 	(double)totalFvsNum / runTimes, minFvs, maxFvs, deviation);
			// fflush(fout2);
			fclose(fout2_prof);

			SA_destroy();
		}
#endif
	}
}
