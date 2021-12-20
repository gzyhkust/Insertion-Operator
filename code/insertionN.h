#ifndef INSERTION_H
#define INSERTION_H

#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "PLF.h"



using namespace std;

struct Request {
	int s, e, com;
	double tim, ddl, pr, len;
};


struct Worker {
	int memory = 0;
	int pid, num, cap, gid, vid;
	double tim;
	int rwn;
	int maxnum = 0;
	vector<int> S;
	vector<double> reach;
	vector<int> picked;
	vector<double> ddl;

	vector<PLF> FunctionTimesegments;
 	vector<tuple<int, double>> trajectory;

	inline void pop() {
		if (!S.empty()) {			
			
			if (FunctionTimesegments.size() == 1)
			{
				FunctionTimesegments.erase(FunctionTimesegments.begin());
			}
			else if(FunctionTimesegments.size() != 0)
			{
				FunctionTimesegments.erase(FunctionTimesegments.begin() + S.size()-1);
				FunctionTimesegments.erase(FunctionTimesegments.begin());				
			}

			trajectory.push_back(tuple<int,double>(S[0], reach[0]));
			S.erase(S.begin());
			reach.erase(reach.begin());
			picked.erase(picked.begin());
			ddl.erase(ddl.begin());
		}
	}	

};



// --------------------------------------------- related to grid index
typedef int NodeID;
typedef int EdgeID;
typedef double EdgeWeight;

struct vertex {
	NodeID id;
	double x;
	double y;
};

struct Grid {
	vector<int> taxi;
};

struct Position {
	double x, y;
};

double rad(double d);
double RealDistance(double lat1,double lng1,double lat2,double lng2);
double gridDist(int a, int b);
int getGridID(int x);
int getGridID(double x, double y);
double centerDist(int gid, int pid);
void updateGrid(int pid, int wid, int tag);
void initGrid();
void insertTaxi(int gid, int tid); 
vector<int> single_search(int s, double ddl);
// ---------------------------------------------------------------
void timeDependentInsertion();
int Pos(int x);
double DDLEndPos(int x);
void updateDriver(int i, double t);
void assignTaix();
void finishTaxi(int i);
void readInput();
void assignTaxi(vector<int>& cars);

void try_insertion(Worker &w, int rid, double &delta, int &i, int &j);
int try_insertion_euclidDist(Worker &w, int rid);

void insertion(Worker &w, int rid, int wid, int i, int j);

void updateDriverArr(Worker& w);

double function_arrival(Worker &w, int i, int j,double arr);
void freeMemory();

void recordTrajectory();

void test();
void testquery();
void testplf();
// ---------------------------------------------------------------
// struct QueryResults{
// 	int s;
// 	int d;
// 	vector<vector<double>> breakpoints;

// 	double Delay(){
// 		return 0;
// 	}
// };
#endif //TDGT_IO_H