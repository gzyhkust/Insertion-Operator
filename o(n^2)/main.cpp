#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include "TDGT.h"
#include "io.h"
#include "TDDijkstra.h"

using namespace std;

// extern string requestFile;
// extern string workerFile;

extern const int MAX_NODE ;
int nV, m, c, n;
double gridL, alpha;

int pos, cnt, mcnt;
//pos = 0;

double EPS = 0.001, INF = 100000;
// extern double gridL, alpha;

int Pos(int x);

TDGT tree;

struct Request {
	int s, e, com;
	double tim, ddl, pr, len;
};
Request* R;

struct Detour
{
	// arrival time at x-th is arx, then the arrival time at y-th in sschedule is arry
	double detour(int x, int y, double arrx, vector<int>& sschedule){
		if( x == y){return arrx;}
		int k = x+1;
		double arry = arrx + tree.tdsp( Pos(sschedule[x]),Pos(sschedule[k]),arrx,false);
		while (k<=y)
		{
			arry += tree.tdsp( Pos(sschedule[k-1]),Pos(sschedule[k]),arry,false);
			k++;
		}
		
		return arry;
	}
};

struct DDL
{	
	// arrival time at i is arri then check if it is feasible or not
	//int ddl_arrival(int i, double arri, vector<int>& schedule, vector<vector<int>>& functionindex, vector<Detour>& detours){
	int ddl_arrival(int i, double arri, vector<int>& schedule, int functionindex[100][100], vector<Detour>& detours){

		if(schedule[i] & 1){
			if(arri > R[(schedule[i])>>1].ddl){return -1;}

			if (i<schedule.size()-1){
				int arr_i_1 = detours[functionindex[i][i+1]].detour(i,i+1,arri,schedule);
				if (ddl_arrival(i+1, arr_i_1, schedule, functionindex, detours)  == -1){return -1;}  
			}
		}
		else
		{
			int arr_i_1 = detours[functionindex[i][i+1]].detour(i,i+1,arri,schedule);
			if (ddl_arrival(i+1, arr_i_1, schedule, functionindex, detours)  == -1){return -1;}  
		}
		
		return 1;
	}
};


struct Worker {
	int pid, num, cap, gid, vid;
	double tim;
	int rwn;
	vector<int> S;
	vector<double> reach;
	vector<int> picked;
	vector<Detour> detours;
 	//vector<vector<int>> functionindex;
	int functionindex[100][100];
	vector<DDL> ddlfuncs;
	//vector<double> arrddl;
	//vector<double> slack;

	inline void pop() {
		if (!S.empty()) {
			S.erase(S.begin());
			reach.erase(reach.begin());
			picked.erase(picked.begin());
			ddlfuncs.erase(ddlfuncs.begin());
			//slack.erase(slack.begin());
		}
	}
};
Worker* W;


// ------------------- related to grid index
typedef int NodeID;
typedef int EdgeID;
typedef double EdgeWeight;

struct vertex {
	NodeID id;
	double x;
	double y;
};

//vector<vertex>& vertices;

struct Grid {
	vector<int> taxi;
};
Grid* gr = NULL;
int* anchor = NULL;
int graph_len, graph_wid, grid_len, grid_wid;
int grid_num_per_row, grid_num_per_col, grid_sz;
double mnx = INF, mxx = -INF, mny = INF, mxy = -INF;
vector<bool> visitGrid;
int gridm = 0;

double gridDist(int a, int b);
int getGridID(int x);
int getGridID(double x, double y);
double centerDist(int gid, int pid);
void updateGrid(int pid, int wid, int tag);
void initGrid();
void insertTaxi(int gid, int tid); 
// ---------------------------------------------------------------
void timeDependentInsertion();
double TDGTPathQuery(int a, int b, double t);
void updateDriver(int i, double t);
void assignTaix();
void finishTaxi(int i);
void readInput();
vector<int> single_search(int s, double ddl);
void assignTaxi(vector<int>& cars);

void try_insertion(Worker &w, int rid, double &delta, int &i, int &j);
void try_insertion_n2(Worker &w, int rid, double &delta, int &i, int &j);

void insertion(Worker &w, int rid, int wid, int i, int j);
void insertion_n2(Worker &w, int rid, int wid, int i, int j);
void updateDriverArr(Worker& w);
void updateDriverArr_n2(Worker& w, int i, int j);
void freeMemory();
// ---------------------------------------------------------------

//
int ServiedRequest = 0;
int InsertionNum = 0;
//

int main(int argc, char **args) {
    
	//read file to initialize the worker/request
    readInput();
	//initGrid();

	cout << "Start sharing..." << endl;
	clock_t tbegin = clock();
	timeDependentInsertion();
	clock_t tend = clock();
	cout << "Finish sharing \t time cost: " + to_string((tend - tbegin) / CLOCKS_PER_SEC) + "s" << endl;

	cout << " ServiedRequest: " << ServiedRequest  <<"  "<< "InsertionNum: "<< InsertionNum <<endl;
	freeMemory();

	return 0;
}

void freeMemory() {
	delete[] R;
	delete[] W;
	// delete[] gr;
	// delete[] anchor;
}

void readInput() {
	
    string requestFile = "/Users/gongcengyang/Desktop/Research/InsertionTimeDepedent/Insertion/TestMain/Data/Fakerequests.txt";
    string workerFile = "/Users/gongcengyang/Desktop/Research/InsertionTimeDepedent/Insertion/TestMain/Data/Fakeworkers.txt";
	//string verticesFile = "/Users/gongcengyang/Desktop/Research/InsertionTimeDepedent/Insertion/TestMain/Data/chengdu_sub_10000/chengdu_sub_co.txt";

    cout << "begin worker" << endl;
	ifstream other;
	other.open(workerFile.c_str());
	if (!other.is_open()) {
		printf("%s does not exist\n", workerFile.c_str());
		exit(0);
	}
	double ddl, pr;
	other >> m >> c >> gridL >> alpha;
	W = new Worker[m];
	for (int i = 0; i < m; ++ i) {
		other >> W[i].pid >> W[i].cap;
		W[i].num = 0;		
	}
	other >> ddl >> pr;
	other.close();
    cout << "end worker" << endl;


	//cout << m << " " << n << " " << c << " " << gridL << " " << alpha << endl;

	string writepath = "/Users/gongcengyang/Desktop/Research/InsertionTimeDepedent/Insertion/TestMain/Data/chengdu_sub_10000/tree.txt";
    ifstream readtree;
    readtree.open(writepath);
    read_TDGT(readtree, tree);

    cout << "begin request" << endl;
	ifstream ifs;
	ifs.open(requestFile.c_str());
	if (!ifs.is_open()) {
		printf("%s does not exist\n", requestFile.c_str());
		exit(0);
	}
	ifs >> n;
	R = new Request[n];
	for (int i = 0; i < n; ++ i) {
		ifs >> R[i].tim >> R[i].s >> R[i].e >> R[i].com;
		//R[i].len = 100;
		R[i].len  = tree.tdsp(R[i].s,R[i].e, R[i].tim, false);
		R[i].ddl = R[i].tim + R[i].len + 600;
		R[i].pr = pr;
	}
	ifs.close();
	cout << "end request" << endl;	

}


void timeDependentInsertion(){

	vector<int> car;
	for(int i = 0; i<m; i++){
		car.push_back(i);
	}
	while (pos<n)
	{	
		for (int i = 0; i < m; ++i){
			updateDriver(i, R[pos].tim);
		}

		// vector<int> car;
		// car = single_search(R[pos].s, R[pos].ddl - R[pos].tim);
		cout << "assignTaxi:" << pos << endl;
		assignTaxi(car);
	}
	for (int i = 0; i < m; ++i){
		finishTaxi(i);
	}
}

void updateDriver(int i, double t){
	Worker& w = W[i];
	while (w.S.size() > 0 && w.tim < t)
	{	
		double tmp = w.reach[0] - w.tim; // == dist(w.pid, Pos(w.S[0]))

		w.tim += tmp;
		w.pid = Pos(w.S[0]);

		if (w.S[0] & 1) {
			w.num -= R[w.S[0] >> 1].com;
			mcnt ++;
		} else {
			w.num += R[w.S[0] >> 1].com;
		}

		//update vector<Detour> detours; vector<vector<int>> functionindex;
		vector<Detour> detours;
		int functionindex[100][100];
		for (int x = 0; x < 100; x++)
		{
			for(int y = 0; y < 100; y++){
				functionindex[x][y] = 0;
			}
		}
		
		int k = 0;
		for(int x = 1; x < w.S.size(); ++x){
			for(int y = x+1; y < w.S.size(); ++y){
				w.functionindex[x-1][y-1] = k;

				detours.push_back(w.detours[w.functionindex[x][y]]);
				k++;
			}
		}
		
		w.detours.clear();
		w.detours = detours;
		//w.functionindex.clear();
		//w.functionindex = functionindex;

		w.pop();
	}
	if (w.tim < t)
	{
		w.tim = t;
	}
	
	
}

int Pos(int x) {
	return (x&1) ? R[x>>1].e : R[x>>1].s;
}


void assignTaxi(vector<int>& car) {
	double opt = INF, fit;
	int sz = car.size(), id = -1;
	int optimal_i, optimal_j, x, y;
	
	//tot_enum = sz, true_enum = 0;
	if (R[pos].len < INF) {
		for (int i = 0; i < sz; ++ i) {
			//true_enum ++;
			fit = INF;

			try_insertion_n2(W[car[i]], pos, fit, x, y);

			if (fit < INF) {
				if (opt > fit) {
					opt = fit;
					id = car[i];

					optimal_i = x;
					optimal_j = y;
				}
			}
		}	
	}
	
	if (id > -1) {
		ServiedRequest ++;
		insertion_n2(W[id], pos, id, optimal_i, optimal_j);
	}
	pos++;
}

/*
----------------------------------------- n^2 method -----------------------------
*/
void try_insertion_n2(Worker &w, int rid, double &delta, int &optimanl_i, int &optimanl_j){
	Request& r = R[rid];
	double opt = INF;
	
	delta = INF;
	if (w.S.empty()) {
		double tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
		tmp += tree.tdsp(r.s, r.e, tmp,false);
		if (tmp < r.ddl + EPS && r.com <= w.cap) {
			delta = tmp, optimanl_i = 0, optimanl_j = 1;
		}
		return;
	}

	vector<int>& picked = w.picked;
	vector<int>& schedule = w.S;
	vector<double>& reach = w.reach;
	double tmp = 0;

	//n^2 enumerate	
	for(int i = 0; i <= w.S.size(); ++i){
		for (int j = i; j <= w.S.size(); ++j){

			int ddl_flag = 1;

			if(i == 0 && j == 0){
				
				tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
				tmp += tree.tdsp(r.s, r.e, tmp, false);
				tmp += tree.tdsp(r.e, Pos(schedule[0]), tmp,false);
				
				ddl_flag = w.ddlfuncs[i].ddl_arrival(i,tmp, w.S, w.functionindex, w.detours);
				if(ddl_flag == -1 || w.num+r.com > w.cap){break;}
				// if feasibel, try to get the detour
				opt = tmp - reach[i];
				if(opt < delta){delta = opt; optimanl_i = i; optimanl_j = j;}
			}
			else if(i == w.S.size()){

				tmp = reach[w.reach.size() - 1] + tree.tdsp(Pos(schedule[w.S.size() - 1]), r.s, reach[w.S.size() - 1],false);
				tmp += tree.tdsp(r.s, r.e, tmp,false);
				if(tmp > r.ddl || picked[i-1]+r.com > w.cap){ddl_flag = -1; break;}
				opt = tmp - reach[w.S.size() - 1];
				if(opt < delta){delta = opt; optimanl_i = i; optimanl_j = j;}

			}
			else
			{
				InsertionNum ++;

				double det_i, det_j;

				if(i == 0){ // arrival time at r.s
					tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
					if(w.num+r.com > w.cap){ddl_flag = -1; break;} 
				}else
				{
					tmp = reach[i-1] + tree.tdsp(Pos(schedule[i-1]), r.s, reach[i-1],false);
					if(picked[i-1]+r.com > w.cap){ddl_flag = -1; break;}
				}
				if (i == j)
				{
					//1. arrival time at j
					tmp += tree.tdsp(r.s,r.e, tmp,false);
					tmp += tree.tdsp(r.e, Pos(schedule[i]), tmp,false);
					ddl_flag = w.ddlfuncs[i].ddl_arrival(i,tmp, w.S, w.functionindex, w.detours);
					if(ddl_flag == -1){break;}
					opt = tmp - reach[i];
					if(opt < delta){delta = opt; optimanl_i = i; optimanl_j = j;}	
				}else
				{
					//1. arrival time at i
					tmp += tree.tdsp(r.s, Pos(schedule[i]), tmp,false);
					ddl_flag = w.ddlfuncs[i].ddl_arrival(i,tmp, w.S, w.functionindex, w.detours);
					if(ddl_flag == -1){break;}

					//2. arrival time at j
					if (j == w.S.size())
					{
						tmp += w.detours[w.functionindex[i][j-1]].detour(i,j-1,tmp, w.S);
						tmp += tree.tdsp(Pos(schedule[j-1]), r.e, tmp, false);
						if( tmp <= r.ddl){
							opt = tmp - reach[j-1];
							if(opt < delta){delta = opt; optimanl_i = i; optimanl_j = j;}							 
						}
					}else
					{
						tmp += w.detours[w.functionindex[i][j-1]].detour(i,j-1,tmp, w.S);
						tmp += tree.tdsp(Pos(schedule[j-1]), r.e, tmp, false);
						tmp += tree.tdsp(r.e, Pos(schedule[j]) , tmp, false);
						ddl_flag = w.ddlfuncs[i].ddl_arrival(j,tmp, w.S, w.functionindex, w.detours);
						if(ddl_flag == -1){break;}
						opt = tmp - reach[j];
						if(opt < delta){delta = opt; optimanl_i = i; optimanl_j = j;}
					}
							
				}
	
			}	
		}			
	}
}
 

void insertion_n2(Worker &w, int rid, int wid, int optimal_i, int optimal_j){
	if (w.S.empty())
	{
		w.S.push_back(rid << 1);
		w.S.push_back(rid << 1 | 1);
		updateDriverArr_n2(w,optimal_i,optimal_j);
		return;
	}

	vector<int> ret;

	ret.clear();

	if(optimal_i<w.S.size()){
		if(optimal_j<w.S.size()){
			for(int k = 0; k < w.S.size(); ++k){
				if (k == optimal_i && k < optimal_j)
				{
					ret.push_back(rid << 1);
					ret.push_back(w.S[k]);
				}
				else if (k == optimal_i && k == optimal_j)
				{
					ret.push_back(rid << 1);
					ret.push_back(rid << 1|1);
					ret.push_back(w.S[k]);
				}
				else if(k > optimal_i && k == optimal_j)
				{
					ret.push_back(rid << 1|1);
					ret.push_back(w.S[k]);
				}else
				{
					ret.push_back(w.S[k]);
				}				
			}
		}
		else
		{
			for(int k = 0; k<w.S.size(); ++k){
				if(k == optimal_i){
					ret.push_back(rid << 1);
					ret.push_back(w.S[k]);							
				}
				else
				{
					ret.push_back(w.S[k]);
				}
			}
			ret.push_back(rid << 1|1);						
		}	
	}
	else
	{
		for(int k = 0; k<w.S.size(); ++k){
			ret.push_back(w.S[k]);
		}
		ret.push_back(rid << 1);
		ret.push_back(rid << 1|1);
	}

	w.S.clear();
	w.S = ret;
	updateDriverArr_n2(w,optimal_i,optimal_j);
}


void updateDriverArr_n2(Worker& w, int optimal_i, int optimal_j){

	//cout << 'updateDriverArr_n2' << endl;
	double tim = w.tim;
	vector<double>& reach = w.reach; 
	vector<Detour>& detours = w.detours;
	vector<DDL>& ddls = w.ddlfuncs;
	
	/////
	reach.clear();
	cout << w.S << endl;
	for(int k = 0;k < w.S.size(); ++k){
		if (k == 0){
			tim += tree.tdsp(w.pid, Pos(w.S[k]), tim,false);
		}
		else{
			tim += tree.tdsp(Pos(w.S[k-1]), Pos(w.S[k]), tim,false);
		}
		reach.push_back(tim);		
	}

	//functionindex.clear();
	for (int i = 0; i < 100; i++)
	{
		for (int j = 0; j < 100; j++)
		{
			w.functionindex[i][j] = 0;
		}
		
	}
	
	detours.clear();
	int k = 0;
	for(int i = 0; i < w.S.size(); ++i){
		for(int j = i+1; j< w.S.size(); ++j){
			w.functionindex[i][j] = k;
			Detour detour;
			detour.detour(i,j, reach[i], w.S);
			detours.push_back(detour);
			k++;
		}
	}

	for(int i = 0; i<w.S.size(); ++i){
		//DDL ddl = new DDL();
		DDL ddl;
		ddl.ddl_arrival(i, reach[i], w.S, w.functionindex, detours);
		ddls.push_back(ddl);
	}
	////
	vector<int>& picked = w.picked;
	picked.clear();
	int cc = w.num;
	for(int k = 0; k < w.S.size(); ++k){
		if (w.S[k] & 1) {
			cc -= R[w.S[k] >> 1].com;
		} else {
			cc += R[w.S[k] >> 1].com;
		}
		picked.push_back(cc);		
	}
}

void finishTaxi(int i) {
	Worker& w = W[i];
	
	for (int i=0; i<w.S.size(); ++i) {
		double tmp = w.reach[i] - w.tim;
		///ans += alpha * tmp;
		w.tim += tmp;
		if (w.S[i] & 1) {
			mcnt ++;
		}
	}
}

// void updateDriverArr_n2(Worker& w){
// 	double tim = w.tim;
// 	vector<double>& reach = w.reach; 
// 	vector<Detour>& detous = w.detours;
// 	vector<vector <int>>& functionindex = w.functionindex;
// 	vector<DDL>& ddl = w.ddlfuncs;
// 	vector<int>& picked = w.picked;
// 	/////

// 	////

// 	reach.clear();
// 	for(int k = 0; k < w.S.size(); ++k){
// 		if (k == 0){
// 			tim += tree.tdsp(w.pid, Pos(w.S[k]), tim,false);
// 		}
// 		else{
// 			tim += tree.tdsp(Pos(w.S[k-1]), Pos(w.S[k]), tim,false);
// 		}
// 		reach.push_back(tim);
// 	}

// 	vector<int>& picked = w.picked;
// 	picked.clear();
// 	int cc = w.num;
// 	for(int k = 0; k < w.S.size(); ++k){
// 		if (w.S[k] & 1) {
// 			cc -= R[w.S[k] >> 1].com;
// 		} else {
// 			cc += R[w.S[k] >> 1].com;
// 		}
// 		picked.push_back(cc);		
// 	}

// 	int k = 0;
// 	//Detour* functions = w.detours;

// 	vector<vector<int>> functionindex_ij;

// 	w.functionindex.clear();
// 	for(int i = 0; i<w.S.size(); ++i){
// 		for(int j = i+1; j<w.S.size(); ++j){
// 			functionindex_ij[i][j] = k;
// 			k++;
// 		}
// 	}
// 	w.functionindex = functionindex_ij;

// 	Detour* detours_ij = new Detour[k];
// 	vector<int> schedule;
// 	k = 0;
// 	for(int i = 0; i<w.S.size(); ++i){
// 		for(int j = i+1; j<w.S.size(); ++j){
// 			detours_ij[k].detour(i,j, reach[i], w.S);
// 		}
// 	}
// 	w.detours.clear();
// 	w.detours = detours_ij;	
// }

/*
  ---------------------------------------  naitive method ----------------------------------
*/
// void insertion(Worker &w, int rid, int wid, int optimal_i, int optimal_j){
// 	if (w.S.empty())
// 	{
// 		w.S.push_back(rid << 1);
// 		w.S.push_back(rid << 1 | 1);
// 		updateDriverArr(w);
// 		return;
// 	}

// 	vector<int>& picked = w.picked;
// 	vector<int> ret = w.S;
// 	vector<double>& reach = w.reach;

// 	ret.clear();
// 	for(int k=0; k < w.S.size(); ++k){
// 		if(k<optimal_i){ret.push_back(w.S[k]);}
// 		else if (k == optimal_i && k < optimal_j)
// 		{
// 			ret.push_back(rid << 1);
// 			ret.push_back(w.S[k]);
// 		}
// 		else if (k == optimal_i && k == optimal_j)
// 		{
// 			ret.push_back(rid << 1);
// 			ret.push_back(rid << 1|1);
// 			ret.push_back(w.S[k]);
// 		}
// 		else if(k > optimal_i && k == optimal_j)
// 		{
// 			ret.push_back(rid << 1|1);
// 		}else
// 		{
// 			ret.push_back(w.S[k]);
// 		}	
// 	}
//	w.S = ret;
// 	updateDriverArr(w);
// }

// void updateDriverArr(Worker& w){
// 	double tim = w.tim;
// 	vector<double>& reach = w.reach;

// 	reach.clear();
// 	for(int k = 0; k < w.S.size(); ++k){
// 		if (k == 0){
// 			tim += tree.tdsp(w.pid, Pos(w.S[k]), tim,false);
// 		}
// 		else{
// 			tim += tree.tdsp(Pos(w.S[k-1]), Pos(w.S[k]), tim,false);
// 		}
// 		reach.push_back(tim);
// 	}

// 	vector<int>& picked = w.picked;
// 	picked.clear();
// 	int cc = w.num;
// 	for(int k = 0; k < w.S.size(); ++k){
// 		if (w.S[k] & 1) {
// 			cc -= R[w.S[k] >> 1].com;
// 		} else {
// 			cc += R[w.S[k] >> 1].com;
// 		}
// 		picked.push_back(cc);		
// 	}

// 	//
// 	cout << "schedule :" << w.S <<endl;
// 	for(int i = 0; i<w.S.size(); ++i){
// 		cout << Pos(w.S[i])<<endl;
// 	}
// 	cout << "reach :" << w.reach <<endl;
// 	cout << "picked :" << w.picked <<endl;
// 	//cout << w.ca
// 	cout << "------------" << endl;
// 	//break;
// 	//		
// }





/////////////////////////////////////////////// basic method backup
/*
int DDLFlagEndK(Worker &w, int k ,double tmp){
	for(int i = 0; i< w.picked.size(); ++i){
		Request& r_k = R[w.picked[i]];
		if (r_k.e == Pos(w.S[k])){
			if (r_k.ddl<tmp){
				return -1;	
			}
		}
	}
	return 1;
}

void try_insertion(Worker &w, int rid, double &delta, int &optimanl_i, int &optimanl_j) {
	Request& r = R[rid];
	double opt = INF;
	
	delta = INF;
	if (w.S.empty()) {
		double tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false) + tree.tdsp(r.s, r.e, w.tim,false);
		if (tmp < r.ddl + EPS && r.com <= w.cap) {
			delta = tmp, optimanl_i = 0, optimanl_j = 1;
		}
		return;
	}

	vector<int>& picked = w.picked;
	vector<int>& schedule = w.S;
	vector<double>& reach = w.reach;
	double tmp = 0;
	for (int i = 0; i <= w.S.size(); ++ i){
		for (int j = i; j<= w.S.size(); ++j){

			int ddl_flag = 1;

			if (i == 0 && j ==0){

				tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false); 
				tmp += tree.tdsp(r.s, r.e, tmp,false);
				if (tmp<=r.ddl && w.num+r.com<=w.cap){ //at least w can pick r
					tmp += tree.tdsp(r.e, schedule[0], tmp,false);
		
					for(int k = 1; k < w.S.size(); ++k){
						tmp += tree.tdsp(schedule[k-1], schedule[k], tmp, false);
						if (w.S[k] & 1){			  // satisfy all request end at k
							ddl_flag = DDLFlagEndK(w, k, tmp);
							if(ddl_flag == -1){break;}
						}
					}
				}else
				{
					ddl_flag = -1;
					break;
				}
				if(tmp<delta){delta = tmp, optimanl_i = i, optimanl_j = j;}
			
			}else if (i == w.S.size()){
				tmp = reach[w.reach.size() - 1] + tree.tdsp(schedule[w.S.size() - 1], r.s, tmp, false);
				tmp = tmp + tree.tdsp(r.s, r.e, tmp,false);
				if (tmp > r.ddl || picked[i-1]+r.com > w.cap){ddl_flag = -1; break;}
				if(tmp<delta){delta = tmp, optimanl_i = i, optimanl_j = j;}
			}
			else{

				if(i == 0){ // arrival time r.s after inserting it 
					tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
					if(w.num+r.com > w.cap){ddl_flag = -1; break;} 
				}else
				{
					tmp = reach[i-1] + tree.tdsp(schedule[i-1], r.s, reach[i-1],false);
					if(picked[i-1]+r.com > w.cap){ddl_flag = -1; break;}
				}

				if(i == j){
					tmp += tree.tdsp(r.s, r.e, tmp,false);
					if(tmp > r.ddl){ddl_flag = -1; break;}
					tmp += tree.tdsp(r.e, schedule[i], tmp,false);
					for (int k = i+1; k < w.S.size(); ++k){
						tmp += tree.tdsp(schedule[k-1], schedule[k],tmp,false);
						if (w.S[k] & 1){			  // satisfy all request end at k
							if(DDLFlagEndK(w, k, tmp) == -1){ddl_flag = -1;break;}  //wait change 2: capacity
						}
					}
					if(ddl_flag == -1){break;}
					if(tmp < delta){delta = tmp, optimanl_i = i, optimanl_j = j;}

				}else
				{
					tmp += tree.tdsp(r.s, schedule[i], tmp,false);
					//if(picked[i - 1] + r.com > w.cap){ddl_flag = -1; break;}
					for(int k = i+1; k < w.S.size(); ++k){

						if(k == j){
							tmp += tree.tdsp(schedule[k-1], r.e, tmp,false);
							if(tmp > r.ddl){ddl_flag = -1;break;} 
							tmp += tree.tdsp(r.e, schedule[k],tmp,false);							
						}
						else
						{
							tmp += tree.tdsp(schedule[k-1], schedule[k],tmp,false);
						}
						if (w.S[k] & 1){			  // satisfy all request end at k
							if(DDLFlagEndK(w, k, tmp) == -1){ddl_flag = -1;break;}
						}	
					}
					if(ddl_flag == -1){break;}
					if (tmp < delta){delta = tmp, optimanl_i = i, optimanl_j = j;}						
				}				
			}
			//if(ddl_flag == -1){break;}				
		}		
	}
}
*/