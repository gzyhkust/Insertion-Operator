# include "insertion.h"

#include "TDGT.h"
#include "io.h"

int nV, m, c, n;
double gridL, alpha;
int pos = 0;
Request* R = NULL;
Worker* W = NULL;
TDGT tree;

int numOfTDSP = 0;
int numofTDSPInitial = 0;

double EPS = 0.001, INF = 100000;

Grid* gr = NULL;
int* anchor = NULL;
int graph_len, graph_wid, grid_len, grid_wid;
int grid_num_per_row, grid_num_per_col, grid_sz;
double mnx = INF, mxx = -INF, mny = INF, mxy = -INF;
vector<bool> visitGrid;
int gridm = 0;

ofstream dump_result(dump_result_path);

void readInput() {
	
    cout << "begin worker" << endl;
	ifstream other;
	other.open(workers_path.c_str());
	if (!other.is_open()) {
		printf("%s does not exist\n", workers_path.c_str());
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

    cout << "begin tree" <<endl;
    ifstream readtree;
    readtree.open(TDGT_tree_path);
    read_TDGT(readtree, tree);
    cout << "end tree" <<endl;

    cout << "begin request" << endl;
	ifstream ifs;
	ifs.open(requests_path.c_str());
	if (!ifs.is_open()) {
		printf("%s does not exist\n", requests_path.c_str());
		exit(0);
	}
	ifs >> n;
	R = new Request[n];
	for (int i = 0; i < n; ++ i) {
		cout << i << endl;
		ifs >> R[i].tim >> R[i].s >> R[i].e >> R[i].com;
		R[i].len = tree.tdsp(R[i].s,R[i].e,R[i].tim,false) * 2;
		R[i].ddl = R[i].tim + R[i].len + ddl;
		R[i].pr = pr;
	}
	ifs.close();
	cout << "end request" << endl;	
	cout << m << " " << n << " " << c << " " << gridL << " " << alpha << endl;

}

void timeDependentInsertion(){
	// no grid index
	vector<int> car;
	for(int i = 0; i<m; i++){
		car.push_back(i);
	}

	while (pos<n)
	{
		for (int i = 0; i < m; ++i){
			updateDriver(i, R[pos].tim);
		}

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
		//ans += alpha * tmp;
		w.tim += tmp;
		w.pid = Pos(w.S[0]);

		if (w.S[0] & 1) {
			w.num -= R[w.S[0] >> 1].com;
			//mcnt ++;
		} else {
			w.num += R[w.S[0] >> 1].com;
		}

		//w.trajectory.push_back(tuple<int,int,double>(Pos(w.S[0]), w.S[0], w.reach[0]));
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
double DDLEndPos(int x){
    return R[x>>1].ddl; 
}

void assignTaxi(vector<int>& car) {
	double opt = INF, fit;
	int sz = car.size(), id = -1;
	int optimal_i, optimal_j, x, y;
	

	clock_t assignbegin = clock();
	if (R[pos].len < INF) {
		for (int i = 0; i < sz; ++ i) {
			fit = INF;

			dump_result<< "For request: " << pos << " assignTaxi: " <<  i <<"ddl: "<<R[pos].ddl<<"\n";
			dump_result << "taxi schedule: " << "size:"<<W[car[i]].S.size() <<"|  "<< W[car[i]].S << "\n";
			dump_result << "taxi reach: " << "size:"<<W[car[i]].reach .size() <<"|  "<< W[car[i]].reach  << "\n";
			try_insertion(W[car[i]], pos, fit, x, y);
			dump_result << "------------------------------Finish taxi------------------------------" <<"\n";

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
		dump_result << "insertion at taxi: "<< id << "at :" <<optimal_i<<"  "<<optimal_j<< "\n";
		//cout << "insertion at taxi: "<< id << "at :" <<optimal_i<<"  "<<optimal_j<< endl;
		cout << "Assigin request: " << pos << endl;
        insertion(W[id], pos, id, optimal_i, optimal_j);
	}

	clock_t afterassignbegin = clock();
	cout << "assighn \t time cost: " + to_string((afterassignbegin - assignbegin) / CLOCKS_PER_SEC) + "s" << endl;

	pos++;
	dump_result << "******************************************Finish request***************************************" << "\n"; 
}

/*
given the arrival time a i, get the arrival time at j 
*/
double function_arrival(Worker &w, int i, int j,double arr){
	if (i == -1)
	{
		return INF;
	}
	else{	
		if(i==j){
			if(arr <= w.ddl[i]){
				return arr;
			}
			else{
				return INF;
			}
		}
		else
		{
			if(arr <= w.ddl[i]){
				return w.FunctionTimesegments[i][j-i].dpt2arr(arr);
			}
			else{
				return INF;
			}
		}
	}
}

void try_insertion(Worker &w, int rid, double &delta, int &optimanl_i, int &optimanl_j) {
	Request& r = R[rid];
	double opt = INF;
	dump_result << "try_insertion:"<< "\n"; 

	if (w.S.empty()) {
		double tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
		numOfTDSP++; 
        tmp += tree.tdsp(r.s, r.e, tmp,false);
		numOfTDSP++;
		if (tmp < r.ddl + EPS && r.com <= w.cap) {
			delta = tmp, optimanl_i = 0, optimanl_j = 1;
		}
		dump_result << "try_insertion case 0, empty taxi: " << delta<< "\n";  
		return;
	}

	vector<int>& picked = w.picked;
	vector<int>& schedule = w.S;
	vector<double>& reach = w.reach;
	vector<double>& ddl = w.ddl;
	int delta_i = -1;

	double det_i = INF, det_j = INF;	// new arrival time at next position after insert r.s before i or j 
	for (int j = 0; j <= w.S.size(); j++)
	{
		if (j == 0)
		{
			double tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
			det_i = tmp;

			if(w.num + r.com <= w.cap){
				det_i += tree.tdsp(r.s, Pos(schedule[j]), det_i, false);
				delta_i = 0;
				if (det_i <= ddl[j])
				{
					tmp += tree.tdsp(r.s, r.e, tmp,false);
					tmp += tree.tdsp(r.e, Pos(schedule[j]), tmp, false);
					opt = function_arrival(w, j, schedule.size()-1, tmp);
					if (opt < delta)
					{
						delta = opt;
						optimanl_i = j;
						optimanl_j = j;
					}
				}
				else{delta_i = -1;det_i = INF;}
			}
			else{delta_i = -1;det_i = INF;}

			// if (w.num + r.com <= w.cap)
			// {
			// 	det_i = tmp + tree.tdsp(r.s, Pos(schedule[0]), tmp,false);			
			// 	delta_i = 0;
			// 	det_j = det_i;

			// 	tmp += tree.tdsp(r.s, r.e, tmp,false);
			// 	if(tmp <= r.ddl){
			// 		tmp += tree.tdsp(r.e, Pos(schedule[0]), tmp,false);
			// 		opt = function_arrival(w,j,schedule.size() - 1, tmp);
			// 		if(opt < delta){
			// 			delta = opt;
			// 			optimanl_i = j;
			// 			optimanl_j = j;
			// 		}			
			// 	}
			// }

		}else if (j == w.S.size())
		{

			//case 1
			double tmp = function_arrival(w, delta_i, j-1, det_i);
			if (tmp < INF)
			{
				tmp += tree.tdsp(Pos(schedule[j-1]), r.e, tmp,false);
				if(tmp <= r.ddl && tmp < delta){
					delta = tmp;
					optimanl_i = delta_i;
					optimanl_j = j;
				}
			}

			//case 2
			tmp = reach[j-1] + tree.tdsp(Pos(schedule[j-1]), r.s, reach[j-1], false);
			tmp += tree.tdsp(r.s, r.e, tmp, false);
			if(tmp <= r.ddl && tmp < delta){
				delta = tmp;
				optimanl_i = j;
				optimanl_j = j;
			}

			
			// //case 1 : insert r.s at delta_i, insert r.e at w.S.size()
			// if(delta_i != -1){
			// 	double tmp = function_arrival(w,delta_i,j-1, det_i);
			// 	if(tmp <= ddl[j-1]){
			// 		tmp += tree.tdsp(Pos(schedule[j-1]), r.e, tmp, false);
			// 	}
			// 	tmp += tree.tdsp(Pos(schedule[j-1]), r.e, tmp, false);
			// 	if(tmp <= r.ddl && tmp < delta){
			// 		delta = tmp;
			// 		optimanl_i = delta_i;
			// 		optimanl_j = j;
			// 	}
			// }
			// //case 2 : insert r.s, r.e at w.S.size()
			// if(picked[j-1] + r.com <= w.cap){
			// 	double tmp = reach[w.S.size() - 1] + tree.tdsp(Pos(schedule[w.S.size()-1]), r.s, reach[w.S.size() - 1] ,false);
			// 	tmp += tree.tdsp(r.s, r.e, tmp ,false);
			// 	if (tmp <= r.ddl && tmp < delta)
			// 	{
			// 		delta = tmp;
			// 		optimanl_i = j;
			// 		optimanl_j = j;
			// 	}
			// }
			
		}else{

			// try to update delta_i
			det_j  = reach[j-1] + tree.tdsp(Pos(schedule[j-1]),r.s, reach[j-1],false);
			double tmp = det_j;
			det_j += tree.tdsp(r.s, Pos(schedule[j]), det_j, false);
			if(picked[j-1]+r.com <= w.cap && det_j <= ddl[j] && det_j < function_arrival(w,delta_i,j, det_i)){
				// update i
				delta_i  = j; 
				det_i = det_j;

				tmp += tree.tdsp(r.s, r.e, tmp,false);
				if(tmp <= r.ddl){
					tmp += tree.tdsp(r.e, Pos(schedule[j]), tmp, false);
					tmp = function_arrival(w, j ,schedule.size()-1, tmp);
					if(tmp < delta){
						delta = tmp;
						optimanl_i = j;
						optimanl_j = j;
					}					
				}
			}
			else{
				tmp = function_arrival(w, delta_i, j-1, det_i);

				tmp += tree.tdsp(Pos(schedule[j-1]), r.e, tmp, false);
				if(tmp <= r.ddl){
					tmp += tree.tdsp(r.e, Pos(schedule[j]), tmp, false);
					tmp = function_arrival(w, j, schedule.size()-1, tmp);
					if(tmp < delta){
						delta = tmp;
						optimanl_i = delta_i;
						optimanl_j = j;
					}
				}
				// if(tmp <= ddl[j-1]){
				// 	tmp += tree.tdsp(Pos(schedule[j-1]), r.e, tmp, false);
				// 	if(tmp <= r.ddl){
				// 		tmp += tree.tdsp(r.e, Pos(schedule[j]), tmp, false);
				// 		tmp = function_arrival(w, j, schedule.size()-1, tmp);
				// 		if(tmp < delta){
				// 			delta = tmp;
				// 			optimanl_i = delta_i;
				// 			optimanl_j = j;
				// 		}
				// 	}
				// }
			}


			/////////////////////////////////////////////////////////////////////////////
			// double tmp = reach[j-1] + tree.tdsp(Pos(schedule[j-1]),r.s, reach[j-1],false);

			// det_j = tmp + tree.tdsp(r.s, Pos(schedule[j]), tmp, false);
			// if (det_j < function_arrival(w,delta_i,j, det_i) && picked[j-1]+r.com <= w.cap){
			// 	delta_i = j;
			// 	det_i = det_j;

			// 	tmp += tree.tdsp(r.s, r.e, tmp ,false);
			// 	tmp += tree.tdsp(r.e, Pos(schedule[j]), tmp, false);
			// 	opt = function_arrival(w, j , schedule.size()-1, tmp);
			// 	if(opt < delta){
			// 		delta = opt;
			// 		optimanl_i = j;
			// 		optimanl_j = j;
			// 	}
			// }
			// else{ // insert r.s at delta_i, insert r.e at j
			// 	tmp = function_arrival(w,delta_i, j-1, det_i);
			// 	tmp += tree.tdsp(Pos(schedule[j-1]), r.e, tmp, false);

			// 	if (tmp <= r.ddl)
			// 	{
			// 		tmp += tree.tdsp(r.e, Pos(schedule[j]), tmp, false);

			// 		opt = function_arrival(w, j ,schedule.size()-1, tmp);
			// 		if(opt < delta){
			// 			delta = opt;
			// 			optimanl_i = delta_i;
			// 			optimanl_j = j;
			// 		}
			// 	}
			// }
		}				
	}
}


void insertion(Worker &w, int rid, int wid, int optimal_i, int optimal_j){
	if (w.S.empty())
	{
		dump_result <<" insertion:" <<"\n";
		dump_result << "assign request:" << rid <<" "<< "into:" << optimal_i<<" "<<optimal_j << "\n";
		dump_result << "before insertion: "<<"\n";
		dump_result <<  w.S << "\n";
		dump_result <<  w.reach << "\n";
		dump_result << w.ddl << "\n";

		w.S.push_back(rid << 1);
		w.S.push_back(rid << 1 | 1);
		updateDriverArr(w);

		dump_result <<"after insertion: " <<"\n";
		dump_result <<  w.S << "\n";
		dump_result <<  w.reach << "\n";
		dump_result << w.ddl << "\n";

		return;
	}

	vector<int>& picked = w.picked;
	vector<double>& reach = w.reach;

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

    cout << "Assign request:" << rid <<" "<< R[rid].ddl << " to worker : " << wid << endl;
    cout << "before insertion: "<<endl;
    cout <<  w.S << endl;
    cout <<  w.reach << endl;

	dump_result <<" insertion:" <<"\n";
	dump_result << "assign request:" << rid <<" "<< "into:" << optimal_i<<" "<<optimal_j << "\n";
    dump_result << "before insertion: "<<"\n";
    dump_result <<  w.S << "\n";
    dump_result <<  w.reach << "\n";
	dump_result << w.ddl << "\n";

	w.S.clear();
	w.S = ret;
	updateDriverArr(w);

    cout << "after insertion: " << optimal_i << " " << optimal_j <<endl;
    cout <<  w.S << endl;
    cout <<  w.reach << endl;
    cout << "------------------------------------" << endl;
	dump_result <<"after insertion: " <<"\n";
    dump_result <<  w.S << "\n";
    dump_result <<  w.reach << "\n";
	dump_result << w.ddl << "\n";
}



void updateDriverArr(Worker& w){
	double tim = w.tim;
	vector<double>& reach = w.reach;
	vector<int>& schedule = w.S;

	vector<vector<PLF>>& FunctionTimesegments = w.FunctionTimesegments;
	//vector<vector<int>>& functionRange = w.functionRange;
	for (auto it = FunctionTimesegments.begin(); it != FunctionTimesegments.end(); it++)
	{
		it->clear();
	}
	
	FunctionTimesegments.clear();
	//functionRange.clear();
	reach.clear();

	// o(n) update reach
	for (int x = 0; x < w.S.size(); ++x)
	{
		double tim_x = 0;
		if(x == 0){
			tim_x = tree.tdsp(w.pid, Pos(schedule[x]), tim, false);
			numofTDSPInitial++;
			tim += tim_x;

			reach.push_back(tim);			
		}
		else{

			tim_x = tree.tdsp(Pos(schedule[x-1]), Pos(schedule[x]), tim, false);
			numofTDSPInitial++;
			tim += tim_x;
			reach.push_back(tim);	
		}
	}


	// update F function
	for (int x = 0; x < w.S.size(); ++x)
	{

		vector<PLF> FunctionS; // from x to y
		vector<int> range;
		for (int y = x; y < w.S.size(); ++y)
		{	
			PLF PLFy;
			if(y == x){
				FunctionS.push_back(PLFy);
				range.push_back(0);
			}
			else if (y == x+1)
			{
				tree.PLCst(Pos(schedule[y-1]), Pos(schedule[y]), reach[y-1], TMAX, PLFy);
				numofTDSPInitial++;
				if(y == x+1){

					int x = 0;
					FunctionS.push_back(PLFy);
					range.push_back(x);
				}				
			}else
			{
				
				tree.PLCst(Pos(schedule[y-1]), Pos(schedule[y]), reach[y-1], TMAX, PLFy);


				int x = 0;
				PLF PLFxy;
				PLFy.compound(FunctionS[FunctionS.size()-1],PLFxy,Pos(schedule[y-1]));
				range.push_back(x);
				FunctionS.push_back(PLFxy);
			}
		}
		
		FunctionTimesegments.push_back(FunctionS);		
	}

	//o(n) uodate ddl
	// cout << "vector<double>& ddl = w.ddl;" <<endl;
	vector<double>& ddl = w.ddl;
	ddl.clear();
	// cout << "after vector<double>& ddl = w.ddl;" <<endl;

	for (int k = w.S.size() - 1; k >= 0; --k)
	{
		if(k == w.S.size() - 1){
			ddl.insert(ddl.begin(), DDLEndPos(w.S[k]));

		}
		else
		{	
			// update ddl of k related to k+1
			double delat = INF;
			double opt = INF;
			if(w.S[k] & 1){
				opt = DDLEndPos(w.S[k]);
			}

			auto s1 = FunctionTimesegments[k][1].f->begin();
			
			if (s1->t + s1->w >= ddl[0])
			{
				delat = ddl[0] - s1->w;
			}
			else{
				while (s1->t + s1->w < ddl[0])
				{
					s1++;
				}
				auto s2 = s1 - 1;

				if (s1 == FunctionTimesegments[k][1].f->end())
				{
					delat = ddl[0] - s1->w;
				}
				else{

					double k = (s1->w - s2->w) / (s1->t - s2->t);
					delat = ( ddl[0] + k*s2->t -s2->w )/(1+k);
				}
			}

			if(delat < opt){
				opt = delat;
			} 
			ddl.insert(ddl.begin(), opt);
			
		}
	}

	// o(n) update picked
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
	cout << "after vector<int>& picked = w.picked;" <<endl;

}

void finishTaxi(int i) {
	Worker& w = W[i];
	
	for (int i=0; i<w.S.size(); ++i) {
		double tmp = w.reach[i] - w.tim;
		///ans += alpha * tmp;
		w.tim += tmp;
		if (w.S[i] & 1) {
			//mcnt ++;
            continue;
		}
	}
}

void freeMemory() {
	delete[] R;
	delete[] W;
}

void recordTrajectory(){

	ofstream dump_trajectory("/homes/zgongae/n/trajectory.txt");
	for (int i = 0; i < m; i++)
	{
		Worker& w = W[i];
		cout << "Taxi: " << i <<endl;
		dump_trajectory << "Taxi: " << i << "  #"<< w.trajectory.size()/2 <<"\n";
		for(auto& tuple: w.trajectory){
			cout << "Arrival node: "<<get<0>(tuple) << "  at time: " << get<2>(tuple) << endl;
			dump_trajectory << "Arrival node: "<<get<0>(tuple)<<" "<<get<1>(tuple) << "  at time: " << get<2>(tuple) <<"\n";
		}
		cout << " --------------------------------------------- " <<endl;
	}
	
	cout << numOfTDSP << endl;
	cout << numofTDSPInitial << endl;

	dump_trajectory.close();
}

