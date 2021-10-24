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
	cout << other.is_open() <<endl;
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
		ifs >> R[i].tim >> R[i].s >> R[i].e >> R[i].com;
		R[i].len = tree.tdsp(R[i].s,R[i].e,R[i].tim,false);
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

		w.trajectory.push_back(tuple<int,int,double>(Pos(w.S[0]), w.S[0], w.reach[0]));
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
	
	if (R[pos].len < INF) {
		for (int i = 0; i < sz; ++ i) {
			fit = INF;

			dump_result<< "For request: " << pos << " assignTaxi:" <<  i <<"\n";
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
		// cout << "insertion" << endl;
		dump_result << "insertion at taxi: "<< id << "at :" <<optimal_i<<"  "<<optimal_j<< "\n";
        insertion(W[id], pos, id, optimal_i, optimal_j);
		// dump_result << "after inseriton taxi schedule: " << W[id].S << "\n";
		// dump_result << "after inseriton taxi reach: " << W[id].reach << "\n"; 
	}
	pos++;
	dump_result << "******************************************Finish request***************************************" << "\n"; 
}

/*
given the arrival time ai i, get the arrival time at j 
*/
double function_arrival(Worker &w, int i, int j,double arr){
	if(i==j){
		return arr;
	}
	else
	{
		// auto psed = w.FunctionTimesegments[i][j].dpt2seg(arr);
		// double weight = w.FunctionTimesegments[i][j].dpt2wgt(arr,psed);
		// return arr+weight;
		//auto psed = w.FunctionTimesegments[i][j-i].dpt2seg(arr);
		//double weight = w.FunctionTimesegments[i][j-i].dpt2wgt(arr,psed);
		//return arr+weight;
		return w.FunctionTimesegments[i][j-i].dpt2arr(arr);
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
	double tmp = 0;
	for (int i = 0; i <= w.S.size(); ++ i){
		for (int j = i; j<= w.S.size(); ++j){

			if (i == 0 && j ==0){ //case 0:  i==j==0

				tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
				numOfTDSP++; 
				tmp += tree.tdsp(r.s, r.e, tmp,false);
				numOfTDSP++;
				if (tmp<=r.ddl && w.num+r.com<=w.cap){ //at least w can pick r
					tmp += tree.tdsp(r.e, Pos(schedule[0]), tmp,false);
					numOfTDSP++;
					if(tmp < ddl[0]){
						opt = function_arrival(w,i,schedule.size() - 1, tmp);
					}
				}else
				{
					break;
				}
				if(opt < delta){
					delta = opt; optimanl_i = i; optimanl_j = j;
					dump_result << "try_insertion case 1, i == j ==0 : " << delta<< "\n"; 				
					} 
			
			}else if (i == w.S.size()){ //case 1:  i==j==w.S.size()
				tmp = reach[reach.size() - 1] + tree.tdsp(Pos(schedule[schedule.size() - 1]), r.s, reach[reach.size() - 1], false);
				numOfTDSP++;
				tmp += tree.tdsp(r.s, r.e, tmp,false);
				numOfTDSP++;
				if (tmp > r.ddl || picked[i-1]+r.com > w.cap){break;}
				opt = tmp;
				if(opt<delta){
					delta = opt, optimanl_i = i, optimanl_j = j;
					dump_result << "try_insertion case 2, i == j ==w.s.size() : " << delta<< "\n"; 
					}
			}
			else{

				if(i == 0){ // arrival time r.s after inserting it 
					tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
					numOfTDSP++;
					if(w.num+r.com > w.cap){break;} 
				}else
				{	
					//cout << "try_insertion :" << rid << " case1" <<endl; 
					tmp = reach[i-1] + tree.tdsp(Pos(schedule[i-1]), r.s, reach[i-1],false);
					numOfTDSP++;
					if(picked[i-1]+r.com > w.cap){break;}
				}

				if(i == j){ //case 2:  i==j
					tmp += tree.tdsp(r.s, r.e, tmp,false);
					numOfTDSP++;
					if(tmp > r.ddl){break;}
					// cout << "if(tmp > r.ddl){break;}" << endl;
					// cout << "try_insertion :" << rid << " case2" <<endl;
					// cout <<"tmp += tree.tdsp(r.e, Pos(schedule[i]), tmp,false);" <<endl;
					tmp += tree.tdsp(r.e, Pos(schedule[i]), tmp,false);
					numOfTDSP++;
					if (tmp <= ddl[i])
					{
						// cout << "opt = function_arrival(w, i,schedule.size() - 1, tmp);" << endl;
						// cout << i << " " << schedule.size() - 1 << endl;
						opt = function_arrival(w, i,schedule.size() - 1, tmp);
						if(opt < delta){
							delta = opt, optimanl_i = i, optimanl_j = j;
							dump_result << "try_insertion case 3, i == j : " << delta<< "\n"; 							
							}
					}
					


				}else //case 3,4: i and j in general ;i in general, j==w.s.size();
				{
					// cout << "try_insertion :" << rid << " case3" <<endl;
					tmp += tree.tdsp(r.s, Pos(schedule[i]), tmp,false);
					numOfTDSP++;

					if (tmp <= ddl[i]){
						// cout << "try_insertion case 3:" << endl;
						// cout << "double arrj_1 = function_arrival(w,i,j-1, tmp);" << endl;
						// cout << i << " " << j-1 <<endl;
						double arrj_1 = function_arrival(w,i,j-1, tmp);
						tmp  = arrj_1 + tree.tdsp(Pos(schedule[j-1]), r.e, arrj_1,false);
						numOfTDSP++;
						if (j == w.S.size())
						{
							if(tmp < r.ddl){
								opt = tmp;
							}else{break;}
						
						}else{
							tmp += tree.tdsp(r.e, Pos(schedule[j]), tmp,false); 
							numOfTDSP++;
							if (tmp < ddl[j])
							{
								// cout << "try_insertion case 3:" << endl;
								// cout << "opt = function_arrival(w, j,schedule.size()-1, tmp);" << endl;
								// cout << j << " " << schedule.size()-1 <<endl;
								opt = function_arrival(w, j,schedule.size()-1, tmp);
								//dump_result << "try_insertion case 5, j == w.s.size() : " << delta<< "\n"; 		
							}else{break;}
													
						}
						if(opt < delta){
							delta = opt, optimanl_i = i, optimanl_j = j;
							dump_result << "try_insertion case 4,5, general case : " << delta<< "\n"; 							
							}						
					}	
				}				
			}			
		}		
	}
}

void insertion(Worker &w, int rid, int wid, int optimal_i, int optimal_j){
	if (w.S.empty())
	{
		w.S.push_back(rid << 1);
		w.S.push_back(rid << 1 | 1);
		updateDriverArr(w);
		return;
	}

	vector<int>& picked = w.picked;
	// vector<int> ret = w.S;
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

	// vector<vector<PLF>> FunctionTimesegments;
	// vector<PLF>  FunctionTimesegment;
	vector<vector<PLF>>& FunctionTimesegments = w.FunctionTimesegments;
	for (auto it = FunctionTimesegments.begin(); it != FunctionTimesegments.end(); it++)
	{
		it->clear();
	}
	
	FunctionTimesegments.clear();
	reach.clear();

	for (int x = 0; x < w.S.size(); ++x)
	{
		double tim_x = 0;
		// cout << "updateDriverArr case1" <<endl;
		// if(x == 0 and w.pid == Pos(schedule[0])){
		// 	reach.push_back(tim);
		// }else if (x == 0)
		// {
		// 	// cout << "updateDriverArr case2" <<endl;
		// 	tim_x = tree.tdsp(w.pid, Pos(schedule[0]), tim, false);
		// 	numofTDSPInitial++;
		// 	tim += tim_x;
		// 	reach.push_back(tim);
		if(x == 0){
			tim_x = tree.tdsp(w.pid, Pos(schedule[x]), tim, false);
			numofTDSPInitial++;
			tim += tim_x;
			// cout << "reach.push_back(tim);" <<endl;
			reach.push_back(tim);			
		}
		else{
			// cout << "updateDriverArr case3" <<endl;
			// cout << schedule <<endl;
			tim_x = tree.tdsp(Pos(schedule[x-1]), Pos(schedule[x]), tim, false);
			numofTDSPInitial++;
			tim += tim_x;
			// cout << "reach.push_back(tim);" <<endl;
			reach.push_back(tim);
			// cout << "after reach.push_back(tim);" <<endl;			
		}
	}


	for (int x = 0; x < w.S.size(); ++x)
	{

		vector<PLF> FunctionS; // from x to y
		for (int y = x; y < w.S.size(); ++y)
		{	
			PLF PLFy;
			if(y == x){
				// cout << "FunctionS.push_back(PLFy);" <<endl;
				FunctionS.push_back(PLFy);
				// cout << "after FunctionS.push_back(PLFy);" <<endl;
			}
			else if (y == x+1)
			{
				double tim_y = tree.tdsp_PLF(Pos(schedule[y-1]),Pos(schedule[y]),reach[y-1],false,PLFy);
				numofTDSPInitial++;
				if(y == x+1){
					FunctionS.push_back(PLFy);
				}				
			}else
			{
				PLF PLFxy;
				//PLF PLFy;
				cout << reach[y-1] << endl;
				double tim_y = tree.tdsp_PLF(Pos(schedule[y-1]),Pos(schedule[y]),reach[y-1],false,PLFy);
				PLFy.compound(FunctionS[FunctionS.size()-1],PLFxy,Pos(schedule[y-1]));
				FunctionS.push_back(PLFxy);
				// cout << "updateDriverArr case5" <<endl;
				// cout << x << " " << y << endl;
				// cout << schedule << endl;
				// cout << PLFy.f->size() << endl;
				// cout << FunctionS.size() <<endl;
				// cout << (*(FunctionS.end()-1)).f->size() << endl;
				// cout << (*(FunctionS.end()-1)).f->size() << endl;
				// (*(FunctionS.end()-1)).compound(PLFy,PLFxy,Pos(w.S[y-1]));

				// FunctionS[FunctionS.size()-1].compound(PLFy,PLFxy,Pos(schedule[y-1]));
				// FunctionS.push_back(PLFxy);
			}
		}
		
		//cout << "FunctionTimesegments.push_back(FunctionS);" <<endl;
		FunctionTimesegments.push_back(FunctionS);
		//cout << "after FunctionTimesegments.push_back(FunctionS);" <<endl;		
	}
	//cout << "w.FunctionTimesegments = FunctionTimesegments;" <<endl;
	//vector<vector<PLF>> WFunctionTimesegment = w.FunctionTimesegments;
	// w.FunctionTimesegments.clear();
	// w.FunctionTimesegments.resize(w.S.size());
	///w.FunctionTimesegments = FunctionTimesegments;

	// cout << "after w.FunctionTimesegments = FunctionTimesegments;" <<endl;


	
	// cout << "vector<double>& ddl = w.ddl;" <<endl;
	vector<double>& ddl = w.ddl;
	ddl.clear();
	// cout << "after vector<double>& ddl = w.ddl;" <<endl;

	for (int k = w.S.size() - 1; k >= 0; --k)
	{
		if(k == w.S.size() - 1){
			ddl.insert(ddl.begin(), DDLEndPos(w.S[k]));
			// cout << ddl.size()<< endl;
			// cout << k << " " << w.S[k] << endl;	
		}
		else
		{	
			// update ddl of k related to k+1
			double delat = INF;
			double opt = INF;
			if(w.S[k] & 1){
				opt = DDLEndPos(w.S[k]);
			}
			// cout << "w.FunctionTimesegments.size()" << endl;
			// cout << w.FunctionTimesegments.size() << endl;
			// cout << w.FunctionTimesegments[k].size() <<endl;
			// cout << k << " " << k+1 <<endl;
			// cout << ddl << endl;

			//auto plc = w.FunctionTimesegments[k][k+1];
			//delat = plc.arr2dpt(ddl[0], plc.f->begin());
			

			//delat = w.FunctionTimesegments[k][1].arr2dpt(ddl[k+1], w.FunctionTimesegments[k][1].f->begin());
			//delat = w.FunctionTimesegments[k][1].arr2dpt(ddl[k+1], w.FunctionTimesegments[k][1].dpt2seg(ddl[ddl.size()-1]));
			
			// auto segment = w.FunctionTimesegments[k][1];
			// for(auto it = segment.f->begin(); it!=segment.f->end()-1; it++){
			// 	if(it->t < ddl[ddl.size()-1] and segment.dpt2arr(it->t)>= ddl[ddl.size()-1]){
			// 		delat = segment.arr2dpt(ddl[ddl.size()-1],it);
			// 		break;
			// 	}
			// }
			auto segment = w.FunctionTimesegments[k][1];
			auto it = segment.f->begin();
			while (it!=segment.f->end()-1 and it->t < ddl[ddl.size()-1])
			{
				if(segment.dpt2arr(it->t) >= ddl[ddl.size()-1]){
					delat = segment.arr2dpt(ddl[ddl.size()-1],it-1);
					//delat = segment.arr2dpt(ddl[ddl.size()-1],it);
					break;
				}
				it++;
			}
			
						
			if(delat < opt){
				opt = delat;
			} 
			ddl.insert(ddl.begin(), opt);
		}
	}

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

	for (int i = 0; i < m; i++)
	{
		Worker& w = W[i];
		cout << "Taxi: " << i <<endl;
		dump_result << "Taxi: " << i <<"\n";
		for(auto& tuple: w.trajectory){
			cout << "Arrival node: "<<get<0>(tuple) << "  at time: " << get<2>(tuple) << endl;
			dump_result << "Arrival node: "<<get<0>(tuple)<<" "<<get<1>(tuple) << "  at time: " << get<2>(tuple) <<"\n";
		}
		cout << " --------------------------------------------- " <<endl;
	}
	
	cout << numOfTDSP << endl;
	cout << numofTDSPInitial << endl;

	dump_result.close();
}


// void test(int s, int e, double t) {
// 	double r1 = 0;
// 	r1 = tree.tdsp(s, e, t, true);
// 	cout << r1 << endl;
// }

	// // ------------------------------- update segments
	// // -------------------------------
	// vector<vector<PLF>> FunctionTimesegments;
	// for(int k = 0; k< w.S.size(); ++k){
	// 	double tim_k = 0;
		
	// 	//vector<tuple<double, double>> segments;
	// 	PLF PLFk;

	// 	if(k ==0){
	// 		//tim_k = tree.tdsp_segments(w.pid, Pos(w.S[k]), tim,false,segments);
	// 		tim_k = tree.tdsp_PLF(w.pid, Pos(w.S[k]), tim,false, PLFk);
	// 	}
	// 	else
	// 	{
	// 		//tim_k = tree.tdsp_segments(Pos(w.S[k-1]), Pos(w.S[k]), tim,false,segments);
	// 		tim_k = tree.tdsp_PLF(Pos(w.S[k-1]), Pos(w.S[k]), tim,false, PLFk);
	// 	}

	// 	tim += tim_k;
	// 	reach.push_back(tim);

	// 	//updateDriverSegment(w, segments, reach, k)	
	// }
	// ///////////////////////////////////////////////////////////

// void updateDriverSegment(Worker& w, vector<tuple<double, double>>& segments, vector<double>& reach, int k){

// 	if(k == 0){
// 		vector<tuple<double, double>> s;
// 		vector<vector<double>> b;

// 		vector< tuple<double, double>> ts2;
// 		ts2.resize(1, 100);
// 		for(auto it = segments.begin(); it!=segments.end(); it++){
// 			if (get<1>(*it) > reach[k])
// 			{
// 				// s.push_back(make_tuple(get<0>(*it),get<1>(*it)));
// 				// b.push_back(get<0>(*it));

// 				//ts2[0].push_back(make_tuple(get<0>(*it),get<1>(*it)));
// 				ts2.push_back(make_tuple(get<0>(*it),get<1>(*it)));
// 			}
			
// 		} 
// 		// w.timesegments[0] = s;
// 		// w.breakpoints[0][0] = b; 

// 		w.timesegments2[0] = ts2;
// 	}

// 	else
// 	{
// 		vector<tuple<double, double>> s;
// 		vector<vector<double>> b;

// 		vector< tuple<double, double>> ts2;

// 		// vector<tuple <double, double>> ts_ = ts2[k-1];
// 		// sort(ts_.begin(), ts_.end());
// 		for (auto it = segments.begin(); it!=segments.end(); it++)
// 		{
// 			if (get<1>(*it) > reach[k])
// 			{
// 				ts2.push_back(make_tuple(get<0>(*it),get<1>(*it)));
// 			}
// 		}
// 		w.timesegments2[k] = ts2;
		

// 		// vector<tuple<double, double>> s = new vector<tuple<double, double>>[k*k];

// 		// for (int x = 0; x <= k; x++)
// 		// {
// 		// 	// x to k
// 		// 	s_ = w.timesegments[x*(k-1) + (k-1)];
// 		// 	b = w.breakpoints[x][k-1];
// 		// 	//update s_ based on segments; generate w.timesegments[x*(k) + (k)]
// 		// 	for(auto it = segments.begin(); it!=segments.end(); it++){
// 		// 		// if (get<1>(*it) > reach[x])
// 		// 		// {
// 		// 		// 	s.push_back(make_tuple(get<0>(*it),get<1>(*it)));
// 		// 		// }
// 		// 	}				
// 		// }

// 		// w.timesegments = s;
// 	}
	
// }

	// ------------------------------- update Functiontimesegments
	// vector<vector <tuple<double, double>>> Functiontimesegments; 
	// for(int k = 0; k< w.S.size(); ++k){
	// 	segments =  w.timesegments2[k];
		
	// 	if(k == 0){
	// 		Functiontimesegments[0][0] = segments;
	// 	}
	// 	else
	// 	{	
	// 		vector <tuple<double, double>> pre_segments = w.Functiontimesegments[][]
	// 		for(auto it = segments.begin(); it!=segments.end(); it++){
				
	// 		}			
	// 	}
		
	// } 
	// -------------------------------


		
	// if( i== j){
	// 	auto pseg = w.FunctionTimesegments[i][w.S.size()-1].dpt2seg(arr);
	// 	double weight = w.FunctionTimesegments[i][w.S.size()-1].dpt2wgt(arr,pseg);
	// 	return arr + weight;		
	// }
	// else
	// {
	// 	auto psegj = w.FunctionTimesegments[i][w.S.size()-1].dpt2seg(arr);
	// }
	
	// return 0;

		// ------------------------------- update ddl
	//vector<double> ddl = w.ddl;
	// vector<double> ddl;
	// ddl.clear();
	// for (int k = w.S.size() - 1; k <= 0; --k){
	// 	if(k == w.S.size() - 1){
	// 		ddl[k] = DDLEndPos(w.S[k]);			
	// 	}
	// 	else
	// 	{	
	// 		// update ddl of k related to k+1
	// 		double delat = INF;
	// 		if(w.S[k] & 1){
	// 			delat = DDLEndPos(w.S[k]);
	// 		}
	// 		//auto segments =  w.timesegments2[k][k+1];
	// 		auto segments =  w.FunctionTimesegments[k][k+1];
	// 		//for(auto it = segments.begin(); it!=segments.end(); it++){
	// 		for(auto it = segments.end(); it>=segments.begin(); it--){
	// 			if(get<0>(*it) >= reach[k] and get<0>(*it) < ddl[k+1]){
	// 				if (ddl[k+1] - get<1>(*it) > get<0>(*it))
	// 				{
	// 					if(ddl[k+1] - get<1>(*it) < delat){
	// 						delat = ddl[k+1] - get<1>(*it);
	// 					}
	// 				}
					
	// 			}
	// 		}
	// 		ddl[k] = delat;
	// 	}
		
	// }
	// w.ddl = dll;


		// for (int x = 0; x < w.S.size(); ++x)
	// {
	// 	vector<PLF> FunctionS;
	// 	for(int y = x; y < w.S.size(); ++y)
	// 	{	
	// 		double tim_k = 0;
	// 		PLF PLFx;
	// 		if(y == x){
				
	// 			if(y ==0){
	// 				tim_k = tree.tdsp_PLF(w.pid, Pos(w.S[y]), tim,false, PLFx);
	// 				//FunctionTimesegments[0][0] = PLFx;
	// 				FunctionS.push_back(PLFx);
	// 			}
	// 			else
	// 			{
	// 				tim_k = tree.tdsp_PLF(Pos(w.S[y-1]), Pos(w.S[y]), tim,false, PLFx);
	// 			}
	// 			tim += tim_k;
	// 			reach.push_back(tim);				
	// 		}
	// 		else
	// 		{
	// 			PLF PLFxy;
	// 			tim_k = tree.tdsp_PLF(Pos(w.S[y-1]), Pos(w.S[y]), tim,false, PLFx);
	// 			// FunctionTimesegments[x][y-1].compound(PLFx,PLFxy,Pos(w.S[y-1]));
	// 			// FunctionTimesegments[x][y] = PLFxy;
	// 			// PLF PLFxy_1 = FunctionTimesegments[x][y-1];

	// 			if(y == x){
	// 				PLFxy_1 = 
	// 			}else{

	// 			}


	// 			PLFxy_1.compound(PLFx,PLFxy,Pos(w.S[y-1]));
	// 			FunctionS.push_back(PLFxy);
	// 		}
			
	// 	}
	// 	FunctionTimesegments.push_back(FunctionS);
	// }

			// for(int y = x; y < w.S.size(); ++y)
		// {	
		// 	double tim_k = 0;
		// 	PLF PLFx;
		// 	if(y == x){
				
		// 		if(y ==0){
		// 			tim_k = tree.tdsp_PLF(w.pid, Pos(w.S[y]), tim,false, PLFx);
		// 			FunctionS.push_back(PLFx);
		// 		}
		// 		else
		// 		{
		// 			tim_k = tree.tdsp_PLF(Pos(w.S[y-1]), Pos(w.S[y]), tim,false, PLFx);
		// 		}
		// 		tim += tim_k;
		// 		reach.push_back(tim);				
		// 	}
		// 	else
		// 	{
		// 		PLF PLFxy;
		// 		tim_k = tree.tdsp_PLF(Pos(w.S[y-1]), Pos(w.S[y]), tim,false, PLFx);
		// 		// FunctionTimesegments[x][y-1].compound(PLFx,PLFxy,Pos(w.S[y-1]));
		// 		// FunctionTimesegments[x][y] = PLFxy;
		// 		// PLF PLFxy_1 = FunctionTimesegments[x][y-1];

		// 		if(y == x){
		// 			PLFxy_1 = 
		// 		}else{

		// 		}


		// 		PLFxy_1.compound(PLFx,PLFxy,Pos(w.S[y-1]));
		// 		FunctionS.push_back(PLFxy);
		// 	}
			
		// }
		// FunctionTimesegments.push_back(FunctionS);

			// /////////////////////////////////////////////////////////////////////
	// vector<vector<PLF>> FunctionTimesegments;

	// reach.clear();
	
	// for (int x = 0; x < w.S.size(); ++x)
	// {
	// 	vector<PLF> FunctionS;

	// 	PLF PLFx;
	// 	double tim_k = 0;
	// 	if(x==0){
	// 		tim_k = tree.tdsp_PLF(w.pid, Pos(w.S[x]), tim,false, PLFx);
	// 		FunctionS.push_back(PLFx);
	// 		for (auto &e: *(PLFx.f)   )
	// 		{
	// 			cout << e << endl;
	// 		}
	// 	}
	// 	else
	// 	{
	// 		tim_k = tree.tdsp_PLF(Pos(w.S[x-1]), Pos(w.S[x]), tim,false, PLFx);
	// 		FunctionS.push_back(PLFx);
	// 	}
	// 	//FunctionS.push_back(PLFx);

	// 	for (auto &e: *(PLFx.f)   )
	// 	{
	// 		cout << e << endl;
	// 	}
	// 	cout << "  ----------------------------  " << endl;
	// 	tim += tim_k;
	// 	reach.push_back(tim);

	// 	for (auto &e: *(FunctionS.begin())->f   )
	// 	{
	// 		cout << e << endl;
	// 	}

	// 	for (int y = x+1; y < w.S.size(); ++y)
	// 	{
	// 		PLF PLFxy;
	// 		PLF PLFy;
	// 		tim_k = tree.tdsp_PLF(Pos(w.S[y-1]), Pos(w.S[y]), reach[y-1],false, PLFy);

	// 		///check point
	// 		// cout << "check poin: " <<  endl;
	// 		// for(auto &e:*PLFy.f){
	// 		// 	cout << e <<endl;
	// 		// }
	// 		// cout << " ---------- " << endl;
	// 		// for(auto &e:*FunctionS[y-1].f){
	// 		// 	cout << e <<endl;
	// 		// }			
	// 		///

	// 		cout << FunctionS.size() <<endl;
	// 		FunctionS.push_back(PLFy);
	// 		cout << FunctionS.size() <<endl;
	// 		cout << (*(FunctionS.begin())).f << endl;

	// 		for (auto &e: *PLFy.f)
	// 		{
	// 			cout << e << endl;
	// 		}
	// 		cout << " +++++++++++++++++++++++++++++++++ " <<endl;
	// 		for (auto &e: *(FunctionS.begin()+1)->f   )
	// 		{
	// 			cout << e << endl;
	// 		}			
			
	// 		(*(FunctionS.begin())).compound(PLFy,PLFxy,Pos(w.S[y-1]));
	// 		//FunctionS[y-1].compound(PLFy,PLFxy,Pos(w.S[y-1]));	

	// 		FunctionS.push_back(PLFxy);		
	// 	}
	// 	FunctionTimesegments.push_back(FunctionS);
		
	// }
	// vector<vector<PLF>> WFunctionTimesegment = w.FunctionTimesegments;
	// w.FunctionTimesegments.clear();
	// w.FunctionTimesegments.resize(w.S.size());
	// w.FunctionTimesegments = FunctionTimesegments;

	// /// check FunctionTimesegments
	// ///
