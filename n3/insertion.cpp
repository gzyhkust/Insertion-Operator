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
        cout << i << endl;
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

		// vector<int> car;
		// car = single_search(R[pos].s, R[pos].ddl - R[pos].tim);
		// cout << "assignTaxi:" << pos << endl;
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
		w.trajectory.push_back(tuple<int,int,double>(Pos(w.S[0]),w.S[0], w.reach[0]));
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

			//cout << "try insertion" << endl;
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
		// cout << "insertion" << endl;
		dump_result << "insertion at taxi: "<< id << "at :" <<optimal_i<<"  "<<optimal_j<< "\n";
        insertion(W[id], pos, id, optimal_i, optimal_j);
		// dump_result << "after inseriton taxi schedule: " << W[id].S << "\n";
		// dump_result << "after inseriton taxi reach: " << W[id].reach << "\n"; 
	}
	pos++;
	dump_result << "******************************************Finish request***************************************" << "\n"; 

}

void try_insertion(Worker &w, int rid, double &delta, int &optimanl_i, int &optimanl_j) {
	Request& r = R[rid];
	double opt = INF;
	dump_result << "try_insertion:"<< "\n"; 
	
	delta = INF;
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
	double tmp = 0;
	for (int i = 0; i <= w.S.size(); ++ i){
		for (int j = i; j<= w.S.size(); ++j){

			int false_flag = 0;

			if (i == 0 && j ==0){ //case 1: i==j==0

				tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
				numOfTDSP++; 
				tmp += tree.tdsp(r.s, r.e, tmp,false);
				numOfTDSP++;
				if (tmp<=r.ddl && w.num+r.com<=w.cap){ //at least w can pick r
					tmp += tree.tdsp(r.e, Pos(schedule[0]), tmp,false);

					if (w.S[0] & 1){			  // satisfy all request end at k
						if(tmp > DDLEndPos(w.S[0])){
							false_flag = 1;
						}
					}

					numOfTDSP++;
					for(int k = 1; k < w.S.size(); ++k){
						tmp += tree.tdsp(Pos(schedule[k-1]), Pos(schedule[k]), tmp, false);
						numOfTDSP++;
						if (w.S[k] & 1){			  // satisfy all request end at k
							if(tmp > DDLEndPos(w.S[k])){
								false_flag = 1;
                                break;
                            }
						}
					}
				}else
				{
					false_flag = 1;
				}
				if(tmp<delta && false_flag == 0){
					delta = tmp, optimanl_i = i, optimanl_j = j;
					dump_result << "try_insertion case 1, i == j ==0 : " << delta<< "\n"; 
				}
			
			}else if (i == w.S.size()){ //case 2: i==j==w.S.size()
				tmp = reach[w.reach.size() - 1] + tree.tdsp(Pos(schedule[w.S.size() - 1]), r.s, reach[w.reach.size() - 1], false);
				numOfTDSP++;
				tmp += tree.tdsp(r.s, r.e, tmp,false);
				numOfTDSP++;
				if (tmp > r.ddl || picked[i-1]+r.com > w.cap){false_flag = 1;}
				if(tmp<delta and false_flag == 0){
					delta = tmp, optimanl_i = i, optimanl_j = j;
					dump_result << "try_insertion case 2, i == j ==w.s.size() : " << delta<< "\n";  				
				}
			}
			else{

				if(i == 0){ // arrival time r.s after inserting it 
					tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
					numOfTDSP++;
					if(w.num+r.com > w.cap){false_flag = 1;} 
				}else
				{
					tmp = reach[i-1] + tree.tdsp(Pos(schedule[i-1]), r.s, reach[i-1],false);
					numOfTDSP++;
					if(picked[i-1]+r.com > w.cap){false_flag = 1;}
				}

				if(i == j){ //case 3: i==j
					tmp += tree.tdsp(r.s, r.e, tmp,false);
					numOfTDSP++;
					if(tmp > r.ddl){false_flag = 1;}
					tmp += tree.tdsp(r.e, Pos(schedule[i]), tmp,false);
                    numOfTDSP++;
					if(w.S[i] & 1){
                        if(tmp > DDLEndPos(w.S[i])){false_flag = 1;}
                    }
					for (int k = i+1; k < w.S.size(); ++k){
						tmp += tree.tdsp(Pos(schedule[k-1]), Pos(schedule[k]),tmp,false);
						numOfTDSP++;
						if (w.S[k] & 1){			
							if(tmp > DDLEndPos(w.S[k])){
								false_flag = 1;
								break;
							} 
						}
					}
					if(tmp < delta && false_flag == 0){
						delta = tmp, optimanl_i = i, optimanl_j = j;
						dump_result << "try_insertion case 3, i == j : " << delta<< "\n";
					}

				}
				
				else//case 4, 5: i and j in general ;i in general, j==w.s.size();
				{
					tmp += tree.tdsp(r.s, Pos(schedule[i]), tmp,false);
                    numOfTDSP++;
					if(w.S[i] & 1){
                        if(tmp > DDLEndPos(w.S[i])){false_flag = 1;}
                    }
					for(int k = i+1; k < w.S.size(); ++k){

						if(k == j){
							tmp += tree.tdsp(Pos(schedule[k-1]), r.e, tmp,false);
							numOfTDSP++;
							if(tmp > r.ddl){
								false_flag = 1;
								break;
							} 
							tmp += tree.tdsp(r.e, Pos(schedule[k]),tmp,false);							
							numOfTDSP++;
						}
						else
						{
							tmp += tree.tdsp(Pos(schedule[k-1]), Pos(schedule[k]),tmp,false);
							numOfTDSP++;
						}
						if (w.S[k] & 1){		
							if(tmp > DDLEndPos(w.S[k])){
								false_flag = 1;
								break;
							}
						}	
					}
					if(j == w.S.size()){//case 5: i in general, j==w.s.size();
						tmp += tree.tdsp(Pos(schedule[schedule.size()-1]), r.e, tmp, false);
						// dump_result << "try_insertion case 5, j ==w.s.size() :" << "\n"; 
						if(tmp>r.ddl){false_flag = 1;}
					}
					if (tmp < delta and false_flag == 0){
						delta = tmp, optimanl_i = i, optimanl_j = j;
						dump_result << "try_insertion case 4,5, general case : " << delta<< "\n"; 
						}						
				}				
			}			
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
		vector<double> ddl;
		for(int i = 0; i<w.S.size(); i++){
			if(w.S[i] & 1){
				ddl.push_back(DDLEndPos(w.S[i]));
			}
			else{ddl.push_back(0);}
		}
		dump_result << ddl << "\n";

		w.S.push_back(rid << 1);
		w.S.push_back(rid << 1 | 1);
		updateDriverArr(w);

		cout << "after insertion: " << optimal_i << " " << optimal_j <<endl;
		cout <<  w.S << endl;
		cout <<  w.reach << endl;
		cout << "------------------------------------" << endl;
		dump_result <<"after insertion: " <<"\n";
		dump_result <<  w.S << "\n";
		dump_result <<  w.reach << "\n";
		vector<double> ddl_;
		for(int i = 0; i<w.S.size(); i++){
			if(w.S[i] & 1){
				ddl_.push_back(DDLEndPos(w.S[i]));
			}
			else{ddl_.push_back(0);}
		}
		dump_result << ddl_ << "\n";

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
	vector<double> ddl;
	for(int i = 0; i<w.S.size(); i++){
		if(w.S[i] & 1){
			ddl.push_back(DDLEndPos(w.S[i]));
		}
		else{ddl.push_back(0);}
	}
	dump_result << ddl << "\n";

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
	vector<double> ddl_;
	for(int i = 0; i<w.S.size(); i++){
		if(w.S[i] & 1){
			ddl_.push_back(DDLEndPos(w.S[i]));
		}
		else{ddl_.push_back(0);}
	}
	dump_result << ddl_ << "\n";
}

void updateDriverArr(Worker& w){
	double tim = w.tim;
	vector<double>& reach = w.reach;

	reach.clear();
	for(int k = 0; k < w.S.size(); ++k){
		if (k == 0){
			tim += tree.tdsp(w.pid, Pos(w.S[k]), tim,false);
			numofTDSPInitial++;
		}
		else{
			tim += tree.tdsp(Pos(w.S[k-1]), Pos(w.S[k]), tim,false);
			numofTDSPInitial++;
		}
		reach.push_back(tim);
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

	//
	// cout << "schedule :" << w.S <<endl;
	// for(int i = 0; i<w.S.size(); ++i){
	// 	cout << Pos(w.S[i])<<endl;
	// }
	// cout << "reach :" << w.reach <<endl;
	// cout << "picked :" << w.picked <<endl;
	// //cout << w.ca
	// cout << "------------" << endl;
	//break;
	//		
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

	//ofstream dump_result(dump_result_path);
	for (int i = 0; i < m; i++)
	{
		Worker& w = W[i];
		cout << "Taxi: " << i <<endl;
		dump_result << "Taxi: " << i <<"\n";
		for(auto& tuple: w.trajectory){
			cout << "Arrival node: "<<get<0>(tuple)<<" "<<get<1>(tuple) << "  at time: " << get<2>(tuple) << endl;
			dump_result << "Arrival node: "<<get<0>(tuple)<<" "<<get<1>(tuple) << "  at time: " << get<2>(tuple) <<"\n";
		}
		cout << " --------------------------------------------- " <<endl;
	}
	
	cout << numOfTDSP << endl;
	cout << numofTDSPInitial << endl;

	dump_result.close();
}