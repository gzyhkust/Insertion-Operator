# include "insertionN2.h"

#include "TDGT.h"
#include "io.h"

//int nV, m, c, n,numOfVertices;
int nV,numOfVertices;
double gridL, alpha;
int pos = 0;
Request* R = NULL;
Worker* W = NULL;
vector<vertex> vertices;
TDGT tree;

//ofstream loginsertion(insertion_log_n2);
//////////////////// result metric
int carnum = 0, cartime = 0;
double assigntime = 0; int assignnum = 0;
double insertiontime = 0; int insertionnum = 0;
int numOfTDSP = 0;
int numofTDSPInitial = 0;
int serverPassenger = 0;

int numOfQuery[10] = {0};
void addNumofQuery()
{
	numOfQuery[9] ++;

	for (int i = 9; i > 0; i--)
	{
		if (numOfQuery[i] >= 1000000)
		{
			numOfQuery[i-1] ++;
			numOfQuery[i] = numOfQuery[i] - 1000000;
		}
	}
}

clock_t assignbegin, assignend, insertionbegin, insertionend;
////////////////////

double EPS = 0.001, INF = 100000;

Grid* gr = NULL;
int* anchor = NULL;
int graph_len, graph_wid, grid_len, grid_wid;
int grid_num_per_row, grid_num_per_col, grid_sz;
double mnx = INF, mxx = -INF, mny = INF, mxy = -INF;
vector<bool> visitGrid;
int gridm = 0;

Position lefttop, leftdown; 
Position righrtop, rightdown;

#define pi 3.1415926535897932384626433832795
#define EARTH_RADIUS 6378.137
double rad(double d)
{
    return d * pi /180.0;
}
double RealDistance(double lat1,double lng1,double lat2,double lng2)
{
	double a;
	double b;
	double radLat1 = rad(lat1);
	double radLat2 = rad(lat2);
	a = radLat1 - radLat2;
	b = rad(lng1) - rad(lng2);
	double s = 2 * asin(sqrt(pow(sin(a/2),2) + cos(radLat1)*cos(radLat2)*pow(sin(b/2),2)));
	s = s * EARTH_RADIUS;
	s = s * 1000;
	return s;
}
double RealDistance(Position p1, Position p2)
{
	double a;
	double b;
	double radLat1 = rad(p1.y);
	double radLat2 = rad(p2.y);
	a = radLat1 - radLat2;
	b = rad(p1.x) - rad(p2.x);
	double s = 2 * asin(sqrt(pow(sin(a/2),2) + cos(radLat1)*cos(radLat2)*pow(sin(b/2),2)));
	s = s * EARTH_RADIUS;
	s = s * 1000;
	return s;
}

void readInput() {
	
    cout << "begin worker" << endl;
	ifstream other;
	other.open(workers_path.c_str());
	if (!other.is_open()) {
		printf("%s does not exist\n", workers_path.c_str());
		exit(0);
	}
	W = new Worker[m];
	for (int i = 0; i < m; ++ i) {
		other >> W[i].pid >> W[i].cap;
		W[i].cap = c;
		W[i].num = 0;		
	}
	other.close();
    cout << "end worker" << endl;

    cout << "begin tree" <<endl;
    ifstream readtree;
    readtree.open(TDGT_tree_path);
    read_TDGT(readtree, tree);
    cout << "end tree" <<endl;

	cout << "begin coordinates" << endl;
	ifstream vertexFile;
    vertexFile.open(vertexes_path.c_str());
	if (!vertexFile.is_open()) {
		printf("%s does not exist\n", vertexes_path.c_str());
		exit(0);
	}
    vertexFile >> numOfVertices;
	vertices.resize(numOfVertices);
    for(int i=0; i<numOfVertices; ++i) {
        vertexFile >> vertices[i].y >> vertices[i].x;
		vertices[i].id = i;
    }
    vertexFile.close();	
	cout << "end coordinates" << endl;

    cout << "begin request" << endl;
	ifstream ifs;
	ifs.open(requests_path.c_str());
	if (!ifs.is_open()) {
		printf("%s does not exist\n", requests_path.c_str());
		exit(0);
	}

	R = new Request[n];
	cout << "R = new Request[n];" << endl;
	for (int i = 0; i < n; ++ i) {
		ifs >> R[i].tim >> R[i].s >> R[i].e >> R[i].com >>R[i].ddl;
		R[i].ddl += ddlconstant;
		//ifs >> R[i].tim >> R[i].s >> R[i].e >> R[i].len >>R[i].com;
		//R[i].ddl = R[i].tim + R[i].len + deltasecond;  // need modify
	}
	ifs.close();
	cout << "end request" << endl;
	gridL = 5000;
	alpha = 1;	
	cout << m << " " << n << " " << c << " " << gridL << " " << alpha << endl;

}

//////////////////////////////////////////////////////////
///Index related
inline int getGridID(int id) {
	
	Position p1, p2;
	p1.x = lefttop.x, p1.y = vertices[id].y;
	double rdis = RealDistance(p1,lefttop);
	int rid = rdis/grid_wid;
	p2.x = vertices[id].x, p2.y = lefttop.y;
	double cdis = RealDistance(p2,lefttop);
	int cid = cdis/grid_len;
	return grid_num_per_row * rid + cid;
}

inline int getGridID(double x, double y) {
	int rid = x / grid_wid;
	int cid = y / grid_len;
	return grid_num_per_row * rid + cid;
}

inline double euclid(int a, int b) {
	//vector<vertex>& vertices = sssp->vertices;
	return sqrt((vertices[a].x-vertices[b].x)*(vertices[a].x-vertices[b].x) + (vertices[a].y-vertices[b].y)*(vertices[a].y-vertices[b].y)) / MAX_speed;
}

inline double euclid(Position &a, Position &b) {
	//return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y)) / MAX_speed;
	return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y)) / MAX_speed;
}

inline double euclidDist(Position &a, Position &b) {
	return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y));
}

inline double euclidDist(int a, int b) {
	return sqrt((vertices[a].x-vertices[b].x)*(vertices[a].x-vertices[b].x) + (vertices[a].y-vertices[b].y)*(vertices[a].y-vertices[b].y));
}

double centerDist(int gid, int pid) {
	Position g, pp;
	int rid = gid / grid_num_per_row;
	int cid = gid % grid_num_per_row;
	double lx = rid * grid_wid, rx = min(1. * graph_wid, lx + grid_wid);
	double ly = cid * grid_len, ry = min(1. * graph_len, ly + grid_len);
	g.x = (lx + rx) * 0.5;
	g.y = (ly + ry) * 0.5;
	pp.x = vertices[pid].x;
	pp.y = vertices[pid].y;
	return euclidDist(g, pp);
}

void insertTaxi(int gid, int tid) {
	gr[gid].taxi.push_back(tid);
	int vid = gr[gid].taxi.size() - 1;
	W[tid].gid = gid;
	W[tid].vid = vid;
}

void deleteTaxi(int gid, int tid) {
	if (gid>=0 && tid>=0 && W[tid].vid>=0 && gr[gid].taxi.size()>W[tid].vid) {
		int tid_ = *gr[gid].taxi.rbegin();
		gr[gid].taxi[W[tid].vid] = tid_;
		W[tid_].vid = W[tid].vid;
		gr[gid].taxi.pop_back();
		W[tid].vid = W[tid].gid = -1;
	}
}

void initGrid() {
	for (int i = 0; i < numOfVertices; ++ i) {
		mnx = min(mnx, vertices[i].x);
		mxx = max(mxx, vertices[i].x);
		mny = min(mny, vertices[i].y);
		mxy = max(mxy, vertices[i].y);
	}
	lefttop.x = mnx, lefttop.y = mxy;
	leftdown.x = mnx,leftdown.y = mny;
	righrtop.x = mxx, righrtop.y = mxy;
	rightdown.x = mxx, rightdown.y = mny;

	// graph_len = RealDistance(mnx, mny, mxx, mny);
	// graph_wid = RealDistance(mnx, mxy, mnx, mny);
	graph_len = RealDistance(leftdown, rightdown);
	graph_wid = RealDistance(leftdown, lefttop);
	cout << graph_len <<" " << graph_wid << endl;
	gridL = grid_value;
	grid_len = gridL; 
	grid_wid = gridL; 

	grid_num_per_row = (graph_len + grid_len) / grid_len;
	grid_num_per_col = (graph_wid + grid_wid) / grid_wid;
	grid_sz = grid_num_per_row * grid_num_per_col;
	visitGrid.resize(grid_sz);
	printf("grid size = %d\n", grid_sz);
	fflush(stdout);

	cout << mnx << " " << mxx <<" " << mny << " " << mxy << endl;
	cout << graph_len << " "<<grid_wid << endl;
	cout << gridL << endl;
	cout << grid_num_per_row<<" "<<grid_num_per_col << endl;
	
	gr = new Grid[grid_sz];
	anchor = new int[grid_sz];
	memset(anchor, -1, sizeof(int)*grid_sz);
	for(int i = 0; i < numOfVertices; ++ i) {
		int id = getGridID(i);
		if (anchor[id] == -1 || centerDist(id, i) < centerDist(id, anchor[id])) {
			anchor[id] = i;
		}
	}

	for (int i = 0; i < m; ++ i) {
		int id = W[i].pid;
		insertTaxi(getGridID(id), i);
	}
}

void updateGrid(int pid, int wid, int tag) {
	int gid = getGridID(pid);
	if (tag == 1) {
		insertTaxi(gid, wid);
	} else {
		deleteTaxi(gid, wid);
	}
}

vector<int> single_search(int s, double ddl) {
	
	//int gid0 = getGridID(vertices[s].x, vertices[s].y);
	int gid0 = getGridID(s);
	int x0 = gid0 / grid_num_per_row;
    int y0 = gid0 % grid_num_per_row;
	double radius = ddl * MAX_speed;
	int rx = ceil(radius / grid_wid), ry = ceil(radius / grid_len);
	vector<int> ret;
	
	fill(visitGrid.begin(), visitGrid.end(), false);
	for (int dx=-rx; dx<=rx; ++dx) {
		for (int dy=-ry; dy<=ry; ++dy) {
			int x = x0 + dx, y = y0 + dy;
			if (x>=0 && x<grid_num_per_col && y>=0 && y<grid_num_per_row) {
				int gid = x*grid_num_per_row + y;
				if (!visitGrid[gid]) {
					visitGrid[gid] = true;
					for (int i = 0; i<gr[gid].taxi.size(); ++ i) {
						ret.push_back(gr[gid].taxi[i]);
					}
				}
			}
		}
	}
	gridm = ret.size();
	return ret;
}
//////////////////////////////////////////////////////////
void timeDependentInsertion(){
	while (pos<n)
	{
		for (int i = 0; i < m; ++i){
			updateDriver(i, R[pos].tim);
		}

		vector<int> car;
		for(int i = 0; i<m; i++){
			car.push_back(i);
		}
		//cout <<" car num: " << car.size() << endl;
		carnum += car.size();
		cartime ++ ;
		assignnum ++;
		// cout << "Assign request: " << pos <<"  with car size: "<< car.size() << endl;
		// //cout << W[0].S << endl;
		// for (int i = 0; i < W[0].S.size(); i++)
		// {
		// 	cout << W[0].S[i] << ": " << W[0].reach[i] << "; ";
		// }
		// cout << endl;
		
		// cout << R[pos].s << " "<< R[pos].e << " " << endl;

		// loginsertion << "Assign request: " << pos <<"  with car size: "<< car.size() << "\n";
		// loginsertion << R[pos].s << " "<< R[pos].e << " " << R[pos].ddl <<  "\n";

		assignbegin = clock();
		assignTaxi(car);
		assignend = clock();
		cout << "**********************************************************************" <<endl;
		//loginsertion << "**********************************************************************" << "\n";
		assigntime += double(assignend - assignbegin) / CLOCKS_PER_SEC;
	}
	for (int i = 0; i < m; ++i){
		finishTaxi(i);
	}
}


void updateDriver(int i, double t){
	Worker& w = W[i];

	while (w.S.size() > 0 && w.tim < t)
	{
		updateGrid(w.pid, i, -1);	
		double tmp = w.reach[0] - w.tim; // == dist(w.pid, Pos(w.S[0]))
		w.tim += tmp;
		w.pid = Pos(w.S[0]);
		updateGrid(w.pid, i, 1);

		if (w.S[0] & 1) {
			w.num -= R[w.S[0] >> 1].com;
		} else {
			w.num += R[w.S[0] >> 1].com;
		}
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

			if(try_insertion_euclidDist(W[car[i]], pos) == 1)
			{
				insertionnum ++;
				insertionbegin = clock();
				try_insertion(W[car[i]], pos, fit, x, y);
				insertionend = clock();
				insertiontime += double(insertionend - insertionbegin) / CLOCKS_PER_SEC;

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
	}
	
	if (id > -1) {
        insertion(W[id], pos, id, optimal_i, optimal_j);
		serverPassenger++;
	}
	pos++;
}

int try_insertion_euclidDist(Worker &w, int rid)
{
	return 1;
}

/*
given the arrival time ai i, get the arrival time at j 
*/
double function_arrival(Worker &w, int i, int j,double arr){

	if(i==j || Pos(w.S[i]) == Pos(w.S[j])){
		return arr;
	}
	else
	{
		return w.FunctionTimesegments[i][j-i].dpt2arr(arr);
	}
}

void try_insertion(Worker &w, int rid, double &delta, int &optimanl_i, int &optimanl_j) {
	Request& r = R[rid];
	double opt = INF;

	if (w.S.empty()) {
		double tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
		numOfTDSP++;
		addNumofQuery(); 
        tmp += tree.tdsp(r.s, r.e, tmp,false);
		numOfTDSP++;
		addNumofQuery();
		if (tmp < r.ddl + EPS && r.com <= w.cap) {
			delta = tmp, optimanl_i = 0, optimanl_j = 1;
		} 
		return;
	}

	vector<int>& picked = w.picked;
	vector<int>& schedule = w.S;
	vector<double>& reach = w.reach;
	vector<double>& ddl = w.ddl;
	double tmp = 0;
	int pickedi = w.num, pickedj = pickedi; 
	for (int i = 0; i <= w.S.size(); ++ i){
		if (i>0)
		{
			pickedi = picked[i-1];
		}
		if (pickedi + r.com > w.cap)
		{
			continue;
		}
		
		
		for (int j = i; j<= w.S.size(); ++j){

			if(j>i)
			{
				pickedj = picked[j-1];
				if (pickedj + r.com > w.cap)
				{
					break;
				}
				
			}

			// loginsertion << "*****************************************" << "\n";
			// loginsertion << "try_insertion: " << i << " " << j << "\n";
			// cout << "*****************************************" << endl;
			// cout << "try_insertion: " << i << " " << j << "   " << schedule[i] << " " << schedule[j] << endl;

			if (i == 0 && j ==0){ //case 0:  i==j==0

				tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
				//loginsertion << "r.s  tree.tdsp(" << w.pid << "," << r.s <<"): " << tmp << "\n";
				numOfTDSP++;
				addNumofQuery();
				tmp += tree.tdsp(r.s, r.e, tmp,false);
				//loginsertion << "r.e  tree.tdsp(" << r.s << "," << r.e <<"): " << tmp << "\n";
				numOfTDSP++;
				addNumofQuery();
				if (tmp<=r.ddl){ //at least w can pick r
					tmp += tree.tdsp(r.e, Pos(schedule[0]), tmp,false);
					//loginsertion << "schedule[i]  tree.tdsp(" << r.e << "," << schedule[0] <<"): " << tmp << "\n";
					numOfTDSP++;
					addNumofQuery();
					if(tmp <= ddl[0]){
						opt = function_arrival(w,i,schedule.size() - 1, tmp);
						//loginsertion << "obj  function_arrival(" << schedule[i] << "," << schedule[schedule.size()-1] <<"): " << tmp << " "<< opt << "\n";
						if(opt < delta){
							//loginsertion <<"update delta"<< opt << " " << i <<"  "<<j<< "\n";
							delta = opt; optimanl_i = i; optimanl_j = j;				
						} 
					}
				}
			
			}else if (i == w.S.size()){ //case 1:  i==j==w.S.size()
				tmp = reach[w.reach.size() - 1] + tree.tdsp(Pos(schedule[w.S.size() - 1]), r.s, reach[w.reach.size() - 1], false);
				//loginsertion << "r.s  tree.tdsp(" << schedule[w.S.size() - 1] << "," << r.s <<"): " << tmp << "\n";
				numOfTDSP++;
				addNumofQuery();
				tmp += tree.tdsp(r.s, r.e, tmp,false);
				//loginsertion << "r.e  tree.tdsp(" << r.s << "," << r.e <<"): " << tmp << "\n";
				numOfTDSP++;
				addNumofQuery();
				if(tmp<delta && tmp <= r.ddl && r.com <= w.cap)
				{
					//loginsertion <<"update delta"<< tmp << " " << i <<"  "<<j<< "\n";
					delta = tmp, optimanl_i = i, optimanl_j = j;  				
				}				
				
				// tmp = reach[reach.size() - 1] + tree.tdsp(Pos(schedule[schedule.size() - 1]), r.s, reach[reach.size() - 1], false);
				// numOfTDSP++;
				// tmp += tree.tdsp(r.s, r.e, tmp,false);
				// numOfTDSP++;
				// //if (tmp > r.ddl || r.com > w.cap)
				// if (tmp > r.ddl)
				// {
				// 	loginsertion <<"ddl&cap false S.size()"<< tmp << " " << r.ddl <<"  "<<picked[i-1]<<" "<<r.com << "\n";
				// 	continue;
				// }
				// opt = tmp;
				// if(opt<delta){
				// 	loginsertion <<"update delta"<< opt << " " << i <<"  "<<j<< "\n";
				// 	delta = opt, optimanl_i = i, optimanl_j = j;
				// }
			}
			else{

				if(i == 0){ // arrival time r.s after inserting it
					tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
					//loginsertion << "r.s  tree.tdsp(" << w.pid << "," << r.s <<"): " << tmp << "\n";
					numOfTDSP++;
					addNumofQuery();
				}else
				{

					tmp = reach[i-1] + tree.tdsp(Pos(schedule[i-1]), r.s, reach[i-1],false);
					//loginsertion << "r.s  tree.tdsp(" << schedule[i-1] << "," << r.s <<"): " << tmp << "\n";
					numOfTDSP++;
					addNumofQuery();
				}

				if(i == j){ //case 2:  i==j
					tmp += tree.tdsp(r.s, r.e, tmp,false);
					//loginsertion << "r.e  tree.tdsp(" << r.s << "," << r.e <<"): " << tmp << "\n";
					numOfTDSP++;
					addNumofQuery();				
					if(tmp > r.ddl)
					{
						//loginsertion <<"ddl false(r.e)"<< tmp << " " << r.ddl << "\n";
						continue;
					}
					tmp += tree.tdsp(r.e, Pos(schedule[i]), tmp,false);
					//loginsertion << "schedule[i]  tree.tdsp(" << r.e << "," << schedule[i] <<"): " << tmp << "\n";
					numOfTDSP++;
					addNumofQuery();
					if (tmp <= ddl[i])
					{
						opt = function_arrival(w, i,schedule.size() - 1, tmp);
						//loginsertion << "obj  function_arrival(" << schedule[i] << "," << schedule[schedule.size()-1] <<"): " << tmp << " "<< opt << "\n";
						if(opt < delta){
							//loginsertion <<"update delta"<< opt << " " << i <<"  "<<j<< "\n";
							delta = opt, optimanl_i = i, optimanl_j = j; 							
						}
					}
					// else
					// {
					// 	loginsertion <<"ddl false (i):  "<< schedule[i] << " " << tmp << " " << ddl[i] <<  "\n";
					// }
				

				}else //case 3,4: i and j in general ;i in general, j==w.s.size();
				{

					tmp += tree.tdsp(r.s, Pos(schedule[i]), tmp,false);
					//loginsertion << "schedule[i]  tree.tdsp(" << r.s << "," << schedule[i] <<"): " << tmp << "\n";
					numOfTDSP++;
					addNumofQuery();

					if (tmp <= ddl[i]){					
						
						// auto s1 = w.FunctionTimesegments[i][j-1-i].f->begin();
						// while (s1 != w.FunctionTimesegments[i][j-1-i].f->end())
						// {
						// 	cout << s1->t << " " << s1->w <<"; ";
						// 	s1 ++;
						// }
						// cout <<endl;						
						double arrj_1 = function_arrival(w,i,j-1, tmp);
						//loginsertion << "schedule[j-1]  function_arrival(" << schedule[i] << "," << schedule[j-1] <<"): " << tmp << " " << arrj_1 << "\n";

						tmp  = arrj_1 + tree.tdsp(Pos(schedule[j-1]), r.e, arrj_1,false);
						//loginsertion << "r.e  tree.tdsp(" << schedule[j-1] << "," << r.e <<"): " << tmp << "\n";
						numOfTDSP++;
						addNumofQuery();
						if(tmp <= r.ddl){
							if (j == w.S.size())
							{
								if(tmp < r.ddl){
									opt = tmp;
								}else
								{
									//loginsertion <<"ddl false (r.e):"<<tmp<< " " << r.ddl << "\n";							
									continue;
								}
							
							}
							else
							{
								tmp += tree.tdsp(r.e, Pos(schedule[j]), tmp,false); 
								//loginsertion << "schedule[j]  tree.tdsp(" << r.e << "," << schedule[j] <<"): " << tmp << "\n";														
								numOfTDSP++;
								addNumofQuery();

								if (tmp <= ddl[j])
								{
									// auto s1 = w.FunctionTimesegments[j][schedule.size()-1-j].f->begin();
									// while (s1 != w.FunctionTimesegments[j][schedule.size()-1-j].f->end())
									// {
									// 	cout << s1->t << " " << s1->w <<"; ";
									// 	s1 ++;
									// }
									// cout <<endl;
									opt = function_arrival(w, j,schedule.size()-1, tmp);
									//cout << "obj  function_arrival(" << schedule[j] << "," << schedule[schedule.size()-1] <<"): " << opt << endl;							

								}else
								{
									//loginsertion <<"ddl false(j):"<<j<<" "<<tmp<< " " << ddl[j] << "\n";
									continue;
								}
														
							}
							
							if (opt < delta)
							{
								//loginsertion <<"update delta"<< opt << " " << i <<"  "<<j<< "\n";
								delta = opt, optimanl_i = i, optimanl_j = j;
							}
							
						}
						else
						{
							//loginsertion <<"ddl false (r.e):"<<tmp<< " " << r.ddl << "\n";
							continue;
						}
											
					}else
					{
						//loginsertion <<"ddl false (i):" <<schedule[i]<<" "<<tmp<< " " << ddl[i] << "\n";
						continue;
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

		cout << "insert at empty car: "<< wid << endl;
		cout << "before insertion: "<<endl;
		cout <<  w.S << endl;
		cout <<  w.reach << endl;

		// loginsertion << "insert at empty car: "<< wid << "\n";
		// loginsertion << "before insertion: " << "\n";
		// for (int k = 0; k < w.S.size(); k++)
		// {
		// 	loginsertion << w.S[k] << " ";
		// }
		// loginsertion << "\n" ;	
		// for (int k = 0; k < w.reach.size(); k++)
		// {
		// 	loginsertion << w.reach[k] << " ";
		// }
		// loginsertion << "\n" ;

		updateDriverArr(w);

		cout << "after insertion: " << optimal_i << " " << optimal_j <<endl;
		cout <<  w.S << endl;
		cout <<  w.reach << endl;
		cout << "------------------------------------" << endl;

		// loginsertion << "after insertion: " << optimal_i << " " << optimal_j << "\n";
		// for (int k = 0; k < w.S.size(); k++)
		// {
		// 	loginsertion << w.S[k] << " ";
		// }
		// loginsertion << "\n" ;
		// for (int k = 0; k < w.reach.size(); k++)
		// {
		// 	loginsertion << w.reach[k] << " ";
		// }
		// loginsertion << "\n" ;
		// loginsertion << "------------------------------------" << "\n";

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

	// loginsertion << "Assign request:" << rid <<" "<< R[rid].ddl << " to worker : " << wid << "\n";
	// loginsertion << "before insertion: " << "\n";
	// for (int k = 0; k < w.S.size(); k++)
	// {
	// 	loginsertion << w.S[k] << " ";
	// }
	// loginsertion << "\n" ;
	// for (int k = 0; k < w.S.size(); k++)
	// {
	// 	loginsertion << w.picked[k] << " ";
	// }
	// loginsertion << "\n" ;	
	// for (int k = 0; k < w.reach.size(); k++)
	// {
	// 	loginsertion << w.reach[k] << " ";
	// }
	// loginsertion << "\n" ;
	// for (int k = 0; k < w.reach.size(); k++)
	// {
	// 	loginsertion << w.ddl[k] << " ";
	// }
	// loginsertion << "\n" ;

	w.S.clear();
	w.S = ret;
	updateDriverArr(w);

    cout << "after insertion: " << optimal_i << " " << optimal_j <<endl;
    cout <<  w.S << endl;
    cout <<  w.reach << endl;
    cout << "------------------------------------" << endl;

	// loginsertion << "after insertion: " << optimal_i << " " << optimal_j << "\n";
	// for (int k = 0; k < w.S.size(); k++)
	// {
	// 	loginsertion << w.S[k] << " ";
	// }
	// loginsertion << "\n" ;
	// for (int k = 0; k < w.S.size(); k++)
	// {
	// 	loginsertion << w.picked[k] << " ";
	// }
	// loginsertion << "\n" ;	
	// for (int k = 0; k < w.reach.size(); k++)
	// {
	// 	loginsertion << w.reach[k] << " ";
	// }
	// loginsertion << "\n" ;
	// for (int k = 0; k < w.reach.size(); k++)
	// {
	// 	loginsertion << w.ddl[k] << " ";
	// }
	// loginsertion << "\n" ;
	// loginsertion << "------------------------------------" << "\n";
	
}


void updateDriverArr(Worker& w){
	double tim = w.tim;
	vector<double>& reach = w.reach;
	vector<int>& schedule = w.S;

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

	for (int x = 0; x < w.S.size(); ++x)
	{

		vector<PLF> FunctionS; // from x to y
		for (int y = x; y < w.S.size(); ++y)
		{	
			PLF PLFy;
			if(y == x){
				FunctionS.push_back(PLFy);
			}
			else if (y == x+1)
			{
				//tree.PLCst(Pos(schedule[y-1]), Pos(schedule[y]), reach[y-1], TMAX, PLFy);
				tree.PLCst(Pos(schedule[y-1]), Pos(schedule[y]), 0, TMAX, PLFy);
				if(y == x+1){
					FunctionS.push_back(PLFy);					
				}				
			}else
			{
				
				//tree.PLCst(Pos(schedule[y-1]), Pos(schedule[y]), reach[y-1], TMAX, PLFy);
				tree.PLCst(Pos(schedule[y-1]), Pos(schedule[y]), 0, TMAX, PLFy);
				
				PLF PLFxy;
				PLFy.compound(FunctionS[FunctionS.size()-1],PLFxy,Pos(schedule[y-1]));
				FunctionS.push_back(PLFxy);
			}
		}
		FunctionTimesegments.push_back(FunctionS);		
	}

	vector<double>& ddl = w.ddl;
	ddl.clear();

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
			
			auto function = FunctionTimesegments[k][1].f;
			auto s1 = function->begin();
			while (s1->t + s1->w < ddl[0] && s1 < function->end())
			{
				s1 ++;
			}
			auto s2 = s1 - 1;
			if (s1 == function->end())
			{
				delat = ddl[0] - s2->w;
			}
			else
			{
				double slop = (s1->w - s2->w) / (s1->t - s2->t);
				delat = (ddl[0] - s2->w + s2->t * slop)/(slop+1);
			}

			if(delat < opt){
				opt = delat;
			} 
			ddl.insert(ddl.begin(), opt);			
		}
	}

	// // //according reach&ddl to segment the F function
	// for (int x = 0; x < schedule.size(); x++)
	// {
	// 	for (int y = x+1; y < schedule.size(); y++)
	// 	{
	// 		double arr,d;
	// 		arr = reach[x], d = ddl[y];	

	// 		while ((FunctionTimesegments[x][y-x].f->begin()+1) -> t < arr & FunctionTimesegments[x][y-x].f->size() >= 5)
	// 		{
	// 			FunctionTimesegments[x][y-x].f->erase(FunctionTimesegments[x][y-x].f->begin());
	// 		}
	// 		double t_, w_;
	// 		auto it1 = FunctionTimesegments[x][y-x].f->begin();
	// 		auto it2 = it1+1;
	// 		double slop = (it2->w - it1->w)/(it2->t - it1->t);
	// 		it1->w = slop*(arr - it1->t) + it1->w;
	// 		it1->t = arr;

	// 		int s = FunctionTimesegments[x][y-x].f->size() - 1;
	// 		while ((FunctionTimesegments[x][y-x].f->begin() + s-1) ->t > d & FunctionTimesegments[x][y-x].f->size() >= 5 )
	// 		{
	// 			FunctionTimesegments[x][y-x].f->erase(FunctionTimesegments[x][y-x].f->begin() + s);
	// 			s--;
	// 		}	
	// 	}
		
	// }

	// calculate memory size of F
	int funmemory = 0;
	for (int x = 0; x < schedule.size(); x++)
	{
		for (int y = x+1; y < schedule.size(); y++)
		{
			auto s1 = FunctionTimesegments[x][y-x].f->begin();
			while (s1 < FunctionTimesegments[x][y-x].f->end())
			{
				funmemory += sizeof(s1->t);
				funmemory += sizeof(s1->w);
				s1++;
			}			
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
		if (cc > w.maxnum)
		{
			w.maxnum = cc;
		}		
	}

	// update memory
	int memo = 0;
	memo += 5 * sizeof(w.pid);
	memo += sizeof(w.tim);
	memo += w.S.size() * sizeof(w.S[0]);
	memo += w.reach.size() * sizeof(w.reach[0]);
	memo += w.picked.size() * sizeof(w.picked[0]);
	memo += w.ddl.size() * sizeof(w.ddl[0]);
	memo += funmemory;
	if(memo > w.memory)
	{
		w.memory = memo;
	}

}

void finishTaxi(int i) {
	Worker& w = W[i];
	
	for (int i=0; i<w.S.size(); ++i) {

		w.trajectory.push_back(tuple<int,double>(w.S[i], w.reach[i]));

		double tmp = w.reach[i] - w.tim;
		w.tim += tmp;
		if (w.S[i] & 1) {
            continue;
		}
	}
}

void freeMemory() {
	delete[] R;
	delete[] W;
	delete[] gr;
}

void recordTrajectory(){

	ofstream dump_trajectory(trajectory_n2_path);
	double totalmemory = 0;
	for (int i = 0; i < m; i++)
	{
		Worker& w = W[i];
		dump_trajectory << "Taxi: " << i << "  #"<< w.maxnum <<"  M"<< w.memory <<"\n";
		totalmemory += w.memory ;
		for(auto& tuple: w.trajectory){
			//cout << "Arrival node: "<<get<0>(tuple) << "  at time: " << get<2>(tuple) << endl;
			dump_trajectory << "Arrival node: "<<get<0>(tuple)<< "  at time: " << get<1>(tuple) <<"\n";
		}
	}
	dump_trajectory.close();


	ofstream dump_result(dump_resultn2_path);
	dump_result << "# of PathQuery: " << numOfTDSP << "\n";
	dump_result << "NumofQuery(): ";
	for (int i = 0; i < 10; i++)
	{
		dump_result << numOfQuery[i] << " ";
	}
	dump_result << "\n";
	dump_result << "Total assgntime: " << assigntime << ", Assignment num: " << assignnum << "Avg assignment time(ms): " << double(1000 * (assigntime/assignnum) ) << "\n";
	dump_result << "Total insertiontime: " << insertiontime << ", Insertion num: " << insertionnum << "Avg insertion time(ms): " << double(1000 * (insertiontime/insertionnum) ) << "\n";
	dump_result << "Memory cost: " << totalmemory<< "   " << totalmemory/1024 <<"(kb)" << "\n";
	dump_result << "Car num: " << carnum << " " << cartime << "Avg car candidate: " << carnum/cartime << "\n";
	dump_result << serverPassenger << "\n";
	dump_result.close();
}

// void updateDriverArr(Worker& w){
// 	double tim = w.tim;
// 	vector<double>& reach = w.reach;
// 	vector<int>& schedule = w.S;

// 	// vector<vector<PLF>> FunctionTimesegments;
// 	// vector<PLF>  FunctionTimesegment;
// 	vector<vector<PLF>>& FunctionTimesegments = w.FunctionTimesegments;
// 	for (auto it = FunctionTimesegments.begin(); it != FunctionTimesegments.end(); it++)
// 	{
// 		it->clear();
// 	}
	
// 	FunctionTimesegments.clear();
// 	reach.clear();
// 	cout <<"all clear"<<endl;
// 	for (int x = 0; x < w.S.size(); ++x)
// 	{
// 		double tim_x = 0;
// 		if(x == 0){
// 			tim_x = tree.tdsp(w.pid, Pos(schedule[x]), tim, false);
// 			numofTDSPInitial++;
// 			tim += tim_x;

// 			reach.push_back(tim);			
// 		}
// 		else{

// 			tim_x = tree.tdsp(Pos(schedule[x-1]), Pos(schedule[x]), tim, false);
// 			numofTDSPInitial++;
// 			tim += tim_x;
// 			reach.push_back(tim);	
// 		}
// 	}
// 	cout <<"get reach"<<endl;

// 	for (int x = 0; x < w.S.size(); ++x)
// 	{

// 		vector<PLF> FunctionS; // from x to y
// 		for (int y = x; y < w.S.size(); ++y)
// 		{	
// 			PLF PLFy;
// 			if(y == x){
// 				FunctionS.push_back(PLFy);
// 			}
// 			else if (y == x+1)
// 			{
// 				tree.PLCst(Pos(schedule[y-1]), Pos(schedule[y]), reach[y-1], TMAX, PLFy);
// 				numofTDSPInitial++;
// 				if(y == x+1){
// 					FunctionS.push_back(PLFy);
// 				}				
// 			}else
// 			{
				
// 				tree.PLCst(Pos(schedule[y-1]), Pos(schedule[y]), reach[y-1], TMAX, PLFy);
				
// 				PLF PLFxy;
// 				PLFy.compound(FunctionS[FunctionS.size()-1],PLFxy,Pos(schedule[y-1]));
// 				FunctionS.push_back(PLFxy);
// 			}
// 		}
// 		FunctionTimesegments.push_back(FunctionS);		
// 	}

// 	vector<double>& ddl = w.ddl;
// 	ddl.clear();

// 	cout << "start ddl:-------------------------------------"<<endl;
// 	for (int k = w.S.size() - 1; k >= 0; --k)
// 	{
// 		cout << k <<" " << w.S[k] << endl;
// 		if(k == w.S.size() - 1){
// 			ddl.insert(ddl.begin(), DDLEndPos(w.S[k]));

// 		}
// 		else
// 		{	
// 			// update ddl of k related to k+1
// 			double delat = INF;
// 			double opt = INF;
// 			if(w.S[k] & 1){
// 				opt = DDLEndPos(w.S[k]);
// 			}
			
// 			cout << " ******************** " <<endl;
// 			cout << ddl[0] << endl;
// 			cout << FunctionTimesegments[k][1] <<endl;
// 			cout << FunctionTimesegments[k][1].f << endl;
// 			cout << " ******************** " <<endl;


// 			auto s1 = FunctionTimesegments[k][1].f->begin();
			
// 			if (s1->t + s1->w >= ddl[0])
// 			{
// 				delat = ddl[0] - s1->w;
// 			}
// 			else{
// 				while (s1->t + s1->w < ddl[0])
// 				{
// 					//cout << "while (s1->t + s1->w < ddl[0])" <<endl;
// 					cout << s1->t << endl;
// 					s1++;
// 				}
// 				cout << " end while (s1->t + s1->w < ddl[0])" <<endl;
// 				auto s2 = s1 - 1;
// 				cout << s2->t <<endl;
// 				cout << s1->t << endl;
// 				cout << "s1 == FunctionTimesegments[k][1].f->end()" <<endl;
// 				if (s1 == FunctionTimesegments[k][1].f->end())
// 				{
// 					delat = ddl[0] - s1->w;
// 				}
// 				else{

// 					double k = (s1->w - s2->w) / (s1->t - s2->t);
// 					delat = ( ddl[0] + k*s2->t -s2->w )/(1+k);
// 				}
// 			}

// 			if(delat < opt){
// 				opt = delat;
// 			} 
// 			ddl.insert(ddl.begin(), opt);
			
// 		}
// 	}

// 	cout << "start picked"<<endl;
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

// }


// void try_insertion(Worker &w, int rid, double &delta, int &optimanl_i, int &optimanl_j) {
// 	Request& r = R[rid];
// 	double opt = INF;

// 	if (w.S.empty()) {
// 		cout <<" if (w.S.empty()) { " << endl;
// 		double tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
// 		numOfTDSP++; 
//         tmp += tree.tdsp(r.s, r.e, tmp,false);
// 		numOfTDSP++;
// 		if (tmp < r.ddl + EPS && r.com <= w.cap) {
// 			delta = tmp, optimanl_i = 0, optimanl_j = 1;
// 		} 
// 		return;
// 	}

// 	vector<int>& picked = w.picked;
// 	vector<int>& schedule = w.S;
// 	vector<double>& reach = w.reach;
// 	vector<double>& ddl = w.ddl;
// 	double tmp = 0;
// 	for (int i = 0; i <= w.S.size(); ++ i){
// 		for (int j = i; j<= w.S.size(); ++j){

// 			if (i == 0 && j ==0){ //case 0:  i==j==0

// 				tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
// 				numOfTDSP++;
// 				tmp += tree.tdsp(r.s, r.e, tmp,false);
// 				numOfTDSP++;
// 				if (tmp<=r.ddl && w.num+r.com<=w.cap){ //at least w can pick r
// 					tmp += tree.tdsp(r.e, Pos(schedule[0]), tmp,false);
// 					numOfTDSP++;
// 					if(tmp <= ddl[0]){
// 						opt = function_arrival(w,i,schedule.size() - 1, tmp);
// 					}
// 					else
// 					{
// 						//////////////////////////////////////////////
// 						cout << "i == 0 && j ==0 one" << endl;
// 						cout << tmp << " " << ddl[0] << endl;
// 						//////////////////////////////////////////////
// 						continue;
// 					}
// 				}else
// 				{
// 					//////////////////////////////////////////////
// 					cout << "i == 0 && j ==0 two" << endl;
// 					cout << tmp << " " << r.ddl << endl;
// 					//////////////////////////////////////////////
// 					continue;
// 				}
// 				if(opt < delta){
// 					delta = opt; optimanl_i = i; optimanl_j = j;				
// 				} 
			
// 			}else if (i == w.S.size()){ //case 1:  i==j==w.S.size()
// 				tmp = reach[reach.size() - 1] + tree.tdsp(Pos(schedule[schedule.size() - 1]), r.s, reach[reach.size() - 1], false);
// 				numOfTDSP++;
// 				tmp += tree.tdsp(r.s, r.e, tmp,false);
// 				numOfTDSP++;
// 				if (tmp > r.ddl || picked[i-1]+r.com > w.cap)
// 				{
// 					//////////////////////////////////////////////
// 					cout << "i==j==w.S.size()" << endl;
// 					cout << tmp << " " << r.ddl << endl;
// 					//////////////////////////////////////////////
// 					continue;
// 				}
// 				opt = tmp;
// 				if(opt<delta){
// 					delta = opt, optimanl_i = i, optimanl_j = j;
// 				}
// 			}
// 			else{

// 				if(i == 0){ // arrival time r.s after inserting it 
// 					tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
// 					numOfTDSP++;
// 					if(w.num+r.com > w.cap){continue;} 
// 				}else
// 				{	
// 					tmp = reach[i-1] + tree.tdsp(Pos(schedule[i-1]), r.s, reach[i-1],false);
// 					numOfTDSP++;
// 					if(picked[i-1]+r.com > w.cap){continue;}
// 				}

// 				if(i == j){ //case 2:  i==j
// 					tmp += tree.tdsp(r.s, r.e, tmp,false);
// 					numOfTDSP++;
// 					if(tmp > r.ddl)
// 					{
// 						//////////////////////////////////////////////
// 						cout << "i==j one" << endl;
// 						cout << tmp << " " << r.ddl << endl;
// 						//////////////////////////////////////////////
// 						continue;
// 					}
// 					tmp += tree.tdsp(r.e, Pos(schedule[i]), tmp,false);
// 					numOfTDSP++;
// 					if (tmp <= ddl[i])
// 					{
// 						opt = function_arrival(w, i,schedule.size() - 1, tmp);
// 						if(opt < delta){
// 							delta = opt, optimanl_i = i, optimanl_j = j; 							
// 						}
// 					}else
// 					{
// 						//////////////////////////////////////////////
// 						cout << "i==j two" << endl;
// 						cout << tmp << " " << ddl[i] << endl;
// 						//////////////////////////////////////////////
// 						continue;
// 					}
				

// 				}else //case 3,4: i and j in general ;i in general, j==w.s.size();
// 				{
// 					tmp += tree.tdsp(r.s, Pos(schedule[i]), tmp,false);
// 					numOfTDSP++;

// 					if (tmp <= ddl[i]){
// 						double arrj_1 = function_arrival(w,i,j-1, tmp);
// 						tmp  = arrj_1 + tree.tdsp(Pos(schedule[j-1]), r.e, arrj_1,false);
// 						numOfTDSP++;
// 						if(tmp <= r.ddl){
// 							if (j == w.S.size())
// 							{
// 								if(tmp < r.ddl){
// 									opt = tmp;
// 								}else
// 								{
// 									//////////////////////////////////////////////
// 									cout << "i and j in general one " << endl;
// 									cout << tmp << " " << r.ddl<< endl;
// 									//////////////////////////////////////////////								
// 									continue;
// 								}
							
// 							}else{
// 								tmp += tree.tdsp(r.e, Pos(schedule[j]), tmp,false); 
// 								numOfTDSP++;
// 								if (tmp <= ddl[j])
// 								{
// 									opt = function_arrival(w, j,schedule.size()-1, tmp);		
// 								}else
// 								{
// 									//////////////////////////////////////////////
// 									cout << "i and j in general two " << endl;
// 									cout << tmp << " " << ddl[j]<< endl;
// 									//////////////////////////////////////////////	
// 									continue;
// 								}
														
// 							}
// 						}
// 						else
// 						{
// 							//////////////////////////////////////////////
// 							cout << "i and j in general three " << endl;
// 							cout << tmp << " " << r.ddl<< endl;
// 							//////////////////////////////////////////////	
// 							continue;
// 						}

// 						if(opt < delta){
// 							delta = opt, optimanl_i = i, optimanl_j = j; 							
// 						}	
											
// 					}else
// 					{
// 						//////////////////////////////////////////////
// 						cout << "i and j in general four " << endl;
// 						cout << tmp << " " <<ddl[i]<< endl;
// 						//////////////////////////////////////////////	
// 						continue;
// 					}	
// 				}				
// 			}			
// 		}		
// 	}
// }