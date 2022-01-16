# include "insertionN.h"

#include "TDGT.h"
#include "io.h"

int nV, m, c, n,numOfVertices;
double gridL, alpha;
int pos = 0;
Request* R = NULL;
Worker* W = NULL;
vector<vertex> vertices;
TDGT tree;

ofstream loginsertion(insertion_log_n);
//////////////////// result metric
int carnum = 0, cartime = 0;
double assigntime = 0; int assignnum = 0;
double insertiontime = 0; int insertionnum = 0;
int numOfTDSP = 0;
int numofTDSPInitial = 0;
int serverPassenger = 0;

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
	other >> m >> c >> gridL >> alpha;
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
	ifs >> n;
	cout  << n << endl;
	R = new Request[n];
	cout << "R = new Request[n];" << endl;
	for (int i = 0; i < n; ++ i) {
		ifs >> R[i].tim >> R[i].s >> R[i].e >> R[i].len >>R[i].com;
		R[i].ddl = R[i].tim + R[i].len * delta;
		//R[i].ddl = R[i].tim + R[i].len + deltasecond;
	}
	ifs.close();
	cout << "end request" << endl;	
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
		car = single_search(R[pos].s, R[pos].ddl-R[pos].tim);
		//cout <<" car num: " << car.size() << endl;
		carnum += car.size();
		cartime ++ ;
		assignnum ++;
		cout << "Assign request: " << pos <<"  with car size: "<< car.size() << endl;
		cout << W[0].S << endl;
		cout << R[pos].s << " "<< R[pos].e << " " << endl;
		
		loginsertion << "Assign request: " << pos <<"  with car size: "<< car.size() << "\n";
		loginsertion << R[pos].s << " "<< R[pos].e << " " << "\n";

		assignbegin = clock();
		assignTaxi(car);
		assignend = clock();
		cout << "**********************************************************************" <<endl;
		loginsertion << "**********************************************************************" << "\n";
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
	//cout << "assignTaxi(vector<int>& car)" <<endl;
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
				//cout << "try_insertion worker:"	<<  i <<endl;
				try_insertion(W[car[i]], pos, fit, x, y);
				//cout << " &&&&&&&&&&&&&&&&&&&&& " << endl;
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

/*
given the arrival time a i, get the arrival time at j 
*/
double function_arrival(Worker &w, int i, int j,double arr){
	if(arr >= INF)
	{
		return INF;
	}
	if (i <= -1)
	{
		return INF;
	}
	else{	
		if(i==j || Pos(w.S[i]) == Pos(w.S[j])){
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
				if (j == i+1)
				{
					return w.FunctionTimesegments[i].dpt2arr(arr);
				}
				else
				{
					return w.FunctionTimesegments[w.S.size()-1+i].dpt2arr(arr);
				}	
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

	if (w.S.empty()) {
		double tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
		numOfTDSP++; 
        tmp += tree.tdsp(r.s, r.e, tmp,false);
		numOfTDSP++;
		if (tmp < r.ddl + EPS && r.com <= w.cap) {
			delta = tmp, optimanl_i = 0, optimanl_j = 0;
		}
		return;
	}
	

	vector<int>& picked = w.picked;
	vector<int>& schedule = w.S;
	vector<double>& reach = w.reach;
	vector<double>& ddl = w.ddl;

	vector<double> scache, ecache, arricache;
	scache.resize(schedule.size()+1, 0), ecache.resize(schedule.size()+1, 0);
	arricache.resize(schedule.size(), 0);

	// ONE: for special case i == j
	for (int j = 0; j <= schedule.size(); j++)
	{
		loginsertion << "*****************************************" << "\n";
		loginsertion << "try_insertion i=j: " << j << "\n";		
		double tmp;
		if (j == 0)
		{
			if(w.num + r.com > w.cap)
			{
				loginsertion <<"cap false (r.s): "<< schedule[j] << "\n";
				continue;
			}
			tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim, false);
			scache[0] = tmp;
			loginsertion << "r.s  tree.tdsp(" << w.pid << "," << r.s <<"): " << w.tim << "\n";
			numOfTDSP++;
		}
		else
		{
			cout << "if(picked[j-1] + r.com > w.cap)" << endl;
			if(picked[j-1] + r.com > w.cap)
			{
				loginsertion <<"cap false (r.s): "<< schedule[j-1] << "\n";
				continue;
			}
			cout << "tmp = reach[j-1] + tree.tdsp(Pos(schedule[j-1]), r.s, reach[j-1], false);" << endl;
			tmp = reach[j-1] + tree.tdsp(Pos(schedule[j-1]), r.s, reach[j-1], false);
			scache[j] = tmp;
			loginsertion << "r.s  tree.tdsp(" << schedule[j-1] << "," << r.s <<"): " << tmp << "\n";
			numOfTDSP++;
		}
		
		cout << " tmp += tree.tdsp(r.s,r.e, tmp, false); " << endl;
		tmp += tree.tdsp(r.s,r.e, tmp, false);
		ecache[j] = tmp;
		loginsertion << "r.e  tree.tdsp(" << r.s << "," << r.e <<"): " << tmp << "\n";
		numOfTDSP++;

		
		if (j == schedule.size())
		{
			cout << " if (tmp <= r.ddl && tmp < delta) " << endl;
			if (tmp <= r.ddl && tmp < delta)
			{
				delta = tmp;
				optimanl_i = j;
				optimanl_j = j;
				cout << " loginsertion << " << endl;
				loginsertion <<"update delta"<< tmp << " " << optimanl_i <<"  "<<optimanl_j<< "\n";
			}	
		}
		else
		{
			if (tmp <= r.ddl)
			{
				tmp += tree.tdsp(r.e, Pos(schedule[j]), tmp, false);
				arricache[j] = tmp;
				loginsertion << "schedule[j]  tree.tdsp(" << r.e << "," << schedule[j] <<"): " << tmp << "\n";
				numOfTDSP++;
				if (tmp <= ddl[j])
				{
					tmp = function_arrival(w,j,schedule.size()-1, tmp);
					loginsertion << "obj  function_arrival(" << schedule[j] << "," << schedule[schedule.size()-1] <<"): " << tmp  << "\n";
					if (tmp < delta)
					{
						delta = tmp;
						optimanl_i = j;
						optimanl_j = j;
						loginsertion <<"update delta"<< tmp << " " << optimanl_i <<"  "<<optimanl_j<< "\n";
					}
				}				
			}		
		}			
	}
	

	cout << " vector<double> arrsj(w.S.size(), INF);  " << endl;
	vector<double> arrsj(w.S.size(), INF);
	vector<int> plc(w.S.size(), INF);

	// TW0: for j = 0;
	loginsertion << "*****************************************" << "\n";
	loginsertion << "try_insertion r.e: " << 0 << "\n";
	cout << " TW0: for j = 0; " << endl;
	if (w.num + r.com <= w.cap)
	{
		double tmp;
		//tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim, false);
		//numOfTDSP++;
		tmp = scache[0];
		loginsertion << "r.s  tree.tdsp(" << w.pid << "," << r.s <<"): " << tmp << "\n";
		if (tmp < r.ddl)
		{
			tmp += tree.tdsp(r.s,Pos(schedule[0]), tmp, false);
			numOfTDSP++;
			loginsertion << "schedule[i]  tree.tdsp(" << r.s << "," << schedule[0] <<"): " << tmp << "\n";
			if (tmp <= ddl[0])
			{
				plc[0] = 0;
				arrsj[0] = tmp;

				loginsertion << "update r.s location " << plc[0] << "\n";

			}
			
		}
		
	}

	// TW0: for j -> 1 to schedule.size()
	for (int j = 1; j <= schedule.size(); j++)
	{
		cout << " TW0: for j -> 1 to schedule.size():" << schedule[j] << endl;
		double tmp=INF, det_i = INF;
		loginsertion << "*****************************************" << "\n";
		loginsertion << "try_insertion r.e: " << j << "\n";

		if(picked[j-1] + r.com <= w.cap)
		{
			//step 1: check insertion r.e before j (ddl constrain & objective)			
			tmp = arrsj[j-1] + tree.tdsp(Pos(schedule[j-1]), r.e, arrsj[j-1],false);
			loginsertion << "r.e  tree.tdsp(" << schedule[j-1] << "," << r.e <<"): " << arrsj[j-1] << " " <<tmp << "\n";
			numOfTDSP++;
			if(tmp <= r.ddl)
			{
				if (j == schedule.size())
				{	
					if (tmp < delta)
					{
						delta = tmp;
						optimanl_i = plc[j-1];
						optimanl_j = j;
						loginsertion <<"update delta"<< tmp << " " << optimanl_i <<"  "<<optimanl_j<< "\n";						
					}
				
					continue;
				}
				
				cout << " before: tree.tdsp(r.e, Pos(schedule[j]), tmp, false)" << schedule[j] << endl;
				tmp += tree.tdsp(r.e, Pos(schedule[j]), tmp, false); ////////****************************
				cout << " after: tree.tdsp(r.e, Pos(schedule[j]), tmp, false)" << schedule[j] << endl;
				loginsertion << "schedule[j]  tree.tdsp(" << r.e << "," << schedule[j] <<"): " << tmp << "\n";
				numOfTDSP++;
				if (tmp <= ddl[j])
				{
					tmp = function_arrival(w,j,schedule.size()-1,tmp);
					loginsertion << "obj  function_arrival(" << schedule[j] << "," << schedule[schedule.size()-1] <<"): " << tmp << "\n";
					if (tmp < delta)
					{
						delta = tmp;
						optimanl_i = plc[j-1];
						optimanl_j = j;
						loginsertion <<"update delta"<< tmp << " " << optimanl_i <<"  "<<optimanl_j<< "\n";
					}
					
				}
				
			}


			if (j == schedule.size())
			{
				continue;
			}
			

			//step 2: update plc[j], insertion r.s before j is better
			plc[j] = plc[j-1];
			arrsj[j]= function_arrival(w,j-1,j,arrsj[j-1]);
			loginsertion << "insert r.s before j schedule[j] function_arrival(" << schedule[j-1] << "," << schedule[j] <<"): " << arrsj[j] << "\n";
			// tmp = reach[j-1] + tree.tdsp(Pos(schedule[j-1]), r.s, reach[j-1], false);
			// numOfTDSP++;
			tmp = scache[j];

			loginsertion << "potential new r.s  tree.tdsp(" << schedule[j-1] << "," << r.s <<"): " << tmp << "\n";
			if (tmp < r.ddl)
			{
				tmp += tree.tdsp(r.s, Pos(schedule[j]), tmp, false);
				loginsertion << "insert r.s at j schedule[j]  tree.tdsp(" << r.s << "," << schedule[j] <<"): " << tmp << "\n";
				numOfTDSP++;
				if (tmp < ddl[j])
				{
					if (tmp < arrsj[j])
					{
						plc[j] = j;
						arrsj[j] = tmp;		
						loginsertion <<"update i for r.s: "<< j << " " << tmp << "\n";				
					}
				}
				
			}
		}
				
	}
}



int try_insertion_euclidDist(Worker &w, int rid)
{
	return 1;
	Request& r = R[rid];
	vector<int>& schedule = w.S;
	vector<double>& reach = w.reach;
	vector<double>& ddl = w.ddl;

	Position pw,ps,pe;
	pw.x = vertices[w.pid].x, pw.y = vertices[w.pid].y;
	ps.x = vertices[r.s].x, ps.y = vertices[r.s].y;
	pe.x = vertices[r.e].x, pe.y = vertices[r.e].y;

	vector<Position> schedulePostions;
	for (int x = 0; x < schedule.size(); x++)
	{
		Position p;
		p.x = vertices[Pos(schedule[x])].x, p.y = vertices[Pos(schedule[x])].y;
		schedulePostions.push_back(p);
	}

	//local cache
	double  dr = RealDistance(ps,pe)/MAX_speed;
	vector<double> scache, ecache;
	scache.resize(schedule.size()+1, 0), ecache.resize(schedule.size()+1, 0);
	scache[0] = RealDistance(pw,ps)/MAX_speed, ecache[0] = RealDistance(pw,pe)/MAX_speed;
	for (int i = 0; i < schedulePostions.size(); i++)
	{
		scache[i+1] = RealDistance(schedulePostions[i], ps);
		ecache[i+1] = RealDistance(schedulePostions[i], pe);
	}

	double tmp = 0, tmpi = 0 ,tmpj = 0;
	if (w.S.empty()) {
		tmp = w.tim + scache[0] + dr;
		if (tmp < r.ddl + EPS && r.com <= w.cap) { return 1;}
	}

	vector<double> detoucache; //detoucache[j] is the minimal value insert i at [0,j]
	for (int j = 0; j <= schedule.size(); j++)
	{
		double det;
		if (j == 0)
		{
			det = scache[j] + scache[j+1];
		}
		else if (j == schedule.size())
		{
			det = scache[j];
			if (det > detoucache[detoucache.size()-1])
			{
				det = detoucache[detoucache.size()-1];
			}
		}
		else
		{
			det = scache[j] + scache[j+1];
			if (det > detoucache[detoucache.size()-1])
			{
				det = detoucache[detoucache.size()-1];
			}
		}
		
		detoucache.push_back(det);
	}

	
	for (int j = 0; j <= schedule.size(); j++)
	{
		if (j == schedule.size())
		{
			if (detoucache[j] + dr + ecache[j] < r.ddl)
			{
				return 1;
			}
		}
		else
		{
			double det = detoucache[j] + ecache[j+1];
			if (j!=0 & ecache[j]!=ecache[j-1]) //r.s also insert atj, i == j
			{
				det += dr;
			}
			else
			{
				det += ecache[j];
			}

			if (det + reach[j] < ddl[j])
			{
				return 1;
			}
						
		}
		
		return -1;
	}
}


void insertion(Worker &w, int rid, int wid, int optimal_i, int optimal_j){
	if (w.S.empty())
	{

		w.S.push_back(rid << 1);
		w.S.push_back(rid << 1 | 1);
		//updateDriverArr(w);

		cout << "insert at empty car: "<< wid << endl;
		cout << "before insertion: "<<endl;
		cout <<  w.S << endl;
		cout <<  w.reach << endl;
		
		loginsertion << "insert at empty car: "<< wid << "\n";
		loginsertion << "before insertion: " << "\n";
		for (int k = 0; k < w.S.size(); k++)
		{
			loginsertion << w.S[k] << " ";
		}
		loginsertion << "\n" ;
		for (int k = 0; k < w.reach.size(); k++)
		{
			loginsertion << w.reach[k] << " ";
		}
		loginsertion << "\n" ;
		

		updateDriverArr(w);

		cout << "after insertion: " << optimal_i << " " << optimal_j <<endl;
		cout <<  w.S << endl;
		cout <<  w.reach << endl;
		cout << w.picked <<endl;
		cout << "------------------------------------" << endl;


		loginsertion << "after insertion: " << optimal_i << " " << optimal_j << "\n";
		for (int k = 0; k < w.S.size(); k++)
		{
			loginsertion << w.S[k] << " ";
		}
		loginsertion << "\n" ;
		for (int k = 0; k < w.reach.size(); k++)
		{
			loginsertion << w.reach[k] << " ";
		}
		loginsertion << "\n" ;
		loginsertion << "------------------------------------" << "\n";

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

	loginsertion << "Assign request:" << rid <<" "<< R[rid].ddl << " to worker : " << wid << "\n";
	loginsertion << "before insertion: " << "\n";
	for (int k = 0; k < w.S.size(); k++)
	{
		loginsertion << w.S[k] << " ";
	}
	loginsertion << "\n" ;
	for (int k = 0; k < w.reach.size(); k++)
	{
		loginsertion << w.reach[k] << " ";
	}
	loginsertion << "\n" ;

	w.S.clear();
	w.S = ret;
	updateDriverArr(w);

    cout << "after insertion: " << optimal_i << " " << optimal_j <<endl;
    cout <<  w.S << endl;
    cout <<  w.reach << endl;
    cout << "------------------------------------" << endl;

	loginsertion << "after insertion: " << optimal_i << " " << optimal_j << "\n";
	for (int k = 0; k < w.S.size(); k++)
	{
		loginsertion << w.S[k] << " ";
	}
	loginsertion << "\n" ;
	for (int k = 0; k < w.reach.size(); k++)
	{
		loginsertion << w.reach[k] << " ";
	}
	loginsertion << "\n" ;
	loginsertion << "------------------------------------" << "\n";
}


void updateDriverArr(Worker& w){
	double tim = w.tim;
	vector<double>& reach = w.reach;
	vector<int>& schedule = w.S;
	vector<PLF>& FunctionTimesegments = w.FunctionTimesegments;
	
	FunctionTimesegments.clear();
	reach.clear();

	//update F function
	for (int x = 1; x < schedule.size(); x++)
	{
		PLF plf;
		tree.PLCst(Pos(schedule[x-1]), Pos(schedule[x]), tim, TMAX, plf);
		FunctionTimesegments.push_back(plf);
	}
	for (int x = 1; x < schedule.size() - 1; x++)
	{
		PLF PLFxn; // from x to n
		PLFxn = FunctionTimesegments[x-1];
		for (int y = x+1; y < schedule.size(); y++)
		{
			PLF PLFy;
			FunctionTimesegments[y-1].compound(PLFxn, PLFy, schedule[y-1]);
			PLFxn = PLFy;			
		}
		FunctionTimesegments.push_back(PLFxn);		
	}
	
	
	// o(n) update reach
	tim += tree.tdsp(w.pid, Pos(schedule[0]), tim, false);
	reach.push_back(tim);
	for (int x = 1; x < schedule.size(); x++)
	{
		tim = FunctionTimesegments[x-1].dpt2arr(tim);
		reach.push_back(tim);
	}

	//o(n) uodate ddl
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
			
			auto function = FunctionTimesegments[k].f;
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


	// //according reach&ddl to segment the F function
	for (int k = 0; k < FunctionTimesegments.size(); k++)
	{
		double arr;
		double d;
		if (k<schedule.size()-1)
		{
			arr = reach[k];
			d = ddl[k+1];
		}
		else
		{
			arr = reach[k - (schedule.size() - 1)];
			d = ddl[ddl.size() - 1];
		}

		while ((FunctionTimesegments[k].f->begin()+1) -> t < arr & FunctionTimesegments[k].f->size() >= 3)
		{
			FunctionTimesegments[k].f->erase(FunctionTimesegments[k].f->begin());
		}

		double t_, w_;
		auto it1 = FunctionTimesegments[k].f->begin();
		auto it2 = it1+1;
		double slop = (it2->w - it1->w)/(it2->t - it1->t);
		it1->w = slop*(arr - it1->t) + it1->w;
		it1->t = arr;

		int s = FunctionTimesegments[k].f->size() - 1;
		while ((FunctionTimesegments[k].f->begin() + s-1) ->t > d & FunctionTimesegments[k].f->size() >= 3 )
		{
			FunctionTimesegments[k].f->erase(FunctionTimesegments[k].f->begin() + s);
			s--;
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
	for (int k = 0; k < FunctionTimesegments.size(); k++)
	{
		auto s1 = FunctionTimesegments[k].f->begin();
		while (s1 < FunctionTimesegments[k].f->end())
		{
			memo += sizeof(s1->t);
			memo += sizeof(s1->w);
			s1 ++;
		}
	}
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

	ofstream dump_trajectory(trajectory_n_path);
	for (int i = 0; i < m; i++)
	{
		Worker& w = W[i];
		dump_trajectory << "Taxi: " << i << "  #"<< w.maxnum <<"  M"<< w.memory <<"\n";
		for(auto& tuple: w.trajectory){
			//cout << "Arrival node: "<<get<0>(tuple) << "  at time: " << get<2>(tuple) << endl;
			dump_trajectory << "Arrival node: "<<get<0>(tuple)<< "  at time: " << get<1>(tuple) <<"\n";
		}
	}
	dump_trajectory.close();


	ofstream dump_result(dump_resultn_path);
	dump_result << numOfTDSP << " " << numofTDSPInitial <<"\n";
	dump_result << carnum << "  " << cartime << "\n";
	dump_result << assigntime << "  " << assignnum << "\n";
	dump_result << insertiontime << "  " << insertionnum << "\n";
	dump_result << serverPassenger << "\n";
	dump_result.close();
}

void test()
{
	ofstream requestsDDL("./data/requestsDisDDL.txt");

    cout << "begin request" << endl;
	ifstream ifs;
	ifs.open(requests_path.c_str());
	if (!ifs.is_open()) {
		printf("%s does not exist\n", requests_path.c_str());
		exit(0);
	}
	ifs >> n;
	double t, len,dis;
	int s,e,cap;
	for (int i = 0; i < n; ++ i) {
		//cout << i << endl;
		ifs >> t >> s >> e >> len >>cap;
		dis = RealDistance(vertices[s].y, vertices[s].x, vertices[e].y, vertices[e].x);
		double nlen = tree.tdsp(s,e,t,false);
		cout << s << " " << e << " " << t << " " << nlen <<" "<< dis << endl;
		requestsDDL << t  << " " << s << " " << e << " " << dis<<" "<<nlen <<" "<< cap <<"\n";

	}
	ifs.close();
	requestsDDL.close();
	cout << "end request" << endl;	
}


void testquery()
{
	for (int i = 0; i < n; ++ i) {
		double x = tree.tdsp(R[i].s,R[i].e,R[i].tim,false);
		//cout << R[i].s <<" "<<R[i].e<<" " <<R[i].tim<<" "<< x <<endl;
	}	
}


void testplf()
{
	PLF PLFy;
	tree.PLCst(1, 1000, 100, TMAX, PLFy);
	cout << PLFy << endl;

	// auto s1 = PLFy.f->begin();
	// while (s1 < PLFy.f->end())
	// {
	// 	cout << s1 -> t << endl;
	// 	s1 ++;
	// }
	int a = 100;
	cout << sizeof (a) <<endl;
	double b = 100;
	cout << sizeof(b) <<endl;
	
}


/*  ----------------------------------------  try_insertion -----------------------------
// void try_insertion(Worker &w, int rid, double &delta, int &optimanl_i, int &optimanl_j) {
// 	Request& r = R[rid];
// 	double opt = INF;

// 	if (w.S.empty()) {
// 		double tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
// 		numOfTDSP++; 
//         tmp += tree.tdsp(r.s, r.e, tmp,false);
// 		numOfTDSP++;
// 		if (tmp < r.ddl + EPS && r.com <= w.cap) {
// 			delta = tmp, optimanl_i = 0, optimanl_j = 0;
// 		}
// 		return;
// 	}
	

// 	vector<int>& picked = w.picked;
// 	vector<int>& schedule = w.S;
// 	vector<double>& reach = w.reach;
// 	vector<double>& ddl = w.ddl;
// 	int delta_i = -1;
	
// 	int capacityflag = 0;

// 	double arrsj = INF; 			// after inserting r.s at delta_i, the new arrival time at j 
// 	for (int j = 0; j <= schedule.size(); j++)
// 	{
// 		double tmp=INF, tmps=INF, det_i = INF;

// 		loginsertion << "*****************************************" << "\n";
// 		loginsertion << "try_insertion: " << j << "\n";

// 		if (j == schedule.size())
// 		{
// 			if(r.com > w.cap){
// 				loginsertion <<"cap false j: "<<r.com << "\n";
// 				continue;
// 			} 

// 			//case 1: i in general, j == schedule.size()
// 			if(delta_i != -1)
// 			{
// 				//cout << "case: i < j, j == schedule.size()" <<endl;
// 				//cout << "tree.tdsp(Pos(schedule[j-1]), r.e, arrsj, false); "<<schedule[j-1] <<" "<<Pos(schedule[j-1])<<" "<< r.e <<" " << arrsj <<endl;
// 				tmp = arrsj + tree.tdsp(Pos(schedule[j-1]), r.e, arrsj, false);
// 				numOfTDSP++;
// 				if (tmp <= r.ddl & tmp < delta)
// 				{
// 					delta = tmp;
// 					optimanl_i = delta_i;
// 					optimanl_j = j;
// 					loginsertion <<"update delta"<< tmp << " " << optimanl_i <<"  "<<optimanl_j<< "\n";
// 				}
// 				else if(tmp > r.ddl)
// 				{
// 					loginsertion <<"ddl&cap false S.size()"<< tmp << " " << r.ddl << "\n";
// 					continue;
// 				}
// 			}

// 			//case 2: i == j == schedule.size()
// 			//cout << "case: i == j == schedule.size()" <<endl;
// 			//cout << "tree.tdsp(Pos(schedule[j-1]), r.s, reach[j-1], false); "<<schedule[j-1] <<" "<<Pos(schedule[j-1])<<" "<< r.s <<" " << reach[j-1] <<endl;
// 			tmp = reach[j-1] + tree.tdsp(Pos(schedule[j-1]), r.s, reach[j-1], false);
// 			//cout << "tree.tdsp(r.s, r.e, tmp, false); " <<r.s<<" "<< r.e <<" " << tmp <<endl;
// 			tmp += tree.tdsp(r.s, r.e, tmp, false);
// 			numOfTDSP++;
// 			numOfTDSP++;
// 			if (tmp <= r.ddl & tmp < delta)
// 			{
// 				delta = tmp;
// 				optimanl_i = j;
// 				optimanl_j = j;
// 				loginsertion <<"update delta, new i for j == s.size()"<< tmp << " " << optimanl_i <<"  "<<optimanl_j<< "\n";
// 			}
			
// 			continue;			
// 		}

// 		// cout << " ----------------------------------- " << endl;
// 		// cout << j << " " << w.num << " " << r.com << " " << w.cap << endl;
// 		if (j == 0)
// 		{	
// 			if(w.num + r.com > w.cap){
// 				loginsertion <<"cap false j: "<<j<<" "<< picked[0] << " " << r.com << "\n";
// 				continue;
// 			}
// 			//cout << "case j == 0, arr[r.s]" <<endl;
// 			tmps = w.tim + tree.tdsp(w.pid, r.s, w.tim, false);
// 			numOfTDSP++; 
// 		}
// 		else
// 		{
// 			if(picked[j-1] + r.com > w.cap){
// 				loginsertion <<"cap false j: "<<j<<" "<< picked[j] << " " << r.com << "\n";
// 				continue;
// 			} 
// 			//cout << "case j in general, arr[r.s]" <<endl;
// 			tmps = reach[j-1] + tree.tdsp(Pos(schedule[j-1]), r.s, reach[j-1], false);
// 			numOfTDSP++;
// 		}
		

// 		// if insert r.s at j donot validate the capacity constrain
// 		// if (tmps < INF)
// 		// {
// 		// cout << "case: r.s -> r.e" <<endl;
// 		// cout << "tree.tdsp(r.s, r.e, tmps, false); "<<r.s <<" "<<r.e<<" "<< tmps <<endl;
// 		tmp = tmps + tree.tdsp(r.s, r.e, tmps, false);
// 		numOfTDSP++;

// 		// cout << "case: r.s -> scedule[j]" <<endl;
// 		// cout << "tree.tdsp(r.s, Pos(schedule[j]), tmps, false); "<<r.s <<" "<<schedule[j]<<" "<< tmps <<endl;
// 		det_i = tmps + tree.tdsp(r.s, Pos(schedule[j]), tmps, false);
// 		numOfTDSP++;
// 		// }
		
// 		// cout << "check new position for i, r.s: " <<endl;
// 		// cout << "function_arrival(w, j-1, j, arrsj): "<<arrsj <<endl;
// 		cout << "j: " << j << endl;
// 		cout << "det_i<=ddl[j] && tmp < r.ddl && ((delta_i == -1)||(delta_i!=-1 && det_i < function_arrival(w, j-1, j, arrsj)))" <<endl;
// 		cout << "delta_i: " << delta_i << " " << "det_i: " << det_i <<" " <<"arrsj: "<< arrsj <<" " << "function_arrival: " <<function_arrival(w, j-1, j, arrsj) <<endl;
// 		loginsertion << "det_i<=ddl[j] && tmp < r.ddl && ((delta_i == -1)||(delta_i!=-1 && det_i < function_arrival(w, j-1, j, arrsj)))" <<"\n";
// 		loginsertion << "delta_i: " << delta_i << " " << "det_i: " << det_i <<" " <<"arrsj: "<< arrsj <<" " << "function_arrival: " <<function_arrival(w, j-1, j, arrsj) <<"\n";
// 		if (det_i<=ddl[j])
// 		{
// 			loginsertion << "true1; ";
// 			cout << "true1; ";
// 		}
// 		else
// 		{
// 			loginsertion << "faluse1; ";
// 			cout << "faluse1; ";
// 		}
// 		if (tmp < r.ddl)
// 		{
// 			loginsertion << "true2; ";
// 			cout << "true2; ";
// 		}
// 		else
// 		{
// 			loginsertion << "faluse2; ";
// 			cout << "faluse2; ";
// 		}	
// 		if (((delta_i == -1)||(delta_i!=-1 && det_i + 0.001 < function_arrival(w, j-1, j, arrsj))))
// 		{
// 			loginsertion << "true3; ";
// 			cout << "true3; ";
// 		}
// 		else
// 		{
// 			loginsertion << "faluse3; ";
// 			cout << "faluse3; ";
// 		}
// 		cout << endl;
// 		loginsertion << "\n ";
		
		 
// 		if (det_i<=ddl[j] && tmp < r.ddl && ((delta_i == -1)||(delta_i!=-1 && det_i + 0.001 < function_arrival(w, j-1, j, arrsj))))
// 		{
// 			// cout << "delta_i = j: " << delta_i << endl;
// 			capacityflag = 1;
// 			if(j-1 >= 0)
// 			{
// 				cout << w.FunctionTimesegments.size() << endl;
// 				loginsertion << w.FunctionTimesegments.size() << " " << j-1 << "\n";
// 				auto s1 = w.FunctionTimesegments[j-1].f->begin();
// 				while (s1 != w.FunctionTimesegments[j-1].f->end())
// 				{
// 					loginsertion << s1->t << " " << s1->w <<"; ";
// 					s1 ++;
// 				}
// 				loginsertion << "\n" ;
// 			}


// 			loginsertion <<"update new i for j =  "<<j<<" "<<det_i<<" "<<ddl[j]<< "\n";
// 			loginsertion << Pos(w.S[j-1]) << " " << Pos(w.S[j]) << " " << arrsj << "\n";
// 			delta_i = j;
// 			arrsj = det_i;
// 			if (tmp < r.ddl)
// 			{
// 				//cout << "update i, i == j:" << endl;
// 				//cout << "tree.tdsp(r.e, Pos(schedule[j]), tmp, false);" << r.e << " " << schedule[j] << " " << tmp << endl;
// 				tmp += tree.tdsp(r.e, Pos(schedule[j]), tmp, false);
// 				numOfTDSP++;

// 				tmp = function_arrival(w,j,schedule.size()-1, tmp);
// 				if (tmp < delta)
// 				{
// 					delta = tmp;
// 					optimanl_i = j;
// 					optimanl_j = j;
// 					loginsertion <<"update delta"<< tmp << " " << optimanl_i <<"  "<<optimanl_j<< "\n";
// 				}				
// 			}
// 		}

// 		else 
// 		{
// 			//cout << " i not change: " <<delta_i<< endl;
// 			if (capacityflag == 0)
// 			{
// 				continue;
// 			}
			

// 			loginsertion <<"donot update i  "<<delta_i<< "\n";
// 			if (delta_i == -1)
// 			{
// 				arrsj = reach[j];
// 			}
// 			else
// 			{
// 				//cout << "case: i<j" << endl;
// 				//cout << "tree.tdsp(Pos(schedule[j-1]), r.e, arrsj, false); " << schedule[j-1] << " " << r.e << " "<<arrsj << endl;
// 				if (picked[j-1] + r.com > w.cap)
// 				{
// 					capacityflag = 0;
// 					continue;
// 				}
				
				
// 				tmp = arrsj + tree.tdsp(Pos(schedule[j-1]), r.e, arrsj, false);
// 				arrsj = function_arrival(w, j-1, j, arrsj);
// 				numOfTDSP++;
// 				if (tmp <= r.ddl)
// 				{
// 					//cout << "case: i<j, " << endl;
// 					//cout << "tree.tdsp(r.e, Pos(schedule[j]), tmp, false); " << r.e << " " << schedule[j] << " "<<tmp << endl;
// 					tmp += tree.tdsp(r.e, Pos(schedule[j]), tmp, false);
// 					numOfTDSP++;
// 					if (tmp > ddl[j])
// 					{
// 						loginsertion <<"ddl false j: "<<j<<" "<<tmp << " "<<ddl[j]<< "\n";
// 						continue;
// 					}
					
// 					//cout << "function_arrival(w,j,schedule.size()-1, tmp); " <<tmp << endl;
// 					tmp = function_arrival(w,j,schedule.size()-1, tmp);
// 					if (tmp < delta)
// 					{
// 						delta = tmp;
// 						optimanl_i = delta_i;
// 						optimanl_j = j;
// 						loginsertion <<"update delta"<< tmp << " " << optimanl_i <<"  "<<optimanl_j<< "\n";
// 					}
					
// 				}
// 				//arrsj = function_arrival(w, j-1, j, arrsj);
// 			}
			
// 		}
		
// 	}
// }
*/













// //try insertion with local cache
// void try_insertion(Worker &w, int rid, double &delta, int &optimanl_i, int &optimanl_j) {
// 	Request& r = R[rid];
// 	vector<int>& picked = w.picked;
// 	vector<int>& schedule = w.S;
// 	vector<double>& reach = w.reach;
// 	vector<double>& ddl = w.ddl;	
	
// 	// local cache
// 	vector<double> scach,cachs, ecach,cache, cachse, secach;
// 	scach.resize(schedule.size()+1, 0), cachs.resize(schedule.size()+1, 0);
// 	ecach.resize(schedule.size()+1, 0), cache.resize(schedule.size()+1, 0);	
// 	cachse.resize(schedule.size()+1, 0), secach.resize(schedule.size()+1, 0);

// 	cachs[0] = tree.tdsp(w.pid, r.s, w.tim,false);
// 	cahce[0] =  tree.tdsp(w.pid, r.e, w.tim,false);
// 	cachse[0] = tree.tdsp(r.s, r.e, cachs[0],false);
// 	for (int i = 0; i < schedule.size(); i++)
// 	{
// 		scach[i] = tree.tdsp(r.s, Pos(schedule[i]),cachs[i],false);
// 		ecach[i] = tree.tdsp(r.e, Pos(schedule[i]),cahce[i],false);
// 		secach[i] = tree.tdsp(r.e, Pos(schedule[i]),cachse[i],false);

// 		cachs[i+1] = tree.tdsp(Pos(schedule[i]),r.s,reach[i],false);
// 		cache[i+1] = tree.tdsp(Pos(schedule[i]),r.e,reach[i],false);
// 		cachse[i+1] = tree.tdsp(r.s, r.e, cachs[i+1], false);
// 	}
	
// 	double opt = INF;
// 	int delta_i = -1;
// 	if (w.S.empty())
// 	{
// 		double tmp = w.tim + cachs[0] + cachse[0];
// 		numOfTDSP++;
// 		numOfTDSP++;
// 		if (tmp < r.ddl + EPS && r.com <= w.cap) {
// 			delta = tmp, optimanl_i = 0, optimanl_j = 0;
// 		}
// 		return;	
// 	}

// 	double arrsj = INF; // after inserting r.s at delta_i, the new arrival time at j 
// 	for (int j = 0; j <= schedule.size(); j++)
// 	{
// 		double tmp=INF, tmps=INF, det_i = INF;
// 		if (j == schedule.size())
// 		{
// 			//case 1: i in general, j == schedule.size()
// 			if(delta_i != -1)
// 			{
// 				tmp = arrsj + tree.tdsp(Pos(schedule[j-1]), r.e, arrsj, false);
// 				numOfTDSP++;
// 				if (tmp <= r.ddl & tmp < delta)
// 				{
// 					delta = tmp;
// 					optimanl_i = delta_i;
// 					optimanl_j = j;
// 				}
// 			}

// 			//case 2: i == j == schedule.size()
// 			tmp = reach[j-1] + cachs[j] + cachse[j];
// 			numOfTDSP++;
// 			numOfTDSP++;
// 			if (tmp <= r.ddl & tmp < delta)
// 			{
// 				delta = tmp;
// 				optimanl_i = j;
// 				optimanl_j = j;
// 			}
			
// 			continue;			
// 		}

// 		if (j == 0 & w.num + r.com <= w.cap)
// 		{ 
// 			tmps = w.tim + cachs[0];
// 			numOfTDSP++; 
// 		}
// 		else
// 		{ 
// 			if (picked[j-1]+r.com <= w.cap)
// 			{
// 				tmps = reach[j-1] + cachs[j];
// 				numOfTDSP++;
// 			}
// 		}
// 		// if insert r.s at j donot validate the capacity constrain
// 		if (tmps < INF)
// 		{
// 			tmp = tmps + cachse[j];
// 			numOfTDSP++;
// 			det_i = tmps + scach[j];
// 			numOfTDSP++;
// 		}
		
	
// 		if (det_i<=ddl[j] & tmps < r.ddl & ((delta_i == -1)|(delta_i!=-1 & det_i < function_arrival(w, j-1, j, arrsj))))
// 		{
// 			delta_i = j;
// 			arrsj = det_i;
// 			if (tmp < r.ddl)
// 			{
// 				//tmp += tree.tdsp(r.e, Pos(schedule[j]), tmp, false);
// 				tmp += ecach[j];
// 				numOfTDSP++;
// 				tmp = function_arrival(w,j,schedule.size()-1, tmp);
// 				if (tmp < delta)
// 				{
// 					delta = tmp;
// 					optimanl_i = j;
// 					optimanl_j = j;
// 				}				
// 			}
// 		}

// 		else 
// 		{
// 			if (delta_i == -1)
// 			{
// 				arrsj = reach[j];
// 			}
// 			else // delta_i != -1 
// 			{
// 				//tmp = arrsj + tree.tdsp(Pos(schedule[j-1]), r.e, arrsj, false);
// 				tmp = arrsj + cache[j];
// 				numOfTDSP++;
// 				if (tmp <= r.ddl)
// 				{
// 					//tmp += tree.tdsp(r.e, Pos(schedule[j]), tmp, false);
// 					tmp += ecach[j];
// 					numOfTDSP++;
// 					tmp = function_arrival(w,j,schedule.size()-1, tmp);
// 					if (tmp < delta)
// 					{
// 						delta = tmp;
// 						optimanl_i = delta_i;
// 						optimanl_j = j;
// 					}
					
// 				}
// 				arrsj = function_arrival(w, j-1, j, arrsj);
// 			}
			
// 		}		
// 	}	
// }


// void updateDriverArr(Worker& w){
// 	double tim = w.tim;
// 	vector<double>& reach = w.reach;
// 	vector<int>& schedule = w.S;
// 	vector<PLF>& FunctionTimesegments = w.FunctionTimesegments;
	
// 	FunctionTimesegments.clear();
// 	reach.clear();

// 	//update F function
// 	for (int x = 1; x < schedule.size(); x++)
// 	{
// 		PLF plf;
// 		tree.PLCst(Pos(schedule[x-1]), Pos(schedule[x]), tim, TMAX, plf);
// 		FunctionTimesegments.push_back(plf);
// 	}
// 	for (int x = 1; x < schedule.size() - 1; x++)
// 	{
// 		PLF PLFxn; // from x to n
// 		PLFxn = FunctionTimesegments[x-1];
// 		for (int y = x+1; y < schedule.size(); y++)
// 		{
// 			PLF PLFy;
// 			FunctionTimesegments[y-1].compound(PLFxn, PLFy, schedule[y-1]);
// 			PLFxn = PLFy;			
// 		}
// 		FunctionTimesegments.push_back(PLFxn);		
// 	}
	
	
// 	// o(n) update reach
// 	tim += tree.tdsp(w.pid, Pos(schedule[0]), tim, false);
// 	reach.push_back(tim);
// 	for (int x = 1; x < schedule.size(); x++)
// 	{
// 		tim = FunctionTimesegments[x-1].dpt2arr(tim);
// 		reach.push_back(tim);
// 	}

// 	//o(n) uodate ddl
// 	vector<double>& ddl = w.ddl;
// 	ddl.clear();
// 	for (int k = w.S.size() - 1; k >= 0; --k)
// 	{
// 		if(k == w.S.size() - 1){
// 			ddl.insert(ddl.begin(), DDLEndPos(w.S[k]));

// 		}
// 		else
// 		{
// 			double delat = INF;
// 			double opt = INF;
// 			if(w.S[k] & 1){
// 				opt = DDLEndPos(w.S[k]);
// 			}

// 			auto s1 = FunctionTimesegments[k].f->begin();
			
// 			if (s1->t + s1->w >= ddl[0])
// 			{
// 				delat = ddl[0] - s1->w;
// 			}
// 			else{
// 				while (s1->t + s1->w < ddl[0])
// 				{
// 					s1++;
// 				}
// 				auto s2 = s1 - 1;

// 				if (s1 == FunctionTimesegments[k].f->end())
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


// 	// //according reach&ddl to segment the F function
// 	for (int k = 0; k < FunctionTimesegments.size(); k++)
// 	{
// 		double arr;
// 		double d;
// 		if (k<schedule.size()-1)
// 		{
// 			arr = reach[k];
// 			d = ddl[k+1];
// 		}
// 		else
// 		{
// 			arr = reach[k - (schedule.size() - 1)];
// 			d = ddl[ddl.size() - 1];
// 		}

// 		while ((FunctionTimesegments[k].f->begin()+1) -> t < arr & FunctionTimesegments[k].f->size() >= 3)
// 		{
// 			cout << " FunctionTimesegments[k].f->begin() + 1 " << endl;
// 			FunctionTimesegments[k].f->erase(FunctionTimesegments[k].f->begin());
// 		}

// 		double t_, w_;
// 		auto it1 = FunctionTimesegments[k].f->begin();
// 		auto it2 = it1+1;
// 		double slop = (it2->w - it1->w)/(it2->t - it1->t);
// 		it1->w = slop*(arr - it1->t) + it1->w;
// 		it1->t = arr;

// 		int s = FunctionTimesegments[k].f->size() - 1;
// 		while ((FunctionTimesegments[k].f->begin() + s-1) ->t > d & FunctionTimesegments[k].f->size() >= 3 )
// 		{
// 			cout << " FunctionTimesegments[k].f->begin() + s " << endl;
// 			FunctionTimesegments[k].f->erase(FunctionTimesegments[k].f->begin() + s);
// 			s--;
// 		}	
// 	}
	
// 	// o(n) update picked
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