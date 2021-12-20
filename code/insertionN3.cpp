# include "insertionN3.h"

#include "TDGT.h"
#include "io.h"

int nV, m, c, n,numOfVertices;
double gridL, alpha;
int pos = 0;
Request* R = NULL;
Worker* W = NULL;
vector<vertex> vertices;
TDGT tree;

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
		//R[i].ddl = TMAX - 1;
		//R[i].ddl = R[i].tim + R[i].len * delta;
		R[i].ddl = R[i].tim + R[i].len + deltasecond;
		if (R[i].ddl > TMAX)
		{
			R[i].ddl  = TMAX - 1;
		}
	}
	ifs.close();
	cout << "end request" << endl;


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
        //vertexFile >> vertices[i].x >> vertices[i].y;
        vertexFile >> vertices[i].y >> vertices[i].x;
		vertices[i].id = i;
    }
    vertexFile.close();	
	cout << "end coordinates" << endl;
	
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
		cout << "updateDriver finished" <<endl;
		vector<int> car;
		car = single_search(R[pos].s, R[pos].ddl-R[pos].tim);
		//cout <<" car num: " << car.size() << endl;
		carnum += car.size();
		cartime ++ ;
		assignnum ++;
		cout << "Assign request: " << pos <<"  with car size: "<< car.size() << endl;
		cout << W[0].S << endl;
		cout << R[pos].s << " "<< R[pos].e << " " << endl;
		assignbegin = clock();
		assignTaxi(car);
		assignend = clock();
		cout << "**********************************************************************" <<endl;
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
				//cout << "try insertion time: " << double(insertionend - insertionbegin) / CLOCKS_PER_SEC << endl;
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
		//cout << "Assigin request: " << pos << endl;
        insertion(W[id], pos, id, optimal_i, optimal_j);
		serverPassenger++;
	}
	pos++;
}

int try_insertion_euclidDist(Worker &w, int rid)
{
	return 1;
	Request& r = R[rid];
	vector<int>& schedule = w.S;
	vector<double>& reach = w.reach;

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
		int flag = 1;  // insert at j is fessible or not
		if (j == schedule.size())
		{
			if (detoucache[j] + ecache[j] < r.ddl)
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
			
			for (int k = j; k < schedule.size(); k++)
			{
				if (schedule[k] & 1)
				{
					if (reach[k] + det > DDLEndPos(schedule[k]))
					{
						flag = -1;
						break;
					}
					
				}
				
			}
			if (flag)
			{
				return flag; //flag not change in k >= j all satisfy
			}

		}
	}
	return -1;
}


void try_insertion(Worker &w, int rid, double &delta, int &optimanl_i, int &optimanl_j) {
	Request& r = R[rid];
	double opt = INF;
	
	delta = INF;
	if (w.S.empty()) {
		double tmp = w.tim + tree.tdsp(w.pid, r.s, w.tim,false);
		numOfTDSP++; 
        tmp += tree.tdsp(r.s, r.e, tmp,false);
		numOfTDSP++;
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
				}
			
			}else if (i == w.S.size()){ //case 2: i==j==w.S.size()
				tmp = reach[w.reach.size() - 1] + tree.tdsp(Pos(schedule[w.S.size() - 1]), r.s, reach[w.reach.size() - 1], false);
				numOfTDSP++;
				tmp += tree.tdsp(r.s, r.e, tmp,false);
				numOfTDSP++;
				if (tmp > r.ddl || picked[i-1]+r.com > w.cap){false_flag = 1;}
				if(tmp<delta and false_flag == 0){
					delta = tmp, optimanl_i = i, optimanl_j = j;  				
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
						numOfTDSP++;
						// dump_result << "try_insertion case 5, j ==w.s.size() :" << "\n"; 
						if(tmp>r.ddl){false_flag = 1;}
					}
					if (tmp < delta and false_flag == 0){
						delta = tmp, optimanl_i = i, optimanl_j = j; 
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
		//updateDriverArr(w);

		cout << "insert at empty car: "<< wid << endl;
		cout << "before insertion: "<<endl;
		cout <<  w.S << endl;
		cout <<  w.reach << endl;

		updateDriverArr(w);

		cout << "after insertion: " << optimal_i << " " << optimal_j <<endl;
		cout <<  w.S << endl;
		cout <<  w.reach << endl;
		cout << "------------------------------------" << endl;
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

    cout << "insert at car: " << wid << endl;
    cout << "before insertion: "<<endl;
    cout <<  w.S << endl;
    cout <<  w.reach << endl;


	w.S.clear();
	w.S = ret;
	updateDriverArr(w);

    cout << "after insertion: " <<"  positions: "<< optimal_i << " " << optimal_j <<endl;
    cout <<  w.S << endl;
    cout <<  w.reach << endl;
    cout << "------------------------------------" << endl;
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

	ofstream dump_trajectory(trajectory_n3_path);
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


	ofstream dump_result(dump_resultn3_path);
	dump_result << numOfTDSP << " " << numofTDSPInitial <<"\n";
	dump_result << carnum << "  " << cartime << "\n";
	dump_result << assigntime << "  " << assignnum << "\n";
	dump_result << insertiontime << "  " << insertionnum << "\n";
	dump_result << serverPassenger << "\n";
	dump_result.close();
}