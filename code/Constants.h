#ifndef TDGT_CONSTANTS_H
#define TDGT_CONSTANTS_H

#include <string>
#include <map>
#include <iostream>
#include <array>
#include <limits>
#include <cmath>
#include <cassert>


#define DE_INTV -2147483648
#define INTV_CNTED -2147483647
#define DE_W 2147483647


static double EPSILON = 1e-5; //1/86400 =1.157e-5

static int TMAX = 86400;
//static int TMAX = 8640000;

// hyper parameters
static unsigned long LEAF_SIZE = 64;
static unsigned long FANOUT = 4;

static int delta = 1.2;
static double deltasecond = 50000;
static double MAX_speed = 26;
static double grid_value = 1000;
static double len = 600;
////////////////////////////////////////////////////////////////
const std::string tgraph_path = "../Data/CAL.tpgr";
const std::string tdgt_idx_path = "../Data/CAL.idx";
const std::string tdsp_query_path = "../Data/CALTDSP.demands";
const std::string tisp_query_path = "../Data/CALTISP.demands";
////////////////////////////////////////////////////////////////

// const std::string workers_path = "/homes/zgongae/data/chengdu/chengduworkersDebug5k.txt";
// const std::string requests_path = "/homes/zgongae/data/chengdu/chengdurequestsDebug5k.txt";
const std::string workers_path = "./data/workers.txt";
const std::string requests_path = "./data/requests.txt";

const std::string vertexes_path = "/homes/zgongae/data/chengdu/Chengdu.co";

const std::string TDGT_tree_path = "/home/data/1TB/zgongae/data/tree_chengdu1220.txt";

const std::string dump_resultn3_path = "./data/dumpResultsn3.txt";
const std::string trajectory_n3_path = "./data/trajectoryn3.txt";

const std::string dump_resultn2_path = "./data/dumpResultsn2.txt";
const std::string trajectory_n2_path = "./data/trajectoryn2.txt";

const std::string dump_resultn_path = "./data/dumpResultsn.txt";
const std::string trajectory_n_path = "./data/trajectoryn.txt";

const std::string insertion_log_n = "./data/logn.txt";
const std::string insertion_log_n2 = "./data/logn2.txt";
const std::string insertion_log_n3 = "./data/logn3.txt";

#endif //TDGT_CONSTANTS_H
