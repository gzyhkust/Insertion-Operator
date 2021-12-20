//
// Created by 王勇 on 2018/7/6.
//

#ifndef TDGT_TDGT_H
#define TDGT_TDGT_H

#include <iostream>
#include <cstdio>
#include <unordered_map>
#include <queue>
#include <utility>
#include <vector>
#include "/homes/zgongae/metis-5.1.0/include/metis.h"
//#include <metis.h>
//#include "metis.h"
#include <algorithm>
#include <ctime>
#include <cmath>
#include <climits>
#include <cstring>
#include <memory>
#include <sstream>
#include "PLF.h"
#include "Constants.h"
#include "TGraph.h"
#include "TNode.h"


using namespace std;

struct TDGT {
    int nodes_num;              //nodes num
    TNode *nodes;                //nodes array
    unordered_map<int, int> gid2leafid;    //mapping global vertex id to leaf node id
    long long tree_build_time;   // index buiding time

    TDGT() : nodes_num(0), nodes(nullptr), tree_build_time(DE_W) {}

    ~TDGT() { delete[] nodes; }

/**
 * synchronize borders generated in upper level to child nodes.
 * @param x node id
 */
    void sync_borders2children(int x) {
        for (int lgid : nodes[x].bid2lgid)
            nodes[x].partion_classes[lgid] = -nodes[x].partion_classes[lgid] - 1;
        vector<int> cnt_son_id(nodes[x].children.size(), 0);// replay lgid gene proc in child node
        for (int i = 0; i < nodes[x].G.n; i++) {
            if (nodes[x].partion_classes[i] < 0) {
                nodes[x].partion_classes[i] = -nodes[x].partion_classes[i] - 1;
                int cid = nodes[x].children[nodes[x].partion_classes[i]];
                nodes[cid].add_border(nodes[x].G.id[i], cnt_son_id[nodes[x].partion_classes[i]]);
                int gid = nodes[x].G.id[i];
                int bid = nodes[x].gid2bid[gid];
                assert(nodes[x].gid2bid.count(gid) && nodes[cid].gid2bid.count(gid));
                nodes[x].bid2cbid[bid] = nodes[cid].gid2bid[gid];// record bid in child
            }
            cnt_son_id[nodes[x].partion_classes[i]]++;
        }
    }

    void build_node(int x) {
        auto child_num = static_cast<int>(nodes[x].children.size());
        nodes[x].id = x;
        nodes[x].deep = nodes[nodes[x].father].deep + 1;
        if (nodes[x].G.n > LEAF_SIZE) { //non-leaf, partition it
            for (int i = 0; i < child_num; i++) {
                int cid = nodes_num + i + 1;
                nodes[x].children[i] = cid;
                nodes[cid].father = x;
            }
            nodes_num += child_num;
            auto **graphs = new TGraph *[child_num];
            for (int i = 0; i < child_num; i++)
                graphs[i] = &nodes[nodes[x].children[i]].G;
            nodes[x].G.Split(nodes[x].partion_classes, child_num);
            nodes[x].G.buildSubgraphs(nodes[x].partion_classes, graphs, child_num);
            for (auto cid:nodes[x].children)
                nodes[cid].initialize();
            delete[] graphs;
        } else { //leaf node, indexing all vertices, alternative choice for space efficiency
            // is indexing between only border and non-border vertices,
            for (int i = 0; i < nodes[x].G.n; i++) {
                nodes[x].partion_classes[i] = i;  // partition all vertices
                gid2leafid[nodes[x].G.id[i]] = x; //label each vertex its leaf node
            }
        }
        nodes[x].add_borders();
        nodes[x].stablize();
        nodes[x].init_matrix();
        if (nodes[x].father) {  //record bid in father
            int f = nodes[x].father;
            for (auto p:nodes[x].gid2bid)
                if (nodes[f].gid2bid.count(p.first))
                    nodes[x].bid2fbid[p.second] = nodes[f].gid2bid[p.first];
        }
        if (nodes[x].G.n > LEAF_SIZE)//sync gid2bid and mat edges to child
            sync_borders2children(x);
    }

/**
 * build Tree structure
 */
    void build_nodes() {
        std::ifstream in(tgraph_path);
        assert(in);
        unsigned int n;
        in >> n;
        in.close();
        //
        int n_ = int(n);
        cout << typeid(n_).name() << endl;
        cout << "num of nodes: " << n_ << endl;
        //
        nodes = new TNode[n];
        nodes[1].G.readGraph(tgraph_path);
        PRINT_BUILD("Building tree structure: vertices# " +
                    to_string(nodes[1].G.n) + "\t edges# " + to_string(nodes[1].G.m))
        nodes[1].initialize();
        nodes_num = 1;
        int x;
        queue<int> nodes_que;
        nodes_que.push(1);
        while (!nodes_que.empty()) {
            x = nodes_que.front();
            nodes_que.pop();
            build_node(x);
            if (nodes[x].G.n > LEAF_SIZE) {
                for (auto &child_id:nodes[x].children)
                    nodes_que.push(child_id);
            }
        }
        PRINT_BUILD("Finish Building Tree Structure with " + to_string(nodes_num) + " nodes")
    }

    //building travel time matrix, i.e., shortest travel time functions between border vertices.
    void build_distMatrix() {
        clock_t tbegin = clock(), tend;
        int disp_freq = 1000;
        for (int x = 1; x <= nodes_num; x++) //build from leaf to root, local optima
            nodes[x].add_edges2mat();// initialize matrix entries with connected edges

        //building buttom up phase
        for (int p = nodes_num; p > 1; p--) {
            nodes[p].floyd();
            nodes[p].sync_entries2father(nodes[nodes[p].father]);
            if (p % disp_freq == 0 || p < 20) {
                tend = clock();
                PRINT_BUILD("Building matrices up, " + node2tring(nodes[p]) +
                            "\ttime cost: " + to_string((tend - tbegin) / CLOCKS_PER_SEC) + "s")
                tbegin = clock();
            }
        }
        tbegin = clock();
        nodes[1].floyd();  //matrix for node[1] to optimal
        tend = clock();
        PRINT_BUILD("Building matrices up, " + node2tring(nodes[1]) +
                    "\ttime cost: " + to_string((tend - tbegin) / CLOCKS_PER_SEC) + "s")
        tbegin = clock();
        int cnt_down_updated_nodes = 0;
        for (int p = 1; p <= nodes_num; p++) {
            if (nodes[p].children[0]) {
                for (auto &pc:nodes[p].children) {
                    if (nodes[pc].sync_entries_from_father(nodes[p])) {// updated then recal
                        nodes[pc].floyd(true);
                        cnt_down_updated_nodes++;
                    }
                    if (pc % disp_freq == 0 || pc < 20) {
                        tend = clock();
                        PRINT_BUILD("Building matrices down, " + node2tring(nodes[pc]) +
                                    "\ttime cost: " + to_string((tend - tbegin) / CLOCKS_PER_SEC))
                        tbegin = clock();
                    }
                }
            }
        }
        string info = to_string(cnt_down_updated_nodes) + "/" + to_string(nodes_num) + "(" +
                      to_string((int) ((double) cnt_down_updated_nodes / nodes_num * 100)) +
                      "%) recomputed node matrices during down phase";
        PRINT_BUILD(info)
    }


    void buildTree() {
        PRINT_BUILD("Start building tree")
        clock_t tbegin = clock();
        build_nodes();
        PRINT_BUILD("Start building build_nodes")
        build_distMatrix();
        PRINT_BUILD("Start building build_distMatrix")
        clock_t tend = clock();
        PRINT_BUILD("Finish Building TDGT \t time cost: " +
                    to_string((tend - tbegin) / CLOCKS_PER_SEC) + "s")
    }


    double find_opt_dptime(PLF &travel_time, double &min_travel_time) {
        min_travel_time = DE_W;
        double t_depart = 0;
        for (auto &e:*travel_time.f) {
            double new_value = e.w;
            if (min_travel_time > new_value) {
                min_travel_time = new_value;
                t_depart = e.t;
            }
        }
        return t_depart;
    }

/**
 * push shortest travel time up from source to LCA
 * @param lfS  leaf node of source vertex
 * @param LCA  least common ancestor node
 * @param fwd_snode child node of LCA, to target node
 */
    void pushBordersUp(int lfS, int LCA, int fwd_snode) {
        //processing shortest path from lfS.father to LCA
        int p = lfS;
        for (; nodes[p].father != LCA; p = nodes[p].father) {
            nodes[nodes[p].father].push_up(nodes[p]);
        }
        //process shared b and record from to LCA, edge cond for up
        nodes[LCA].arr2b = std::vector<double>(nodes[LCA].gid2bid.size(), DE_W);
        nodes[LCA].pre_vex = std::vector<int>(nodes[LCA].gid2bid.size());
        int cnt1 = 0, cnt2 = 0;
        vector<int> from(nodes[p].bid2fbid.size()), to(nodes[fwd_snode].bid2fbid.size());
        for (auto &e:nodes[p].bid2fbid) {
            from[cnt1++] = e.second;
            nodes[LCA].arr2b[e.second] = nodes[p].arr2b[e.first];
            nodes[LCA].pre_vex[e.second] = e.second;
        }
        for (auto &e:nodes[fwd_snode].bid2fbid)
            if (nodes[LCA].arr2b[e.second] == DE_W)
                to[cnt2++] = e.second;
        //handle forwarding on LCA
        for (int i = 0; i < cnt1; i++) {
            for (int j = 0; j < cnt2; j++) {
                int b1 = from[i], b2 = to[j];
                double arr_new = nodes[LCA].matrix[b1][b2].dpt2arr(nodes[LCA].arr2b[b1]);
                if (lt(arr_new, nodes[LCA].arr2b[b2])) {
                    nodes[LCA].arr2b[b2] = arr_new;
                    nodes[LCA].pre_vex[b2] = b1;
                }
            }
        }
    }

    /**
     * push down from LCA to target
     * @param lfT leaf node of target vertex
     * @param LCA least common ancestor
     * @param down ordered tree node path from LCA to lfT
     * @param T_ target vertex
     */
    void pushBordersDown(int lfT, int LCA, stack<int> &down, int T_) {
        int pre = LCA;
        int p = down.top();
        down.pop();
        while (!down.empty()) {
            nodes[p].push_down(nodes[pre], nodes[down.top()]);
            pre = p;
            p = down.top();
            down.pop();
        }
        //process lfT based on arr to borders
        nodes[p].arr2b = std::vector<double>(nodes[p].gid2bid.size(), DE_W);
        nodes[p].pre_vex = std::vector<int>(nodes[p].gid2bid.size());
        for (auto &e:nodes[p].bid2fbid) {
            nodes[p].arr2b[e.first] = nodes[nodes[p].father].arr2b[e.second];
            nodes[p].pre_vex[e.first] = e.first;
            double arr_new = nodes[p].matrix[e.first][T_].dpt2arr(nodes[p].arr2b[e.first]);
            if (lt(arr_new, nodes[p].arr2b[T_])) {
                nodes[p].arr2b[T_] = arr_new;
                nodes[p].pre_vex[T_] = e.first;
            }
        }
    }


    int findLCA(int x, int y, stack<int> &down) {
        if (nodes[x].deep < nodes[y].deep)
            while (nodes[y].deep > nodes[x].deep) {
                down.push(y);
                y = nodes[y].father;
            }
        else
            while (nodes[x].deep > nodes[y].deep)
                x = nodes[x].father;
        while (x != y) {
            down.push(y);
            x = nodes[x].father;
            y = nodes[y].father;
        }
        return x;
    }

/**
 * time dependent shortest path
 * @param S source vertex
 * @param T target vertex
 * @param td departure time
 * @param need_path whether need path
 * @return
 */
    double tdsp(int S, int T, double td, bool need_path) {
        if (td >= TMAX) return 100000;
        if (S == T) return 0;
        int lfS = gid2leafid[S], lfT = gid2leafid[T];
        int S_ = nodes[lfS].gid2bid[S], T_ = nodes[lfT].gid2bid[T];
        nodes[lfS].matrix[S_][S_] = PLF({Segment(0, 0, INTV_CNTED)});
        nodes[lfT].matrix[T_][T_].f = nodes[lfS].matrix[S_][S_].f;
        int LCA;
        stack<int> down;
        nodes[lfS].arr2b = vector<double>(nodes[lfS].gid2bid.size());
        nodes[lfS].pre_vex = vector<int>(nodes[lfS].gid2bid.size());
        if (lfS == lfT) {//the same leaf node, calculate arr[S_, T_](t_d)
            LCA = lfS;
            PLF &cur_plf = nodes[lfS].matrix[S_][T_];
            nodes[lfS].arr2b[T_] = cur_plf.dpt2arr(td);
            nodes[lfS].pre_vex[T_] = S_;
            //assert(ge(nodes[lfS].arr2b[T_], td));
        } else {// in different nodes, push arr from lfS to lfT
            for (auto &e:nodes[lfS].bid2fbid) {
                PLF &cur_plf = nodes[lfS].matrix[S_][e.first];
                nodes[lfS].arr2b[e.first] = cur_plf.dpt2arr(td);
                nodes[lfS].pre_vex[e.first] = S_;
                //assert(ge(nodes[lfS].arr2b[e.first], td));
            }
            //get LCA and pushing down tree nodes path
            LCA = findLCA(lfS, lfT, down);
            //pushing arrs along passing tree nodes.
            pushBordersUp(lfS, LCA, down.top());
            pushBordersDown(lfT, LCA, down, T_);
        }
        if (need_path)
            path_recovery(S, T, td, LCA);
        return nodes[lfT].arr2b[T_] - td;
    }


/**
 * retrieve detailed path (recursive)
 * @param p  tree node id
 * @param S  retrieving source vertex
 * @param T  retrieving target vertex
 * @param path  path container
 * @param start_time  cur time
 */

    void find_path_recur(int p, int S, int T, vector<int> &path, double &start_time) {
        int S_ = nodes[p].gid2bid[S];
        int T_ = nodes[p].gid2bid[T];
        PLF &plf = nodes[p].matrix[S_][T_];
        auto p_seg = plf.dpt2seg(start_time);
//        assert(p_seg->intv == -1 || p_seg->nidx > 0);
        assert(p_seg->w != DE_W);
        if (p_seg->intv == INTV_CNTED) {//directly connected
            path.push_back(T);
            start_time = plf.dpt2arr(start_time, p_seg);//
        } else if (p_seg->intv < 0) {//connected in node -intv
            find_path_recur(-p_seg->intv, S, T, path, start_time);
        } else {// connected in cur node
            int K = nodes[p].bid2gid[p_seg->intv]; // intermediate vertex global id
            find_path_recur(p, S, K, path, start_time);
            return find_path_recur(p, K, T, path, start_time);
//            int k = nodes[p_seg->nidx].bid2gid[p_seg->intv]; //use global id
//            find_path_recur(p_seg->nidx, S, k, path, start_time);
//            return find_path_recur(p_seg->nidx, k, T, path, start_time);
        }
    }

    void path_recovery(int S, int T, double td, int LCA) {
        vector<int> path;
        vector<unsigned long> pieceidx;
        pieceidx.push_back(0);
        int lfS = gid2leafid[S], lfT = gid2leafid[T];
        if (lfS == lfT) {
            if (S != T) {
                find_path_recur(lfS, S, T, path, td);
                pieceidx.push_back(path.size());
            }
        } else {
            int q = lfT;
            int v2 = nodes[lfT].gid2bid[T];
            int v1 = nodes[q].pre_vex[v2];
            if (v1 != v2) {
                find_path_recur(q, nodes[q].bid2gid[v1], T, path, nodes[q].arr2b[v1]);
                pieceidx.push_back(path.size());
            }
            do {
                v2 = nodes[q].bid2fbid[v1]; // round to new turn in father node
                q = nodes[q].father;
                v1 = nodes[q].pre_vex[v2];
                if (v1 == v2) //shared bor v1=v2 in cur
                    continue;
                find_path_recur(q, nodes[q].bid2gid[v1], nodes[q].bid2gid[v2], path, nodes[q].arr2b[v1]);
                pieceidx.push_back(path.size());
            } while (q != LCA);
            stack<int> down;
            int p = gid2leafid[S];
            while (p != LCA) {
                down.push(p);
                p = nodes[p].father;
            }

            v2 = nodes[q].bid2cbid[v1]; //node pre LCA
            q = down.top();
            down.pop();
            while (q != lfS) {
                v1 = nodes[q].pre_vex[v2];
                if (v1 != v2) {
                    find_path_recur(q, nodes[q].bid2gid[v1], nodes[q].bid2gid[v2], path, nodes[q].arr2b[v1]);
                    pieceidx.push_back(path.size());
                }
                v2 = nodes[q].bid2cbid[v1]; //node pre LCA
                q = down.top();
                down.pop();
            }
            assert(down.empty());
            v1 = nodes[lfS].pre_vex[v2]; //lfS
            if (v1 != v2) {
                find_path_recur(lfS, S, nodes[lfS].bid2gid[v2], path, td);
                pieceidx.push_back(path.size());
            }
        }
        cout << S;
        for (auto iter = pieceidx.rbegin(); iter + 1 != pieceidx.rend(); iter++) {
            for (unsigned long i = *(iter + 1); i != *iter; i++)
                cout << " " << path[i];
        }
        cout << endl;
    }



/////////////////////////////////////////// functions to find PLC ///////////////////////////////////
/**
 * time interval shortest path
 * @param S source vertex
 * @param T target vertex
 * @param ts earlist dpt time of interval
 * @param te latest dpt time of interval
 */
    void PLCst(int S, int T, double ts, double te, PLF& plc) {
        int lfS = gid2leafid[S], lfT = gid2leafid[T];
        int S_ = nodes[lfS].gid2bid[S], T_ = nodes[lfT].gid2bid[T];
        nodes[lfS].matrix[S_][S_] = PLF({Segment(0, 0, INTV_CNTED)});
        nodes[lfT].matrix[T_][T_].f = nodes[lfS].matrix[S_][S_].f;
        int LCA;
        stack<int> down;
        nodes[lfS].arrf2b = vector<PLF>(nodes[lfS].gid2bid.size());
        nodes[lfS].pre_vex = vector<int>(nodes[lfS].gid2bid.size());
        if (lfS == lfT) {
            LCA = lfS;
            PLF &cur_plf = nodes[lfS].matrix[S_][T_];
            nodes[lfS].arrf2b[T_].f = cur_plf.getSlice(ts, te);
            nodes[lfS].arrf2b[T_].set_pre_vexs(S_);

            plc = cur_plf;
        } else {
            for (auto &e:nodes[lfS].bid2fbid) {
                PLF &cur_plf = nodes[lfS].matrix[S_][e.first];
                nodes[lfS].arrf2b[e.first].f = cur_plf.getSlice(ts, te);
                nodes[lfS].arrf2b[e.first].set_pre_vexs(S_);
            }
            LCA = findLCA(lfS, lfT, down);
            pushBordersUpPLC(lfS, LCA, ts, te, down.top());
            pushBordersDownPLC(lfT, LCA, ts, te, down, T_, plc);
            
        }
    }


//push travel time function up from source to LCA
    void pushBordersUpPLC(int lfS, int LCA, double ts, double te, int fwd_snode) {
        int p = lfS;
        for (; nodes[p].father != LCA; p = nodes[p].father)
            nodes[nodes[p].father].push_up(ts, te, nodes[p]);
        //process shared entries to LCA
        PLF plf_inf{Segment(ts, DE_W)};
        nodes[LCA].arrf2b = std::vector<PLF>(nodes[LCA].gid2bid.size(), plf_inf);
        int cnt1 = 0, cnt2 = 0;
        vector<int> from(nodes[p].bid2fbid.size()), to(nodes[fwd_snode].bid2fbid.size());
        for (auto &e:nodes[p].bid2fbid) {
            from[cnt1++] = e.second;
            nodes[LCA].arrf2b[e.second].f = nodes[p].arrf2b[e.first].f;
        }
        for (auto &e:nodes[fwd_snode].bid2fbid)
            if (nodes[LCA].arrf2b[e.second].f->front().w == DE_W)
                to[cnt2++] = e.second;
        for (int i = 0; i < cnt1; i++) //forwarding
            for (int j = 0; j < cnt2; j++) {
                int b1 = from[i], b2 = to[j];
                PLF PLFc;
                nodes[LCA].matrix[b1][b2].compound(ts, te, nodes[LCA].arrf2b[b1], PLFc, b1);
                nodes[LCA].arrf2b[b2].minimize(PLFc);
            }
    }

// push travel time function from LCA to target
    void pushBordersDownPLC(int lfT, int LCA, double ts, double te, stack<int> &down, int T_, PLF& plc) {
        int pre = LCA;
        int p = down.top();
        down.pop();
        while (!down.empty()) {//process based on tree node path from LCA to lfT
            nodes[p].push_down(ts, te, nodes[pre], nodes[down.top()]);
            pre = p;
            p = down.top();
            down.pop();
        }
        //process lfT based on arr to borders
        PLF plf_inf{Segment(ts, DE_W)};
        nodes[p].arrf2b = std::vector<PLF>(nodes[p].gid2bid.size(), plf_inf);
        for (auto &e:nodes[p].bid2fbid) {
            nodes[p].arrf2b[e.first].f = nodes[nodes[p].father].arrf2b[e.second].f;
            PLF plf_c;
            nodes[p].matrix[e.first][T_].compound(ts, te, nodes[p].arrf2b[e.first], plf_c, e.first);
            nodes[p].arrf2b[T_].minimize(plf_c);
        }
        plc = nodes[p].arrf2b[T_];
    }
/////////////////////////////////////////////////////////////////////////////////////////////
/**
 * time interval shortest path
 * @param S source vertex
 * @param T target vertex
 * @param ts earlist dpt time of interval
 * @param te latest dpt time of interval
 */
    void tisp(int S, int T, double ts, double te) {
        int lfS = gid2leafid[S], lfT = gid2leafid[T];
        int S_ = nodes[lfS].gid2bid[S], T_ = nodes[lfT].gid2bid[T];
        nodes[lfS].matrix[S_][S_] = PLF({Segment(0, 0, INTV_CNTED)});
        nodes[lfT].matrix[T_][T_].f = nodes[lfS].matrix[S_][S_].f;
        int LCA;
        stack<int> down;
        nodes[lfS].arrf2b = vector<PLF>(nodes[lfS].gid2bid.size());
        nodes[lfS].pre_vex = vector<int>(nodes[lfS].gid2bid.size());
        if (lfS == lfT) {
            LCA = lfS;
            PLF &cur_plf = nodes[lfS].matrix[S_][T_];
            nodes[lfS].arrf2b[T_].f = cur_plf.getSlice(ts, te);
            nodes[lfS].arrf2b[T_].set_pre_vexs(S_);
        } else {
            for (auto &e:nodes[lfS].bid2fbid) {
                PLF &cur_plf = nodes[lfS].matrix[S_][e.first];
                nodes[lfS].arrf2b[e.first].f = cur_plf.getSlice(ts, te);
                nodes[lfS].arrf2b[e.first].set_pre_vexs(S_);
            }
            LCA = findLCA(lfS, lfT, down);
            pushBordersUp(lfS, LCA, ts, te, down.top());
            pushBordersDown(lfT, LCA, ts, te, down, T_);
        }
//        find_opt_dptime()
//        if (need_path) {
//            path_recovery(S, T,)
//        }
    }


//push travel time function up from source to LCA
    void pushBordersUp(int lfS, int LCA, double ts, double te, int fwd_snode) {
        int p = lfS;
        for (; nodes[p].father != LCA; p = nodes[p].father)
            nodes[nodes[p].father].push_up(ts, te, nodes[p]);
        //process shared entries to LCA
        PLF plf_inf{Segment(ts, DE_W)};
        nodes[LCA].arrf2b = std::vector<PLF>(nodes[LCA].gid2bid.size(), plf_inf);
        int cnt1 = 0, cnt2 = 0;
        vector<int> from(nodes[p].bid2fbid.size()), to(nodes[fwd_snode].bid2fbid.size());
        for (auto &e:nodes[p].bid2fbid) {
            from[cnt1++] = e.second;
            nodes[LCA].arrf2b[e.second].f = nodes[p].arrf2b[e.first].f;
        }
        for (auto &e:nodes[fwd_snode].bid2fbid)
            if (nodes[LCA].arrf2b[e.second].f->front().w == DE_W)
                to[cnt2++] = e.second;
        for (int i = 0; i < cnt1; i++) //forwarding
            for (int j = 0; j < cnt2; j++) {
                int b1 = from[i], b2 = to[j];
                PLF PLFc;
                nodes[LCA].matrix[b1][b2].compound(ts, te, nodes[LCA].arrf2b[b1], PLFc, b1);
                nodes[LCA].arrf2b[b2].minimize(PLFc);
            }
    }

// push travel time function from LCA to target
    void pushBordersDown(int lfT, int LCA, double ts, double te, stack<int> &down, int T_) {
        int pre = LCA;
        int p = down.top();
        down.pop();
        while (!down.empty()) {//process based on tree node path from LCA to lfT
            nodes[p].push_down(ts, te, nodes[pre], nodes[down.top()]);
            pre = p;
            p = down.top();
            down.pop();
        }
        //process lfT based on arr to borders
        PLF plf_inf{Segment(ts, DE_W)};
        nodes[p].arrf2b = std::vector<PLF>(nodes[p].gid2bid.size(), plf_inf);
        for (auto &e:nodes[p].bid2fbid) {
            nodes[p].arrf2b[e.first].f = nodes[nodes[p].father].arrf2b[e.second].f;
            PLF plf_c;
            nodes[p].matrix[e.first][T_].compound(ts, te, nodes[p].arrf2b[e.first], plf_c, e.first);
            nodes[p].arrf2b[T_].minimize(plf_c);
        }
    }





};





// /**
//  * time dependent shortest path
//  * @param S source vertex
//  * @param T target vertex
//  * @param td departure time
//  * @param need_path whether need path
//  * @return
//  */
//     double tdsp_segments(int S, int T, double td, bool need_path, vector<tuple<double, double>>& segments) {
//         if (S == T) return 0;
//         int lfS = gid2leafid[S], lfT = gid2leafid[T];
//         int S_ = nodes[lfS].gid2bid[S], T_ = nodes[lfT].gid2bid[T];
//         nodes[lfS].matrix[S_][S_] = PLF({Segment(0, 0, INTV_CNTED)});
//         nodes[lfT].matrix[T_][T_].f = nodes[lfS].matrix[S_][S_].f;
//         int LCA;
//         stack<int> down;
//         nodes[lfS].arr2b = vector<double>(nodes[lfS].gid2bid.size());
//         nodes[lfS].pre_vex = vector<int>(nodes[lfS].gid2bid.size());
//         if (lfS == lfT) {//the same leaf node, calculate arr[S_, T_](t_d)
//             LCA = lfS;
//             PLF &cur_plf = nodes[lfS].matrix[S_][T_];

//             /// record the breakpoints
//             for (auto &e:cur_plf.f)
//             {
//                 double t = e.t;
//                 segments.push_back(make_tuple(e.t, e.w));
//                 // if(t > T){
//                 //     segments.push_back(make_tuple(e.t, e.w));
//                 // }
//             }
//             ///

//             nodes[lfS].arr2b[T_] = cur_plf.dpt2arr(td);
//             nodes[lfS].pre_vex[T_] = S_;
//             assert(ge(nodes[lfS].arr2b[T_], td));
//         } else {// in different nodes, push arr from lfS to lfT
//             for (auto &e:nodes[lfS].bid2fbid) {
//                 PLF &cur_plf = nodes[lfS].matrix[S_][e.first];

//                 /// record the breakpoints
//                 for (auto &e:cur_plf.f)
//                 {
//                     double t = e.t;
//                     segments.push_back(make_tuple(e.t, e.w));
//                     // if(t > T){
//                     //     segments.push_back(make_tuple(e.t, e.w));
//                     // }
//                 }                
//                 ///
//                 nodes[lfS].arr2b[e.first] = cur_plf.dpt2arr(td);
//                 nodes[lfS].pre_vex[e.first] = S_;
//                 assert(ge(nodes[lfS].arr2b[e.first], td));
//             }
//             //get LCA and pushing down tree nodes path
//             LCA = findLCA(lfS, lfT, down);
//             //pushing arrs along passing tree nodes.
//             pushBordersUp(lfS, LCA, down.top());
//             pushBordersDown(lfT, LCA, down, T_);
//         }
//         if (need_path)
//             path_recovery(S, T, td, LCA);
        
//         sort(segments.begin(), segments.end());
//         return nodes[lfT].arr2b[T_] - td;
//     }

#endif //TDGT_TDGT_H
