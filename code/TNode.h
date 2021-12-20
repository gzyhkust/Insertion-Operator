#ifndef TDGT_TNODE_H
#define TDGT_TNODE_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>
#include "PLF.h"
#include "Segment.h"
#include "Constants.h"
#include "TGraph.h"
#include "TMatrix.h"


using namespace std;

/**
 * Tree node structure
 */
struct TNode {
    int id, father, deep; //current node id，father node id
    std::vector<int> children;    //children nodes' ids
    TMatrix matrix;  //matrix of travel time functions(PLFs)，and associated intermediate vertices info.
    std::unordered_map<int, int> gid2bid; // global graph vertex id to border id
    TGraph G;  // associated sub-graph
    std::vector<int> bid2cbid, bid2gid; //border id in child node, global graph id
    //following  only used for init matrix
    std::vector<int> partion_classes;  //partitioning classes of vertices, i.e., 0,1,...FANOUT-1
    std::vector<int> bid2lgid;
    std::unordered_map<int, int> bid2fbid;

    std::vector<PLF> arrf2b;        //arrival time functions from source to borders for TISP querying
    std::vector<double> arr2b;      //arrival time from source to borders for TDSP querying
    std::vector<int> pre_vex;    //recording preceding vertices before borders for TDSP querying

    TNode() {
        id = 0;
        father = 0;
        deep = 0;
        children = std::vector<int>(FANOUT, 0);
    }

    void initialize() { //allocate enough space for containers
        assert(G.n > 0);
        partion_classes = std::vector<int>(G.n, -1);
        bid2cbid = std::vector<int>(G.n, -1);
        bid2gid = std::vector<int>(G.n, -1);
        bid2lgid = std::vector<int>(G.n, -1);
    }

    void stablize() {   // reclaim extra space
        bid2cbid.resize(gid2bid.size());
        bid2gid.resize(gid2bid.size());
        bid2lgid.resize(gid2bid.size());
    }

    void add_edges2mat() {     //initialize matrix entries with connected edges
        for (int bid1 = 0; bid1 < bid2lgid.size(); bid1++) {
            int lgid1 = bid2lgid[bid1];
            for (int j = G.head[lgid1]; j != -1; j = G.next[j]) {
                int lgid2 = G.adjv[j];
                if (gid2bid.count(G.id[lgid2])) { //another v is border
                    int bid2 = gid2bid[G.id[lgid2]];
                    matrix[bid1][bid2].f = G.weights[j].f; // modified to no copy
                }
            }
        }
    }

    void sync_entries2father(TNode &fnode) {
        for (auto &e1:bid2fbid)
            for (auto &e2: bid2fbid) {
                fnode.matrix[e1.second][e2.second] = matrix[e1.first][e2.first];
                fnode.matrix[e1.second][e2.second].set_pre_nodes(-id);
            }
    }

    bool sync_entries_from_father(TNode &fnode) {
        bool changed = false;
        for (auto &e1:bid2fbid)
            for (auto &e2: bid2fbid)
                if (e1.first != e2.first &&
                    matrix[e1.first][e2.first].min4sync(fnode.matrix[e1.second][e2.second], -fnode.id))
                    changed = true;
        return changed;
    }

    void add_border(int gid_global, int gid_local) {
        auto iter = gid2bid.find(gid_global);
        if (iter == gid2bid.end()) {
            int bid_new = static_cast<int>(gid2bid.size());
            gid2bid[gid_global] = bid_new;
            bid2gid[bid_new] = gid_global;
            bid2lgid[bid_new] = gid_local;
        }
    }

    void add_borders()                //add cut vertices as borders
    {
        for (int i = 0; i < G.n; i++) {
            int id = G.id[i];
            for (int j = G.head[i]; j != -1; j = G.next[j]) {
                if (partion_classes[i] != partion_classes[G.adjv[j]]) {
                    add_border(id, i);
                    break;
                }
            }
        }
    }

    void init_matrix() {// initialize matrix entries
        matrix.init(gid2bid.size());
    }


    // index shortest travel time functions
    void floyd(bool down = false) {
        for (int k = 0; k < gid2bid.size(); k++) {
            for (int i = 0; i < gid2bid.size(); i++) {
                for (int j = 0; j < gid2bid.size(); j++) {
                    if (down && bid2fbid.count(i) && bid2fbid.count(j))//sync from fnode, already optimal
                        continue;
                    if (i != j && j != k && i != k && neq(matrix[i][k].f->front().w, DE_W) &&
                        neq(matrix[k][j].f->front().w, DE_W)) {//relax current PLF
                        PLF PLFc; //compounded PLF
                        matrix[k][j].compound(matrix[i][k], PLFc, k);
                        assert(PLFc.f->front().w != DE_W);
                        matrix[i][j].minimize(PLFc);
                    }
                }
            }
        }
    }


    void push_up(TNode &pre_node) {
        arr2b = std::vector<double>(gid2bid.size(), DE_W);
        pre_vex = std::vector<int>(gid2bid.size());
        int cnt1 = 0, cnt2 = 0;
        vector<int> from(pre_node.bid2fbid.size()), to(bid2fbid.size());
        for (auto &e:pre_node.bid2fbid) { //init arr to shared b in cur
            from[cnt1++] = e.second;
            arr2b[e.second] = pre_node.arr2b[e.first];
            pre_vex[e.second] = e.second; //no pre in cur
        }
        for (auto &e:bid2fbid)
            if (arr2b[e.first] == DE_W)
                to[cnt2++] = e.first;
        for (int i = 0; i < cnt1; i++) {
            for (int j = 0; j < cnt2; j++) {
                int b1 = from[i], b2 = to[j];
                double arr_new = matrix[b1][b2].dpt2arr(arr2b[b1]);
                if (lt(arr_new, arr2b[b2])) {
                    arr2b[b2] = arr_new;
                    pre_vex[b2] = b1;
                }
            }
        }
    }

    void push_down(TNode &pre_node, TNode &fwd_node) {
        arr2b = std::vector<double>(gid2bid.size(), DE_W);
        pre_vex = std::vector<int>(gid2bid.size());
        int cnt1 = 0, cnt2 = 0;
        vector<int> from(bid2fbid.size()), to(fwd_node.bid2fbid.size());
        for (auto &e:bid2fbid) {
            from[cnt1++] = e.first;
            arr2b[e.first] = pre_node.arr2b[e.second];
            pre_vex[e.first] = e.first;
        }
        for (auto &e:fwd_node.bid2fbid)
            if (arr2b[e.second] == DE_W)
                to[cnt2++] = e.second;
        for (int i = 0; i < cnt1; i++) {
            for (int j = 0; j < cnt2; j++) {
                int b1 = from[i], b2 = to[j];
                double arr_new = matrix[b1][b2].dpt2arr(arr2b[b1]);
                if (lt(arr_new, arr2b[b2])) {
                    arr2b[b2] = arr_new;
                    pre_vex[b2] = b1;
                }
            }
        }
    }

    void push_up(double ts, double te, TNode &pre_node) {
        PLF plf_inf{Segment(ts, DE_W)};
        arrf2b = std::vector<PLF>(gid2bid.size(), plf_inf);
        int cnt1 = 0, cnt2 = 0;
        vector<int> from(pre_node.bid2fbid.size()), to(bid2fbid.size());
        for (auto &e:pre_node.bid2fbid) {//get shared borders from pre
            from[cnt1++] = e.second;
            arrf2b[e.second].f = pre_node.arrf2b[e.first].f;
        }
        for (auto &e:bid2fbid)
            if (arrf2b[e.first].f->front().w == DE_W)
                to[cnt2++] = e.first;
        for (int i = 0; i < cnt1; i++) //forwarding
            for (int j = 0; j < cnt2; j++) {
                int b1 = from[i], b2 = to[j];
                PLF PLFc;
                matrix[b1][b2].compound(ts, te, arrf2b[b1], PLFc, b1);
                arrf2b[b2].minimize(PLFc);
            }
    }

    void push_down(double ts, double te, TNode &pre_node, TNode &fwd_node) {
        PLF plf_inf{Segment(ts, DE_W)};
        arrf2b = std::vector<PLF>(gid2bid.size(), plf_inf);
        int cnt1 = 0, cnt2 = 0;
        vector<int> from(bid2fbid.size()), to(fwd_node.bid2fbid.size());
        for (auto &e:bid2fbid) {
            from[cnt1++] = e.first;
            arrf2b[e.first].f = pre_node.arrf2b[e.second].f;
        }
        for (auto &e:fwd_node.bid2fbid)
            if (arrf2b[e.second].f->front().w == DE_W)
                to[cnt2++] = e.second;
        for (int i = 0; i < cnt1; i++)
            for (int j = 0; j < cnt2; j++) {
                int b1 = from[i], b2 = to[j];
                PLF PLFc;
                matrix[b1][b2].compound(ts, te, arrf2b[b1], PLFc, b1);
                arrf2b[b2].minimize(PLFc);
            }
    }

    friend std::ostream &operator<<(std::ostream &out, TNode &);

};

std::ostream &operator<<(std::ostream &out, TNode &node) {
    if (node.children[0])
        out << "\tnleaf\t";
    else
        out << "\tleaf\t";
    out << "\tdeep: " << node.deep << "\t";
    out << "\tvexs num: " << node.G.n << "\t";
    out << "\tborders num: " << node.gid2bid.size() << "\t";
    out << "\tmat size: " << node.matrix.n * node.matrix.n << "\t";
    int cnt_segs = 0;
    for (int i = 0; i < node.matrix.n; i++)
        for (int j = 0; j < node.matrix.n; j++)
            cnt_segs += node.matrix[i][j].f->size();
    out << "segs num: " << cnt_segs;
    return out;
}

string node2tring(TNode &node) {
    string nodeinfo;
    nodeinfo += "node id:" + to_string(node.id) + "\t";
    if (node.children[0])
        nodeinfo += "nonleaf\t";
    else
        nodeinfo += "leaf\t";
    nodeinfo += "\tdeep: " + to_string(node.deep) + "\t";
    nodeinfo += "\tvexs num: " + to_string(node.G.n) + "\t";
    nodeinfo += "\tborders num: " + to_string(node.gid2bid.size()) + "\t";
    nodeinfo += "\tmat size: " + to_string(node.matrix.n * node.matrix.n) + "\t";
    int cnt_segs = 0;
    for (int i = 0; i < node.matrix.n; i++)
        for (int j = 0; j < node.matrix.n; j++)
            cnt_segs += node.matrix[i][j].f->size();
    nodeinfo += "segs num: " + to_string(cnt_segs);
    return nodeinfo;
}

#endif //TDGT_TNODE_H
