#ifndef TDGT_IO_H
#define TDGT_IO_H

#include <fstream>
#include <unordered_map>
#include "Segment.h"
#include "TDGT.h"

using namespace std;

template<typename T>
inline void write_value(std::ofstream &os, const T &value) {
    os.write((const char *) &value, sizeof value);
}

template<typename T>
inline void read_value(std::ifstream &is, T &value) {
    is.read((char *) &value, sizeof value);
}


//Segment
void write_value(std::ofstream &os, const Segment &seg) {
    assert(neq(seg.w, INT_MAX) && seg.intv != DE_INTV);
    write_value(os, seg.t);
    write_value(os, seg.w);
    write_value(os, seg.intv);
}

void read_value(std::ifstream &is, Segment &seg) {
    read_value(is, seg.t);
    read_value(is, seg.w);
    read_value(is, seg.intv);
    assert(neq(seg.w, DE_W) && seg.intv != DE_INTV);
}


//plf
void write_value(std::ofstream &os, const PLF &plf) {
    write_value(os, *plf.f);
}

void read_value(std::ifstream &is, PLF &plf) {
    assert(plf.f && plf.f->empty());
    read_value(is, *plf.f);
}


//vector
template<typename T>
void write_value(std::ofstream &os, const std::vector<T> &vec) {
    write_value(os, vec.size());
    for (auto &e:vec)
        write_value(os, e);
}

template<typename T>
void read_value(std::ifstream &is, std::vector<T> &vec) {
    typename std::vector<T>::size_type vsize;
    read_value(is, vsize);
    vec = std::vector<T>(vsize);
    for (auto &e:vec)
        read_value(is, e);
}

template<typename A, typename B>
void write_value(std::ofstream &os, const std::unordered_map<A, B> &m) {
    write_value(os, m.size());
    for (auto &e:m) {
        write_value(os, e.first);
        write_value(os, e.second);
    }
}

template<typename A, typename B>
void read_value(std::ifstream &is, std::unordered_map<A, B> &m) {
    typename std::unordered_map<A, B>::size_type msize;
    read_value(is, msize);
    for (int i = 0; i < msize; i++) {
        std::pair<A, B> e;
        read_value(is, e.first);
        read_value(is, e.second);
        m.insert(e);
    }
}

void write_matrix(std::ofstream &os, TMatrix &mat) {
    write_value(os, mat.n);
    for (int i = 0; i < mat.n; i++)
        for (int j = 0; j < mat.n; j++)
            if (i != j) {
                write_value(os, mat[i][j].f->size());
                for (auto &e:*mat[i][j].f)
                    write_value(os, e);
            }
}

void read_matrix(std::ifstream &is, TMatrix &mat) {
    read_value(is, mat.n);
    mat.init(mat.n);
    for (int i = 0; i < mat.n; i++)
        for (int j = 0; j < mat.n; j++)
            if (i != j) {
                vector<Segment>::size_type size;
                read_value(is, size);
                mat[i][j].f = make_shared<vector<Segment>>(size);
                for (auto &e:*mat[i][j].f)
                    read_value(is, e);
            }
}


void write_tree_node(std::ofstream &os, TNode &node) {
    write_value(os, node.id);
    write_value(os, node.father);
    write_value(os, node.deep);
    write_value(os, node.children);
    write_value(os, node.gid2bid);
    write_value(os, node.bid2fbid);
    write_value(os, node.bid2cbid);
    write_value(os, node.bid2gid);
    write_matrix(os, node.matrix);
}

void read_tree_node(std::ifstream &is, TNode &node) {
    read_value(is, node.id);
    read_value(is, node.father);
    read_value(is, node.deep);
    read_value(is, node.children);
    read_value(is, node.gid2bid);
    read_value(is, node.bid2fbid);
    read_value(is, node.bid2cbid);
    read_value(is, node.bid2gid);
    read_matrix(is, node.matrix);
}

void write_TDGT(std::ofstream &os, const TDGT &tree) {

    write_value(os, FANOUT);
    write_value(os, LEAF_SIZE);
    write_value(os, tree.nodes_num);
    write_value(os, tree.tree_build_time);
    for (int i = 1; i <= tree.nodes_num; i++)
        write_tree_node(os, tree.nodes[i]);
    write_value(os, tree.gid2leafid);
    os.close();
}

void read_TDGT(std::ifstream &is, TDGT &tree) {
    read_value(is, FANOUT);
    read_value(is, LEAF_SIZE);
    read_value(is, tree.nodes_num);
    read_value(is, tree.tree_build_time);
    tree.nodes = new TNode[tree.nodes_num + 1];
    for (int i = 1; i <= tree.nodes_num; i++)
        read_tree_node(is, tree.nodes[i]);
    read_value(is, tree.gid2leafid);
    is.close();
}

void print_shortcut(std::ostream &out, std::vector<Segment> &f) {
    for (auto &e:f)
        out << e << " ";
}


void print_matrix(std::ostream &out, TNode &tnode) {
    for (auto iteri = tnode.gid2bid.begin(); iteri != tnode.gid2bid.end(); iteri++) {
        for (auto iterj = tnode.gid2bid.begin(); iterj != tnode.gid2bid.end(); iterj++) {
            out << "Edge (" << iteri->first << "," << iterj->first << "): ";
            print_shortcut(out, *tnode.matrix[iteri->second][iterj->second].f);
            out << std::endl;
        }
    }
}

void print_graph(std::ostream &out, TGraph &G) {
    assert(G.head.size() == G.n);
    for (int i = 0; i < G.n; i++) {
        for (int j = G.head[i]; j != -1; j = G.next[j]) {
            int adj = G.adjv[j];
            out << "Edge (" << G.id[i] << "," << G.id[adj] << "): ";
            print_shortcut(out, *G.weights[j].f);
            out << std::endl;
        }
    }
}

template<typename T>
std::ostream &operator<<(std::ostream &out, std::vector<T> &vec) {
    for (auto &e:vec)
        out << e << ' ';
    return out;
}

#endif //TDGT_IO_H
