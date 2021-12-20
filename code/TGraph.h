#ifndef TDGT_TGRAPH_H
#define TDGT_TGRAPH_H

#include <memory>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
//#include <metis.h>
#include "/homes/zgongae/metis-5.1.0/include/metis.h"
#include "PLF.h"
#include "Segment.h"
#include "misc.h"

struct TGraph {
    unsigned long n, m;
    int mCnt; //count current added edges
    std::vector<int> id, head, next, adjv;  //id maintain global id of vertices
    std::vector<PLF> weights;   //time dependent edge weights

    TGraph() {
        n = 0;
        m = 0;
        mCnt = 0;
    }

    ~TGraph() {
        head.clear();
        id.clear();
        next.clear();
        adjv.clear();
        weights.clear();
    }

    void init(unsigned long _n, unsigned long _m) {
        this->n = _n;
        this->m = _m;
        head = std::vector<int>(n, -1);
        id = std::vector<int>(n, -1);
        next = std::vector<int>(m, -1);
        adjv = std::vector<int>(m, -1);
        weights = std::vector<PLF>(m);
    }

    // n,m    s,t weight piece num \n every piece by coordinate
    void readGraph(const std::string &path) {
        std::ifstream in(path);
        assert(in.is_open());
        unsigned int n, m, interNum, Period;
        in >> n >> m >> interNum >> Period;
        init(n, m);
        int vs, vt, weight_piece_num;
        while (in >> vs >> vt >> weight_piece_num) {//input edges
            auto f = std::make_shared<std::vector<Segment>>(weight_piece_num);
            double t, w;
            for (int i = 0; i < weight_piece_num; i++) {
                in >> t >> w;
                (*f)[i] = {t, w, INTV_CNTED};
            }
            addEdge(vs, vt, f);
        }
        for (int i = 0; i < n; i++)
            id[i] = i;
    }

    void addEdge(int s, int t, std::shared_ptr<std::vector<Segment>> &f) {
        adjv[mCnt] = t;
        next[mCnt] = head[s];
        head[s] = mCnt;
        weights[mCnt].f = f;
        assert(weights[mCnt].f->front().intv == INTV_CNTED);
        mCnt++;
    }

    //partition algorithm
    void Split(std::vector<int> &PClasses, int nparts = FANOUT)  //mark partition result in partion_classes
    {
        idx_t options[METIS_NOPTIONS];
        memset(options, 0, sizeof(options));
        {
            METIS_SetDefaultOptions(options);
            options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; // _RB
            options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // _VOL
            options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; // _RM  _SHEM
            options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_RANDOM; // _GROW _EDGE _NODE
            options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM; // _GREEDY _SEP2SIDED _SEP1SIDED
            options[METIS_OPTION_NCUTS] = 5;  // generate how many versions, choose the least cut one, 1,1
//            options[METIS_OPTION_NSEPS] = 10;
//            options[METIS_OPTION_NITER] = 10; //default
            /* balance factor, used to be 500 */
            options[METIS_OPTION_UFACTOR] = 600; //balance factor can be adjusted.
            // options[METIS_OPTION_MINCONN];
            options[METIS_OPTION_CONTIG] = 1;
            /////////options[METIS_OPTION_CONTIG] = 0;
            // options[METIS_OPTION_SEED];
            options[METIS_OPTION_NUMBERING] = 0;
            // options[METIS_OPTION_DBGLVL] = 0;
        }
        idx_t nvtxs = static_cast<idx_t>(n);
        idx_t ncon = 1;
        auto *xadj = new idx_t[n + 1];
        auto *adjncy = new idx_t[mCnt];
        //auto *adjwgt = new idx_t[mCnt]; changing by zengyang
        int *adjwgt = new idx_t[mCnt];
        auto *part = new idx_t[n];

        int xadj_pos = 1;
        int adjncy_pos = 0;

        xadj[0] = 0;
        for (int i = 0; i < n; i++) {
            for (int j = head[i]; j != -1; j = next[j]) {
                adjncy[adjncy_pos] = adjv[j];
                adjncy_pos++;
            }
            xadj[xadj_pos++] = adjncy_pos;
        }
        for (int i = 0; i < adjncy_pos; i++) // partition does not concern edge weight
            adjwgt[i] = 1;
        int objval = 0;

        ////////////////changing by zengyang
        int64_t l_nparts = nparts;
        int64_t l_objval = objval;
        ////////////////

        // METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, nullptr, nullptr, adjwgt, &l_nparts, nullptr, nullptr, options,
        //                     &l_objval,
        //                     part);
        cout << &nvtxs << " " << &ncon << " " << xadj << " " << adjncy << " " <<endl;
        METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, nullptr, nullptr, adjwgt, &nparts, nullptr, nullptr, options,
                    &objval,
                    part);
        
        // METIS_API(int) METIS_PartGraphKway(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, 
        //                 idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, 
        //                 idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, 
        //                 idx_t *edgecut, idx_t *part);
        for (int i = 0; i < n; i++) {
            PClasses[i] = part[i];
        }
        delete[] xadj;
        delete[] adjncy;
        delete[] adjwgt;
        delete[] part;


    }

    void buildSubgraphs(std::vector<int> &par_cls_idx, TGraph *Gs[], unsigned long cls_num = FANOUT) {
        // Partition
        std::vector<int> count1(cls_num, 0);
        std::vector<int> m(cls_num, 0);
        std::vector<int> new_id(n);
        //cnt vertices and edges for building subgraphs
        for (int i = 0; i < n; i++) {
            new_id[i] = count1[par_cls_idx[i]]++;
            for (int j = head[i]; j != -1; j = next[j]) {
                if (par_cls_idx[i] == par_cls_idx[adjv[j]])
                    m[par_cls_idx[i]]++;
            }
        }
        //add edges of subgraphs
        for (int t = 0; t < cls_num; t++) {
            (*Gs[t]).init(count1[t], m[t]);
            for (int i = 0; i < n; i++)
                if (par_cls_idx[i] == t)
                    for (int j = head[i]; j != -1; j = next[j])
                        if (par_cls_idx[i] == par_cls_idx[adjv[j]]) {
                            (*Gs[t]).addEdge(new_id[i], new_id[adjv[j]], weights[j].f);
                        }
        }
        //record global id for each vertex
        for (int i = 0; i < n; i++)
            (*Gs[par_cls_idx[i]]).id[new_id[i]] = id[i];
    }

};


#endif //TDGT_TGRAPH_H
