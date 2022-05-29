#ifndef TREE_DECOMP_H
#define TREE_DECOMP_H

#include <iostream>
#include <sstream>
#include <map>
#include <queue>
#include <unordered_set>
#include <vector>
#include <fstream>
#include <math.h>
#include <string>
#include <algorithm> 
#include "PLF.h"

using namespace std;

typedef unsigned int NodeId;
typedef unsigned int EdgeId;

struct vEdge {
    unsigned int dst;
    PLF weights;
};


struct TreeNode {
    int height;
    int pnodeid;
    vector<unsigned int> cnodeid;
    map<unsigned int, PLF> edges;

    map<unsigned int, PLF> rev_edges;

    map<unsigned int, PLF> labels;
};

vector<string> split (string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}
////////////////////////////////////////////////////////////////////////////


class TDTree
{
    public:
    int vorder = 0;
    NodeId n;
    EdgeId m;
    vector<vector<vEdge> > graph;
    vector<TreeNode> tnodes;
    vector<int> vertexOrder;
    vector<NodeId> sccid;
    //for 2-hop index:
    vector<NodeId> core_vertexes; 

    /////////////////////////////////////////////////// Debug Functions ///////////////////////////////////////////////////
    void report_shortcutes(vector<map<unsigned int, PLF> > &shortcuts)
    {
        cout << "Report Shortcurts ****************************: " << endl;
        for (int ii = 0; ii < graph.size(); ii++)
        {
            auto &vv = shortcuts[ii];
            cout << "From " << ii << " to: " << endl;
            for (map<unsigned int, PLF>::iterator it = vv.begin(); it!=vv.end(); ++it)
            {
                cout << it->first <<" ,";
                cout << it->second << endl;
            }
        }      
    }

    void report_degree2nodeQ(vector<vector<NodeId> > &degree2nodeQ)
    {
        cout << " Report degree2nodeQ ****************************: " << endl;
        for (int i = 0; i < degree2nodeQ.size(); i++)
        {
            cout << "degree: " << i <<":";
            for (int j = 0; j < degree2nodeQ[i].size(); j++)
            {
                cout << degree2nodeQ[i][j] << ",";
            }
            cout << endl;                
        }
    }

    void report_vPosition(vector<pair<NodeId, NodeId> > &vPosition)
    {
        cout << " Report vPosition ****************************: " << endl;
        for (int i = 0; i < graph.size(); i++)
        {
            cout << "vertex:" << i << " degree: " << vPosition[i].first << " index:" << vPosition[i].second << endl;
        }
        cout << endl;
    }    

    void report_vertexorder()
    {
        cout << "Report VertexOrder: " << endl;
        for (int ii = 0; ii < graph.size(); ii++)
        {
            cout << ii << ":" << vertexOrder[ii] << ", ";
        }
        cout << endl;
    }

    ///////////////////////////////////////////////////                 ///////////////////////////////////////////////////

    void load_graph_zy()
    {
        clock_t tbegin, tend;
        tbegin = clock();

        ifstream readgraph(tgrapgh_path);
        assert(readgraph.is_open());
        int a,b;
        readgraph >> n >> m >> a >> b;
        cout << "vertexes: " << n << " edges: " << m <<endl;
        graph.resize(n);

        int vs, vt, weight_piece_num;
        while (readgraph >> vs >> vt >> weight_piece_num) 
        {//input edges
            vEdge ve;
            ve.dst = vt;
            for (int i = 0; i < weight_piece_num; i++)
            {
                double t,w;
                readgraph >> t >> w;
                Segment seg(t,w);
                ve.weights.f->push_back(seg);
            }
            
            graph[vs].push_back(ve);
        }

        tend = clock();
        cout << "Finish loading graph \t time cost: " + to_string((tend - tbegin) / CLOCKS_PER_SEC) + "s" << endl;
        // cout << "============================report graph========================" << endl;
        // for (int i = 0; i < graph.size(); i++)
        // {
        //     for (int j = 0; j < graph[i].size(); j++)
        //     {
        //         auto &e = graph[i][j];
        //         cout << i << " " << e.dst << endl;
        //         cout << e.weights << endl;
        //         cout << " -------------------------------- " << endl;
        //     }
            
        // }
        // cout << "==========================report end======================" << endl;        
    }

    void report() 
    {
        cout << "============================report========================" << endl;
        cout << "---tree index---" << endl;

        for (unsigned int i = 0; i < tnodes.size(); i++)
        {
            cout << "vid: " << i << endl;
            cout << "order: " << vertexOrder[i] << endl;
            cout << "pnode: " << tnodes[i].pnodeid << endl;
            cout << "height: " << tnodes[i].height << endl;
            cout << "cnodes: ";
            for (auto c:tnodes[i].cnodeid)cout << c << ", ";
            cout << endl;
            cout << "X(node): ";
            for(map<unsigned int, PLF>::iterator it = tnodes[i].edges.begin(); it!=tnodes[i].edges.end(); ++it)
            {
                cout << it->first << ", ";
            }
            cout << endl;
            cout << " --------------------------------------" << endl;
        }
        

        // for (unsigned int i = 0; i < tnodes.size(); ++i) 
        // {
        //     cout << "vid(" << i << ") ord(" << vertexOrder[i] << ") pnode(" << tnodes[i].pnodeid << ") tid(" << sccid[i]
        //         << ") height("
        //         << tnodes[i].height << ") cnodes: ";
        //     for (auto c:tnodes[i].cnodeid)cout << c << ", ";
        //     cout << " edges: ";
        //     auto &edges = tnodes[i].edges;
        //     for (auto &e:edges) 
        //     {
        //         cout << " (" << e.first << "->" << e.second << ") ";
        //     }
        //     cout << endl;
        // }

        cout << "==========================report end======================" << endl;
    }    

    unsigned int report_idx_size() 
    {

        int intNum = 0, doubleNum = 0;
        intNum += vertexOrder.size();
        intNum += core_vertexes.size();

        for (int i = 0; i < graph.size(); i++)
        {
            for (int j = 0; j < graph[i].size(); j++)
            {
                auto &e = graph[i][j];
                intNum += 1;
                intNum += e.weights.f->size();
                doubleNum += 2*e.weights.f->size();
            }
            
        }

        for (int i = 0; i < tnodes.size(); i++)
        {
            intNum += 2;
            intNum += tnodes[i].cnodeid.size();

            for (map<unsigned int, PLF>::iterator it = tnodes[i].edges.begin(); it!=tnodes[i].edges.end(); ++it)
            {
                intNum+=1;
                intNum += it->second.f->size();
                doubleNum += it->second.f->size()*2;                
            }            
        }
        
        unsigned int idxsize = (intNum*4 + doubleNum*8)/1000000;
        cout << "idx_size= " <<idxsize << "MB" << endl;
        return idxsize;    
    }


    bool restore_idx()
    {
        return false;
    }

    void store_index()
    {
        ofstream index_write(index_path);
        for (unsigned int i = 0; i < tnodes.size(); i++)
        {
            index_write << i << " " << vertexOrder[i] << " " << tnodes[i].pnodeid << " " << tnodes[i].height << "\n";
            
            index_write << tnodes[i].cnodeid.size() << "\n";
            for (auto c:tnodes[i].cnodeid)index_write << c << "\n";
            //index_write << "\n";

            index_write  << tnodes[i].edges.size() << "\n";
            for (map<unsigned int, PLF>::iterator it = tnodes[i].edges.begin(); it!=tnodes[i].edges.end(); ++it)
            {
                int segsize = 0;
                for (auto f : *it->second.f) segsize ++;
                index_write << it->first << " " << segsize << "\n";
                for (auto f : *it->second.f)
                {
                    index_write << f.t << " " << f.w << " " << f.intv << "\n";
                }
            }

            index_write  << tnodes[i].rev_edges.size() << "\n";
            for (map<unsigned int, PLF>::iterator it = tnodes[i].rev_edges.begin(); it!=tnodes[i].rev_edges.end(); ++it)
            {
                int segsize = 0;
                for (auto f : *it->second.f) segsize ++;
                index_write << it->first << " " << segsize << "\n";
                for (auto f : *it->second.f)
                {
                    index_write << f.t << " " << f.w << " " << f.intv << "\n";
                }
            }            
        }
    }

    void restore_index()
    {
        ifstream index_read;
        index_read.open(index_path.c_str());
        if (!index_read.is_open()) {
            printf("%s does not exist\n", index_path.c_str());
            exit(0);
        }

        vertexOrder.resize(graph.size());
        tnodes.resize(graph.size());
        string line;
        string delimiter = " ";
        while(getline(index_read, line))
        {
            //cout << "vector<string> nodes = split(line, delimiter);" << endl;
            vector<string> nodes = split(line, delimiter);
            //cout << nodes[0] << " " << nodes[1] << " " << nodes[2] << " " << nodes[3]<<endl;
            //unsigned int vid = stoi(nodes[0]);
            int vid = stoi(nodes[0]);
            //cout << vid << endl;
            //cout << vertexOrder[vid] << endl;
            vertexOrder[vid] = stoi(nodes[1]);
            //cout << "vertexOrder[vid] = stoi(nodes[1]);" << endl;
            tnodes[vid].pnodeid = stoi(nodes[2]);
            //cout << "tnodes[vid].pnodeid = stoi(nodes[2]);" << endl;
            tnodes[vid].height = stoi(nodes[3]);
            //cout << "tnodes[vid].height = stoi(nodes[3]);" << endl;

            string cnodeidNumline;
            getline(index_read, cnodeidNumline);
            //cout << "vector<string> cnodeidNumstring = split(cnodeidNumline, delimiter);" << endl;
            vector<string> cnodeidNumstring = split(cnodeidNumline, delimiter);        
            for (int i = 0; i < stoi(cnodeidNumstring[0]); i++)
            {
                string cnodeidline;
                getline(index_read, cnodeidline);
                //cout << "vector<string> cnodeidstring = split(cnodeidline, delimiter);" << endl;
                vector<string> cnodeidstring = split(cnodeidline, delimiter);
                tnodes[vid].cnodeid.push_back(stoi(cnodeidstring[0]));            
            }

            string edgesNumline;
            getline(index_read, edgesNumline);
            //cout << "vector<string> edgesNumlinestring = split(edgesNumline, delimiter); " << endl;
            vector<string> edgesNumlinestring = split(edgesNumline, delimiter); 
            for (int i = 0; i < stoi(edgesNumlinestring[0]); i++)
            {
                string segsNumline;
                getline(index_read, segsNumline);
                //cout << "vector<string> segsNumlinetring = split(segsNumline, delimiter); " << endl;
                vector<string> segsNumlinetring = split(segsNumline, delimiter);
                int std =  stoi(segsNumlinetring[0]);
                for (int j = 0; j < stoi(segsNumlinetring[1]); j++)
                {
                    string weightline;
                    getline(index_read, weightline);
                    //cout << "vector<string> weightsstring = split(weightline, delimiter); " << endl;
                    vector<string> weightsstring = split(weightline, delimiter);
                    double t = stod(weightsstring.at(0));
                    double w = stod(weightsstring.at(1));
                    int intv = stoi(weightsstring.at(2));
                    Segment seg(t,w,intv);
                    tnodes[vid].edges[std].f->push_back(seg);
                }
                                
            }


            string rev_edgesNumline;
            getline(index_read, rev_edgesNumline);
            //cout << "vector<string> edgesNumlinestring = split(edgesNumline, delimiter); " << endl;
            vector<string> rev_edgesNumlinestring = split(rev_edgesNumline, delimiter); 
            for (int i = 0; i < stoi(rev_edgesNumlinestring[0]); i++)
            {
                string segsNumline;
                getline(index_read, segsNumline);
                //cout << "vector<string> segsNumlinetring = split(segsNumline, delimiter); " << endl;
                vector<string> segsNumlinetring = split(segsNumline, delimiter);
                int std =  stoi(segsNumlinetring[0]);
                for (int j = 0; j < stoi(segsNumlinetring[1]); j++)
                {
                    string weightline;
                    getline(index_read, weightline);
                    //cout << "vector<string> weightsstring = split(weightline, delimiter); " << endl;
                    vector<string> weightsstring = split(weightline, delimiter);
                    double t = stod(weightsstring.at(0));
                    double w = stod(weightsstring.at(1));
                    int intv = stoi(weightsstring.at(2));
                    Segment seg(t,w,intv);
                    tnodes[vid].rev_edges[std].f->push_back(seg);
                }
                                
            }                 
        }
    }

    void cal_parentchild_relationship_of_tree() 
    {
        for (unsigned int vid = 0; vid < tnodes.size(); ++vid) 
        {
            auto &node = tnodes[vid];
            auto &edges = node.edges;
            if (edges.size() > 0) 
            {   //not in the core or root
                NodeId neighbor = edges.begin()->first;
                NodeId idx = neighbor;
                int order = vertexOrder[neighbor];

                for (auto &e:edges) 
                {
                    neighbor = e.first;
                    if (vertexOrder[neighbor]!=-1) 
                    {
                        if (order==-1) 
                        {
                            idx = neighbor;
                            order = vertexOrder[neighbor];
                        } 
                        else if (vertexOrder[neighbor] < order) 
                        {
                            idx = neighbor;
                            order = vertexOrder[neighbor];
                        }
                    }
                }
                node.pnodeid = static_cast<int>(idx);
                tnodes[idx].cnodeid.push_back(vid);
            } 
            else 
            {
                node.pnodeid = -1;
            }
        }
    }

    void calculate_height()
    {
        queue<NodeId> Q;
        for (unsigned int i = 0; i < tnodes.size(); i++)
        {
            if (tnodes[i].pnodeid == -1)
            {
                Q.push(i);
                tnodes[i].height = 1;
            }   
        }

        while (!Q.empty()) 
        {
            auto v = Q.front();
            Q.pop();
            for (auto a:tnodes[v].cnodeid) 
            {
                Q.push(a);
                tnodes[a].height = tnodes[v].height + 1;
            }
        }
    }


    void graph_reduction_zy()
    {
        vector<vector<NodeId> > degree2nodeQ;//degree->a list of vertex of this degree
        vector<pair<NodeId, NodeId> > vPosition(graph.size());//(degree,idx)
        for (NodeId v = 0; v < graph.size(); ++v) 
        {
            NodeId degree = graph[v].size();
            if (degree >= degree2nodeQ.size()) {
                degree2nodeQ.resize(degree + 1);
            }
            vPosition[v] = make_pair(degree, degree2nodeQ[degree].size());
            degree2nodeQ[degree].push_back(v);
        }  


        vector<map<unsigned int, PLF> > shortcuts(graph.size());
        for (NodeId s = 0; s < graph.size(); s++)
        {
            for (unsigned int i = 0; i < graph[s].size(); ++i)
            {
                auto &e = graph[s][i];
                shortcuts[s][e.dst] = e.weights;
            }
        }  


        vertexOrder.resize(graph.size(), -1);
        NodeId mindegree = 0;
        tnodes.resize(graph.size());

        for (int i = 0; i < graph.size() - 1; i++)
        {
            mindegree = 0;
            while (mindegree < degree2nodeQ.size() && degree2nodeQ[mindegree].empty()) 
            {
                mindegree++;
            }

            NodeId vid = degree2nodeQ[mindegree].back();
            degree2nodeQ[mindegree].pop_back();
            vertexOrder[vid] = vorder++;  


            // reduing on node vid
            cout << "  ===================================== Reducing node =====================================: " << vid << endl;
            cout << "  ===================================== Mnidegree =====================================: " << mindegree << endl;
            auto &v = shortcuts[vid];

            // Reducing step 1: find nbr(vid, G_{i-1})
            vector<unsigned int> valid_neighbor_index;
            for (map<unsigned int, PLF>::iterator it = v.begin(); it!=v.end(); ++it)
            {
                // vertexOrder[id]!=-1 => vertex id has been eliminated 
                if (vertexOrder[it->first]==-1) 
                {
                    valid_neighbor_index.push_back(it->first);
                }  
            }

            // Reducing step 2: find u,w in nbr(vid, G_{i-1})
            vector<int> neighbor_degree_increase_cnt(valid_neighbor_index.size(), -1);
            for (unsigned int ii = 0; ii < valid_neighbor_index.size(); ++ii) 
            {
                for (unsigned int jj = ii + 1; jj < valid_neighbor_index.size(); ++jj) 
                {
                    NodeId ivid  = valid_neighbor_index[ii];
                    NodeId jvid  = valid_neighbor_index[jj];
                    //cout << "build shortcuts between: " << ivid << " " << jvid << endl;
                    
                    PLF PLFij, PLFji;
                    // shortcuts[ivid][vid].compound(shortcuts[vid][jvid],PLFij,vid);
                    // shortcuts[jvid][vid].compound(shortcuts[vid][ivid],PLFji,vid);
                    shortcuts[vid][jvid].compound(shortcuts[ivid][vid],PLFij,vid);
                    shortcuts[jvid][vid].compound(shortcuts[vid][ivid],PLFji,vid);
                    if (shortcuts[ivid].find(jvid) == shortcuts[ivid].end())
                    {
                        neighbor_degree_increase_cnt[ii]++;
                        neighbor_degree_increase_cnt[jj]++;
                        shortcuts[ivid][jvid] = PLFij;
                        shortcuts[jvid][ivid] = PLFji;
                    }
                    else
                    {
                        shortcuts[ivid][jvid].minimize(PLFij);
                        shortcuts[jvid][ivid].minimize(PLFji);
                    }



                    // map<unsigned int, PLF>::iterator it = shortcuts[ivid].begin();
                    // while (it!=shortcuts[ivid].end())
                    // {
                    //     if (it->first == jvid)
                    //     {
                    //         break;
                    //     }
                    //     it++;
                    // }
                    // if ((it == shortcuts[ivid].end()) and (it->first != jvid))
                    // {
                    //     neighbor_degree_increase_cnt[ii]++;
                    //     neighbor_degree_increase_cnt[jj]++;

                    //     //cout << ivid << " not connect to " << jvid << endl;
                    //     //cout << ivid << " :" << neighbor_degree_increase_cnt[ii] << "; " << jvid << " :" << neighbor_degree_increase_cnt[jj] <<endl;   
                    // }
                    // PLF PLFij, PLFji;
                    // shortcuts[ivid][vid].compound(shortcuts[vid][jvid],PLFij,vid);
                    // shortcuts[ivid][jvid] = PLFij;
                    // shortcuts[jvid][vid].compound(shortcuts[vid][ivid],PLFji,vid);
                    // shortcuts[jvid][ivid] = PLFji;
                    
                    // cout << shortcuts[ivid].find(jvid) << endl;
                    // cout << shortcuts[jvid].end() << endl;
                    // if (shortcuts[ivid].find(jvid) == shortcuts[jvid].end())
                    // {
                    //     neighbor_degree_increase_cnt[ii]++;
                    //     neighbor_degree_increase_cnt[jj]++;

                    //     cout << ivid << " not connect to " << jvid << endl;
                    //     cout << ivid << " :" << valid_neighbor_index[ii] << "; " << jvid << " :" << valid_neighbor_index[jj] <<endl;
                    // }
                }

            }
            
            // report_shortcutes(shortcuts);
            // report_vertexorder();


            // Reducing step 3: update vPosition of compound vertexes
            for (unsigned int i = 0; i < valid_neighbor_index.size(); ++i) 
            { 
                //if (neighbor_degree_increase_cnt[i] > 0) 
                if (neighbor_degree_increase_cnt[i] != 0) 
                {
                    //cout << "1. swap and delete " << endl;
                    int x = valid_neighbor_index[i];
                    //cout << " pair<NodeId, NodeId> &p = vPosition[x]; " << vPosition[x].first << " " << vPosition[x].second << endl;
                    pair<NodeId, NodeId> &p = vPosition[x];
                    //1. swap and delete
                    //cout << " degree2nodeQ[p.first][p.second] = degree2nodeQ[p.first].back(); " << degree2nodeQ[p.first].size() << " " << p.second << endl;
                    degree2nodeQ[p.first][p.second] = degree2nodeQ[p.first].back();
                    //cout << " vPosition[degree2nodeQ[p.first].back()].second = p.second; " << degree2nodeQ[p.first].back() << " " << p.second << endl;
                    vPosition[degree2nodeQ[p.first].back()].second = p.second;
                    //cout << " degree2nodeQ[p.first].pop_back(); " << p.first << " " << p.second << endl;
                    degree2nodeQ[p.first].pop_back();
                    //2. place in a new position
                    //cout << "place in a new position" << endl;
                    //p.first += static_cast<unsigned int>(neighbor_degree_increase_cnt[i]);
                    //cout << "2. place in a new position " << endl;
                    //cout << " test :  p.first += neighbor_degree_increase_cnt[i];"<< p.first << " " << neighbor_degree_increase_cnt[i] << endl;
                    p.first += neighbor_degree_increase_cnt[i];
                    if (p.first == 0)
                    {
                        continue;
                    }
                    
                    //cout << " after :  p.first += neighbor_degree_increase_cnt[i];"<< p.first << endl;
                    if (p.first >= degree2nodeQ.size()) 
                    {
                        degree2nodeQ.resize(p.first + 1);
                    }
                    //cout << "mindegree = min(mindegree, p.first);" << endl;
                    //mindegree = min(mindegree, p.first);
                    //cout << "p.second = degree2nodeQ[p.first].size();" << endl;
                    //cout << "p.first: " << p.first << " " << "mindegree: " << mindegree << endl;
                    //cout << "degree2nodeQ.size(): " << degree2nodeQ.size() << endl;
                    p.second = degree2nodeQ[p.first].size();
                    //cout << " degree2nodeQ[p.first].push_back(x);" << endl;
                    degree2nodeQ[p.first].push_back(x);
                }
                //cout << " ------------------ finish neighbor: ------------ " << i << " " << valid_neighbor_index.size() << endl;
            }

            // report_degree2nodeQ(degree2nodeQ);
            // report_vPosition(vPosition);

            // Reducing step 4: building tnodes

            //cout << "Reducing step 4: building tnodes" << endl;
            for(auto i:valid_neighbor_index) 
            {
                //cout << i << " " << v.size() <<endl;
                //cout << v[i] << endl;
                //cout << tnodes.size() << endl;
                tnodes[vid].edges[i] = v[i];
                //cout << tnodes[vid].edges[i] << " &&&&&&&&&&&&&&&& " << endl;
                //tnodes[vid].edges.insert(std::move(v.idx(i)));
                // reverse
                // tnodes[i].rev_edges[vid] = v[i];
            }
            if (tnodes[vid].edges.size()==0) {
                //cout << "if (tnodes[vid].edges.size()==0)" << endl; 
                core_vertexes.emplace_back(vid);
            }          
        }
        
        //cout << "cal_parentchild_relationship_of_tree();" << endl;
        cal_parentchild_relationship_of_tree();
        //cout << "calculate_height();" << endl;
        calculate_height();   

        // cout <<"check reverse ......." << endl;
        // for (unsigned int i = 0; i < tnodes.size(); i++)
        // {
        //     cout << tnodes[i].rev_edges.size() << endl;
        // }
         

    }

    pair<unsigned int, unsigned int> qLCA(NodeId s, NodeId t) 
    {
        NodeId sAncestor = s, tAncestor = t;
        while (true) 
        {
            if (tnodes[sAncestor].height > tnodes[tAncestor].height) 
            {
                if (static_cast<int>(tAncestor)==tnodes[sAncestor].pnodeid) return make_pair(tAncestor, tAncestor);
                sAncestor = static_cast<unsigned int>(tnodes[sAncestor].pnodeid);
            } 
            else if (tnodes[sAncestor].height < tnodes[tAncestor].height) 
            {
                if (static_cast<int>(sAncestor)==tnodes[tAncestor].pnodeid) return make_pair(sAncestor, sAncestor);
                tAncestor = static_cast<unsigned int>(tnodes[tAncestor].pnodeid);
            } 
            else 
            {
                if (tnodes[sAncestor].pnodeid==-1 || sAncestor==tAncestor)return make_pair(sAncestor, tAncestor);
                sAncestor = static_cast<unsigned int>(tnodes[sAncestor].pnodeid);
                tAncestor = static_cast<unsigned int>(tnodes[tAncestor].pnodeid);
            }
        }
    }

    int exist_in_tnode(NodeId id, NodeId i)
    {
        map<unsigned int, PLF>::iterator it;
        it = tnodes[i].edges.find(id);
        if(it == tnodes[i].edges.end())
        {
            return 0;
        }
        else
        {
            return 1;
        }
        
    }

    vector<unsigned int> uVertexes(NodeId s, NodeId p)
    {
        vector<NodeId> vertexes;
        for (map<unsigned int, PLF>::iterator it = tnodes[p].edges.begin(); it!=tnodes[p].edges.end(); it++)
        {
            if (exist_in_tnode(it->first, s) == 0)
            {
                vertexes.push_back(it->first);
            }
        }
        return vertexes;  
    }

    vector<NodeId> vVetexes(NodeId s, NodeId p)
    {
        vector<NodeId> vertexes;
        for (map<unsigned int, PLF>::iterator it = tnodes[p].edges.begin(); it!=tnodes[p].edges.end(); it++)
        {
            if (exist_in_tnode(it->first, s) == 1)
            {
                vertexes.push_back(it->first);
            }
        }
        vertexes.push_back(p);
        return vertexes;  
    }

    map<NodeId, PLF> Path_computing_bottom2up_zy(NodeId bottom, NodeId up)
    {
        map<unsigned int, PLF> results;
        NodeId s = bottom;
        if (bottom != up)
        {
            map<unsigned int, PLF> dismap = tnodes[bottom].edges;
            NodeId anc = tnodes[bottom].pnodeid;
            while (anc!=up)
            {
                vector<NodeId> uVertex = uVertexes(s,anc);
                vector<NodeId> vVetex = vVetexes(s, anc);

                // cout << " ********************************************* " << endl;
                // cout << "check node: " << anc << endl;           

                for (auto &uid : uVertex)
                {
                    //cout << "new vertex " << uid << endl;
                    for (auto &vid : vVetex)
                    {
                        PLF plfvu;
                        if (exist_in_tnode(uid, vid) == 0)
                        {
                            continue;
                        }
                        
                        // cout << "compound " << vid << "->" << uid << endl;
                        // cout << tnodes[vid].edges[uid] << " compound to :" << dismap[vid] << endl;
                        tnodes[vid].edges[uid].compound(dismap[vid],plfvu,vid);
                    
                        map<unsigned int, PLF>::iterator it;
                        it = dismap.find(uid);
                        if(it == dismap.end())
                        {
                            // cout << "first compound: " << endl;
                            // cout << plfvu << endl;
                            dismap[uid] = plfvu;
                        } 
                        else
                        {
                            // cout << "anpther compound: " << endl;
                            // cout << dismap[uid] << endl;
                            // cout << plfvu << endl;
                            dismap[uid].minimize(plfvu);

                        }
                    }
                    // cout << "dismap[uid] ->" << uid << endl;
                    // cout << dismap[uid] << endl;
                    // cout << " --------------------------------- " << endl;
                }

                s = anc;
                anc = tnodes[anc].pnodeid;
            }

            for (map<unsigned int, PLF>::iterator it = dismap.begin(); it!=dismap.end(); ++it)
            {
                results[it->first] = it->second;
            }
        }
        else
        {
            for (map<unsigned int, PLF>::iterator it = tnodes[up].edges.begin(); it!=tnodes[up].edges.end(); ++it)
            {
                results[it->first] = it->second;
            } 
        }

        return results;
    }

    //vector<PLF> query(NodeId s, NodeId t) 
    PLF query(NodeId s, NodeId t) 
    {
        vector<PLF> result;
        auto lca = qLCA(s, t);
        //cout << lca.first << " " << lca.second << endl;
        //auto s_result = Path_computing_bottom2up(s, lca.first);
        auto s_result = Path_computing_bottom2up_zy(s, lca.first);
        auto t_result = Path_computing_bottom2up_zy(t, lca.second);
        //auto t_result = Path_computing_bottom2up(t, lca.second);
        // cout << "s_result : " << s_result.size() <<endl;
        // for (map<unsigned int, PLF>::iterator it = s_result.begin(); it!=s_result.end(); ++it)
        // {
        //     cout << it->first <<" ,";
        //     cout << it->second << endl;
        // }
        // cout << "t_result : " << endl;
        // for (map<unsigned int, PLF>::iterator it = t_result.begin(); it!=t_result.end(); ++it)
        // {
        //     cout << it->first <<" ,";
        //     cout << it->second << endl;
        // }
        
        PLF result_;
        if (lca.first == t)
        {
            //result.push_back();
            result_ = s_result[t];
        }
        else if (lca.first == s)
        {
            //result.push_back();
            result_ = s_result[s];
        }
        else
        {
            PLF plfvu;
            t_result[lca.first].compound(s_result[lca.first],plfvu,lca.first);
            result_ = plfvu;
        }
        return result_;
        
        
        
        // for (map<NodeId, PLF>::iterator its = s_result.begin(); its!=s_result.end(); ++its)
        // {
        //     cout << its->first << endl;
        //     map<NodeId, PLF>::iterator itt = t_result.find(its->first);
        //     if (itt!=t_result.end())
        //     {
        //         PLF plfvu;
        //         s_result[its->first].compound(t_result[its->first],plfvu,its->first);
        //         result.push_back(plfvu);                
        //     }
            
        // }
        // return result;
    }


    void creat_index() 
    {

        //graph_reduction();
        graph_reduction_zy();

        //report();
    }

    void Init() 
    {
        load_graph_zy();

        vector<int> dgr;

        // for (int i = 0; i < graph.size(); i++)
        // {
        //     dgr.push_back(graph[i].size());
        // }
        // auto maxPosition = max_element(dgr.begin(), dgr.end());
        // auto minPosition = min_element(dgr.begin(), dgr.end());
        // cout << *maxPosition << " at the position of " << maxPosition - dgr.begin() << endl;
        // cout << *minPosition << " at the position of " << minPosition - dgr.begin() << endl;

        clock_t tbegin = clock();
        if (!restore_idx()) 
        {
            // creat_index();
            // store_index();
            
            restore_index();
        }
        unsigned int idxsiz = report_idx_size();
        clock_t tend = clock();
        cout << "Finish building time cost: " + to_string((tend - tbegin) / CLOCKS_PER_SEC) + "s" << endl;

        clock_t querybegin = clock();
        int num = 0;        
        for (int x = 0; x < 1000; x++)
        {
            int i = rand()%400000;
            int j = rand()%400000;
            if (i == j)
            {
                j++;
            }
            //cout << i << " " << j << endl;
            PLF result = query(i,j);       
        }
        
        clock_t queryend = clock();
        cout << "Finish querying time cost: " + to_string((queryend - querybegin) / CLOCKS_PER_SEC) + "s" << endl;                
        

        // ofstream of("./result.txt");
        // of << (tend - tbegin) / CLOCKS_PER_SEC << "\n"; 

        // ************************ debug query
        // cout << "============================= check shorcuts================================" << endl;
        // cout << "7->5" << endl;
        // cout << tnodes[7].edges[5] << endl;

        // cout << " ---------------------------- " << endl;
        // //vector<PLF> result1 = query(7,8);
        // cout << 7 << " " << 2 << endl;
        // vector<PLF> result1 = query(7,2);

        // cout << " ****************************** " << endl;
        // cout << 7 << " " << 0 << endl;
        // //vector<PLF> result2 = query(7,0);
        // PLF result2 = query(7,0);

        // cout << " ****************************** " << endl;
        // cout << 7 << " " << 8 << endl;
        // //vector<PLF> result3 = query(7,8);        
        // PLF result3 = query(7,8);
        // cout << result3 << endl;     

        // for (int i = 0; i < 9; i++)
        // {
        //     for (int j = 0; j < 9; j++)
        //     {
        //         if(i == j)
        //         {
        //             continue;
        //         }
        //         cout << i << " " << j << endl;
        //         cout << query(i,j) << endl;
        //         cout << "--------------------------------------" << endl;
        //     }
            
        // }
        

        // cout << " ****************************** " << endl;
        // cout << 7 << " " << 4 << endl;
        // vector<PLF> result3 = query(7,4);
        //cout << tnodes[7].edges[5] << endl;
        // ************************



        // vector<int> cnodeNum;
        // vector<int> pNum;
        // int avg_cnodenum=0, avgxnodenum=0;
        // int maxheight = -100;
        // for (int i = 0; i < tnodes.size(); i++)
        // {
        //     //cout << tnodes[i].height << endl;
        //     if(tnodes[i].height > maxheight)
        //     {
        //         maxheight = tnodes[i].height;
        //     }

        //     cnodeNum.push_back(tnodes[i].cnodeid.size());

        //     int psize = 0;
        //     for (map<unsigned int, PLF>::iterator it = tnodes[i].edges.begin(); it!=tnodes[i].edges.end(); ++it)
        //     {
        //         psize++;
        //     }
            
        //     pNum.push_back(psize);    

        //     if (!tnodes[i].cnodeid.empty())
        //     {
        //         avg_cnodenum += tnodes[i].cnodeid.size();
        //         //cout << tnodes[i].cnodeid.size() << endl;
        //     }
        //     avgxnodenum += tnodes[i].edges.size();    
        // }
        
        // cout << "maxheight: " << maxheight << endl;
        // auto maxcnodeNum = max_element(cnodeNum.begin(), cnodeNum.end());
        // auto mincnodeNum = min_element(cnodeNum.begin(), cnodeNum.end()); 
        // cout << "conde: " << *maxcnodeNum << " " << *mincnodeNum << endl;       

        // auto maxpNum = max_element(pNum.begin(), pNum.end());
        // auto mincpNum = min_element(pNum.begin(), pNum.end());
        // cout << "pNum: " << *maxpNum << " " << *mincpNum << endl;

        // cout << " avg_cnodenum " << avg_cnodenum/tnodes.size() << "  " << " avgxnodenum " << avgxnodenum/tnodes.size() <<endl;

        // ofstream of("results.txt");
        // of << (tend - tbegin) / CLOCKS_PER_SEC << "\n";
        // of << idxsiz << "\n";


        //report();

        // vector<PLF> result1 = query(7,8);

        // cout << endl;
        // cout << endl;

        // vector<PLF> result2 = query(2,3);
        // clock_t querybegin = clock();
        // int num = 0;        
        // for (int i = 1; i < 10; i++)
        // {
        //     for (int j = 40000; j < 40010; j++)
        //     {
        //         cout << i << " " << j << endl;
        //         vector<PLF> result = query(i,j);
        //         //auto lca = qLCA(i, j);
        //         for(auto &plf:result)
        //         {
        //             cout << plf.dpt2arr(10) << " ";
        //         }
        //         cout << "-----------------------------" << endl;
        //     }
        // }
        // for (int x = 0; x < 100; x++)
        // {
        //     int i = rand()%400000;
        //     int j = rand()%400000;
        //     if (i == j)
        //     {
        //         j++;
        //     }
        //     cout << i << " " << j << endl;
        //     vector<PLF> result = query(i,j);
        //     for(auto &plf:result)
        //     {
        //         plf.dpt2arr(10);
        //     }           
        // }
        
        // clock_t queryend = clock();
        // cout << "Finish querying time cost: " + to_string((queryend - querybegin) / CLOCKS_PER_SEC) + "s" << endl;                
        

        // int Xnum = 0;
        // for (int i = 0; i < tnodes.size(); i++)
        // {
        //     Xnum+=tnodes[i].edges.size();
        // }
        // cout << Xnum << " " << tnodes.size() << endl;    
        











        // /// test 1
        // clock_t querybegin = clock();
        // int num = 0;
        // for (int i = 100; i < 200; i++)
        // {
        //     for (int j = 20000; j < 20100; j++)
        //     {
        //         num ++;
        //         cout << i << " " << j << endl;
        //         vector<PLF> result = query(i,j);
        //         //auto lca = qLCA(i, j);
        //     }
            
        // }
        // clock_t queryend = clock();
        // cout << "Finish querying time cost: " + to_string((queryend - querybegin) / CLOCKS_PER_SEC) + "s" << endl;
        // cout << num << endl;

        // /// test 2
        // querybegin = clock();
        // num = 0;
        // for (int i = 1000; i < 1001; i++)
        // {
        //     for (int j = 10000; j < 10004; j++)
        //     {
        //         num ++;
        //         cout << i << " " << j << endl;
        //         vector<PLF> result = query(i,j);
        //     }
            
        // }
        // queryend = clock();
        // cout << "Finish querying time cost: " + to_string((queryend - querybegin) / CLOCKS_PER_SEC) + "s" << endl;
        // cout << num << endl;

        // /// test 3
        // querybegin = clock();
        // num = 0;
        // for (int i = 2000; i < 2100; i++)
        // {
        //     for (int j = 400000; j < 400003; j++)
        //     {
        //         num ++;
        //         cout << i << " " << j << endl;
        //         vector<PLF> result = query(i,j);
        //     }
            
        // }
        // queryend = clock();
        // cout << "Finish querying time cost: " + to_string((queryend - querybegin) / CLOCKS_PER_SEC) + "s" << endl;
        // cout << num << endl;        

    }

};




#endif //TREE_DECOMP_H


// void graph_reduction() 
// {
//     // STEP1: build order
//     clock_t tbegin, tend;
//     tbegin = clock();
//     vector<unsigned int> order2id;
//     vector<vector<unsigned int> > deg2vid;
//     for (uint i = 0; i < graph.size(); ++i) 
//     {
//         if (deg2vid.size() <= graph[i].size()) 
//         {
//             deg2vid.resize(graph[i].size() + 1);
//         }
//         deg2vid[graph[i].size()].push_back(i);
//     }
//     vertexOrder.resize(graph.size());
//     for (auto &d:deg2vid) 
//     {
//         for (auto vid:d) 
//         {
//             order2id.push_back(vid);
//             vertexOrder[vid] = vorder++;
//         }
//     }

//     // STEP2: edge shortcut in original graph
//     vector<map<unsigned int, PLF> > shortcuts(graph.size());
//     for (NodeId s = 0; s < graph.size(); s++)
//     {
//         for (unsigned int i = 0; i < graph[s].size(); ++i)
//         {
//             auto &e = graph[s][i];
//             if (vertexOrder[s] < vertexOrder[e.dst])
//             {
//                 shortcuts[s][e.dst].push_back(e.weights);
//             }
            
//         }
        
//     }
    
//     // STEP3: build tnodes and add some shortcuts
//     tnodes.resize(graph.size());

//     NodeId reduced_cnt = 0;

//     // enumerate each vertex to try to reduce based on the order
//     for (unsigned int order = 0; order < order2id.size(); order++)
//     {
//         reduced_cnt ++;
//         NodeId vid = order2id[order];
//         auto &e = shortcuts[vid];

//         for (int i = 0; i < e.size(); i++)
//         {
//             for (int j = i+1; j < e.size(); j++)
//             {
                
//             }
            
//         }
        
//     }

//     for (unsigned int order = 0; order < order2id.size(); ++order) 
//     {
//         reduced_cnt++;
//         NodeId vid = order2id[order];
//         auto &v = shortcuts[vid];

//         unsigned int i, j;
//         //add short cuts
//         for (i = v._begin(); i!=v._end(); v._next(i)) 
//         {
//             for (j = i, v._next(j); j!=v._end(); v._next(j)) 
//             {
//                 auto &ipair = v.idx(i);
//                 auto &jpair = v.idx(j);
//                 if (vertexOrder[ipair.first] > vertexOrder[jpair.first]) 
//                 {
//                     shortcuts[jpair.first][ipair.first].combine(ipair.second + jpair.second);
//                 } else 
//                 {
//                     shortcuts[ipair.first][jpair.first].combine(ipair.second + jpair.second);
//                 }
//             }
//         }
//         tnodes[vid].edges = std::move(v);

//         if (tnodes[vid].edges.size()==0) 
//         {
//             core_vertexes.emplace_back(vid);
//         }
//     }
//     cal_parentchild_relationship_of_tree();
//     tend = clock();
//     // mlog("finish graph reduction in %f seconds.", etime - stime)
//     // mlog("there are %lu core vertexes.", core_vertexes.size())
// }



    // void dfs(NodeId i, vector<bool> &visited, int &size, const NodeId ccid) 
    // {
    //     visited[i] = true;
    //     sccid[i] = ccid;
    //     size++;
    //     for (const auto &a:graph[i]) 
    //     {
    //         if (!visited[a.dst]) {
    //             dfs(a.dst, visited, size, ccid);
    //         }
    //     }
    // }

    // void scc_cnt() 
    // {
    //     vector<bool> visited(graph.size(), false);
    //     sccid.resize(graph.size());

    //     int cnt = 0;
    //     NodeId ccid = 0;
    //     map<int, int> size_summary;
    //     for (NodeId i = 0; i < graph.size(); i++) 
    //     {
    //         if (!visited[i]) 
    //         {
    //             cnt++;
    //             int size = 0;
    //             dfs(i, visited, size, ccid++);
    //             if (size_summary.find(size)!=size_summary.end()) {
    //                 size_summary[size]++;
    //             } else {
    //                 size_summary[size] = 1;
    //             }
    //         }
    //     }
    // }



    // void graph_reduction()
    // {
    //     vector<vector<NodeId> > degree2nodeQ;//degree->a list of vertex of this degree
    //     vector<pair<NodeId, NodeId> > vPosition(graph.size());//(degree,idx)

    //     for (NodeId v = 0; v < graph.size(); ++v) 
    //     {
    //         NodeId degree = graph[v].size();
    //         if (degree >= degree2nodeQ.size()) {
    //             degree2nodeQ.resize(degree + 1);
    //         }
    //         vPosition[v] = make_pair(degree, degree2nodeQ[degree].size());
    //         degree2nodeQ[degree].push_back(v);
    //     }
        
    //     /////////////////////////////////////////////////////////////
    //     // for (int i = 0; i < degree2nodeQ.size(); i++)
    //     // {
    //     //     cout << "degree: " << i << endl;
    //     //     for (int j = 0; j < degree2nodeQ[i].size(); j++)
    //     //     {
    //     //         cout << degree2nodeQ[i][j] << ", " ;
    //     //     }
    //     //     cout << endl;
    //     //     cout << " ====================================== " << endl;
            
    //     // }
    //     /////////////////////////////////////////////////////////////

    //     vector<map<unsigned int, PLF> > shortcuts(graph.size());
    //     for (NodeId s = 0; s < graph.size(); s++)
    //     {
    //         for (unsigned int i = 0; i < graph[s].size(); ++i)
    //         {
    //             auto &e = graph[s][i];
    //             shortcuts[s][e.dst] = e.weights;
    //         }
    //     }  

    //     vertexOrder.resize(graph.size(), -1);
    //     NodeId mindegree = 0;
    //     tnodes.resize(graph.size());

    //     NodeId reduced_cnt = 0;
    //     unsigned int threshold_k = 0;
    //     vector<pair<unsigned int, unsigned int> > twidth_threshold2reduced_vcnt;
    //     mindegree = 0;  

    //     while (true) 
    //     {
    //         while (mindegree < degree2nodeQ.size() && degree2nodeQ[mindegree].empty()) 
    //         {
    //             // if (mindegree==threshold_k) 
    //             // {
    //             //     if (reduced_cnt!=0) 
    //             //     {
    //             //         cout << "reducting degree: " << mindegree << endl;
    //             //         //twidth_threshold2reduced_vcnt.emplace_back(mindegree, reduced_cnt);
    //             //         twidth_threshold2reduced_vcnt.push_back(make_pair(mindegree, reduced_cnt));
    //             //         reduced_cnt = 0;
    //             //     }
    //             //     threshold_k++;
    //             // }
    //             mindegree++;
    //         }

    //         // cout << "mindegree: " << mindegree <<endl;
    //         // break;
            
    //         //if (k <= mindegree || degree2nodeQ.size()==mindegree)
    //         if (degree2nodeQ.size()==mindegree)
    //             break;
    //         reduced_cnt++;
    //         NodeId vid = degree2nodeQ[mindegree].back();
    //         degree2nodeQ[mindegree].pop_back();
    //         vertexOrder[vid] = vorder++;

    //         //adding short cuts
    //         auto &v = shortcuts[vid];

    //         cout << "for nodeid: " << vid << endl;

    //         vector<unsigned int> valid_neighbor_index;
    //         for (map<unsigned int, PLF>::iterator it = v.begin(); it!=v.end(); ++it)
    //         {
    //             // cout << it->first << endl;
    //             // cout << it->second << endl;
    //             if (vertexOrder[it->first] == -1)
    //             {
    //                 valid_neighbor_index.push_back(it->first);
    //             }
                
    //         }
    //         //cout << " --------------------------------------- " << endl;

    //         vector<int> neighbor_degree_increase_cnt(valid_neighbor_index.size(), -1);
    //         //add short cuts
    //         cout << "add shortcuts" << endl;
    //         for (unsigned int ii = 0; ii < valid_neighbor_index.size(); ++ii) 
    //         {
    //             for (unsigned int jj = ii + 1; jj < valid_neighbor_index.size(); ++jj) 
    //             {
    //                 NodeId ivid  = valid_neighbor_index[ii];
    //                 NodeId jvid  = valid_neighbor_index[jj];

    //                 if (shortcuts[ivid].find(jvid) == shortcuts[jvid].end())
    //                 {
    //                     neighbor_degree_increase_cnt[ii]++;
    //                     neighbor_degree_increase_cnt[jj]++;
    //                 }


    //                 PLF PLFij, PLFji;
    //                 shortcuts[ivid][vid].compound(shortcuts[vid][jvid],PLFij,vid);
    //                 shortcuts[ivid][jvid] = PLFij;
    //                 shortcuts[jvid][vid].compound(shortcuts[vid][ivid],PLFji,vid);
    //                 shortcuts[jvid][ivid] = PLFji;
    //             }

    //         }

    //         cout << "update vPosition" << endl;
    //         //update vPosition
    //         for (unsigned int i = 0; i < valid_neighbor_index.size(); ++i) 
    //         {

    //             if (neighbor_degree_increase_cnt[i]!=0) 
    //             {
    //                 //unsigned int &x = v.idx(valid_neighbor_index[i]).first;
    //                 //cout << "unsigned int &x = v.idx(valid_neighbor_index[i]).first;" << endl;
    //                 //unsigned int &x = valid_neighbor_index[i];
    //                 int x = valid_neighbor_index[i];
    //                 pair<NodeId, NodeId> &p = vPosition[x];
    //                 //swap and delete
    //                 //cout << "swap and delete" << endl;
    //                 degree2nodeQ[p.first][p.second] = degree2nodeQ[p.first].back();
    //                 vPosition[degree2nodeQ[p.first].back()].second = p.second;
    //                 degree2nodeQ[p.first].pop_back();
    //                 //place in a new position
    //                 cout << "place in a new position" << endl;
    //                 //p.first += static_cast<unsigned int>(neighbor_degree_increase_cnt[i]);
    //                 p.first += neighbor_degree_increase_cnt[i];
    //                 if (p.first >= degree2nodeQ.size()) 
    //                 {
    //                     degree2nodeQ.resize(p.first + 1);
    //                 }
    //                 cout << "mindegree = min(mindegree, p.first);" << endl;
    //                 mindegree = min(mindegree, p.first);
    //                 cout << "p.second = degree2nodeQ[p.first].size();" << endl;
    //                 cout << "p.first: " << p.first << " " << "mindegree: " << mindegree << endl;
    //                 cout << "degree2nodeQ.size(): " << degree2nodeQ.size() << endl;
    //                 p.second = degree2nodeQ[p.first].size();
    //                 cout << " degree2nodeQ[p.first].push_back(x);" << endl;
    //                 degree2nodeQ[p.first].push_back(x);
    //             }
    //         }


    //         //cout << "building tnodes" << endl;
    //         //building tnodes
    //         for(auto i:valid_neighbor_index) 
    //         {
    //             tnodes[vid].edges[i] = v[i];
    //             //tnodes[vid].edges.insert(std::move(v.idx(i)));
    //         }
    //         if (tnodes[vid].edges.size()==0) {
    //             core_vertexes.emplace_back(vid);
    //         }
    //     }

    //     // //cout << "collect remaining shortcuts, as root" << endl;
    //     // for(auto &l1:degree2nodeQ) 
    //     // {
    //     //     //collect remaining shortcuts, as root
    //     //     for (NodeId &vid:l1) 
    //     //     {
    //     //         core_vertexes.push_back(vid);
    //     //     }
    //     // }
    //     // //cout << "cal_parentchild_relationship_of_tree();" << endl;
    //     // cal_parentchild_relationship_of_tree();
    //     // //cout << "calculate_height();" << endl;
    //     // calculate_height();        
    // }


    // vector<pair<unsigned int, PLF>> BU_qcomputingV2(NodeId bottom, NodeId up) 
    // /// calculating all possible shortest path among these vertexes.
    // {
    //     vector<pair<unsigned int, PLF> > result;
    //     if (bottom!=up) 
    //     {
    //         //flatten_hashmap<unsigned int, SCAttr> distmap = tnodes[bottom].edges;
    //         map<unsigned int, PLF> dismap = tnodes[bottom].edges;
    //         //auto anc = static_cast<NodeId>(tnodes[bottom].pnodeid);
    //         NodeId anc = tnodes[bottom].pnodeid;

    //         while (anc!=up) 
    //         {
    //             flatten_hashmap<unsigned int, SCAttr> tmp;
    //             for (auto &e:tnodes[anc].edges) {
    //                 distmap[e.first].combine(distmap[anc] + e.second);
    //                 swap(tmp[e.first], distmap[e.first]);
    //             }
    //             swap(tmp, distmap);
    //             anc = static_cast<NodeId>(tnodes[anc].pnodeid);
    //         }
    //         for (auto &p:distmap)
    //             result.emplace_back(p);
    //     } 
    //     else //if bottom is the up
    //     {   
    //         for (auto &e:tnodes[up].edges)
    //         {
    //             result.emplace_back(e);
    //         }
    //         result.emplace_back(up, SCAttr::self_loop());
    //     }
    //     return result;
    // }


    // PLF query(NodeId s, NodeId t) 
    // {
    //     PLF result;
    //     auto lca = qLCA(s, t);
    //     auto s_result = BU_qcomputingV2(s, lca.first);
    //     auto t_result = BU_qcomputingV2(t, lca.second);

    //     PLF min_edge_attr;
    //     for (auto &a:s_result) 
    //     {
    //         if (a.second.route_cnt()==0) continue;
    //         for (auto &b:t_result) 
    //         {
    //             if (b.second.route_cnt()==0) continue;
    //             if (a.first==b.first)
    //             result.combine(a.second + b.second);
    //             else { result.combine(a.second + b.second + cp_distv2(a.first, b.first)); }
    //         }
    //     }
    //     return result;
    // }


    //store index
            // for (unsigned int i = 0; i < tnodes.size(); i++)
        // {
        //     index_write << "vid: " << i << "order: " << vertexOrder[i] << "pnode: " << tnodes[i].pnodeid << "height: " << tnodes[i].height << "\n";
        //     index_write << "cnodes: ";
        //     for (auto c:tnodes[i].cnodeid)index_write << c << " ";
        //     index_write << "\n";

        //     index_write << "X(node) size: " << tnodes[i].edges.size() << "\n";
        //     for (map<unsigned int, PLF>::iterator it = tnodes[i].edges.begin(); it!=tnodes[i].edges.end(); ++it)
        //     {
        //         int segsize = 0;
        //         for (auto f : *it->second.f) segsize ++;
        //         index_write << it->first << " " << segsize << "\n";
        //         for (auto f : *it->second.f)
        //         {
        //             index_write << f.t << " " << f.w << " " << f.intv << "\n";
        //         }
        //     }
        // }


    // void calculate_height() 
    // {
    //     cout << "void calculate_height() " << endl;
    //     queue<NodeId> Q;
    //     for (auto v:core_vertexes) 
    //     {
    //         cout << "for (auto v:core_vertexes) " << endl;
    //         tnodes[v].height = 0;
    //         for (auto v1:tnodes[v].cnodeid) 
    //         {
    //             Q.push(v1);
    //             tnodes[v1].height = 1;
    //         }
    //     }

    //     while (!Q.empty()) 
    //     {
    //         cout << " while (!Q.empty())" << endl;
    //         auto v = Q.front();
    //         Q.pop();
    //         for (auto a:tnodes[v].cnodeid) 
    //         {
    //             Q.push(a);
    //             tnodes[a].height = tnodes[v].height + 1;
    //         }
    //     }
    // }

    // void load_graph() 
    // {
    //     clock_t tbegin, tend;
    //     tbegin = clock();

    //     ifstream readgraph;
    //     readgraph.open(tgrapgh_path.c_str());
    //     if (!readgraph.is_open()) {
    //         printf("%s does not exist\n", tgrapgh_path.c_str());
    //         exit(0);
    //     }
    //     readgraph >> n >> m;
    //     cout << "vertexes: " << n << " edges: " << m <<endl;
    //     graph.resize(n);

    //     string line;
    //     string delimiter = " ";
    //     getline(readgraph, line);
    //     while(getline(readgraph, line))
    //     {
    //         NodeId sid;
    //         vEdge ve;

    //         vector<string> edge = split(line, delimiter);
    //         sid = stoi(edge[0]);
    //         ve.dst = stoi(edge[1]);
    //         int numInterpolation = stoi(edge[2]);

    //         string weightline;
    //         getline(readgraph, weightline);
    //         vector<string> weightsstring = split(weightline, delimiter);
    //         int i = 0;

    //         while (i < weightsstring.size()-1)
    //         {
    //             double t = stod(weightsstring.at(i));
    //             double w = stod(weightsstring.at(i+1));
    //             Segment seg(t,w);
    //             ve.weights.f->push_back(seg);
    //             i = i + 2;
    //         }

    //         graph[sid].push_back(ve);
    //     }

    //     readgraph.close();

    //     tend = clock();
    //     cout << "Finish loading graph \t time cost: " + to_string((tend - tbegin) / CLOCKS_PER_SEC) + "s" << endl;
    //     // cout << "============================report graph========================" << endl;
    //     // for (int i = 0; i < graph.size(); i++)
    //     // {
    //     //     for (int j = 0; j < graph[i].size(); j++)
    //     //     {
    //     //         auto &e = graph[i][j];
    //     //         cout << i << " " << e.dst << endl;
    //     //         cout << e.weights << endl;
    //     //         cout << " -------------------------------- " << endl;
    //     //     }
            
    //     // }
    //     // cout << "==========================report end======================" << endl;
        
    // }


        // map<NodeId, PLF> Path_computing_bottom2up(NodeId bottom, NodeId up)
    // {
    //     map<unsigned int, PLF> results;
    //     NodeId s = bottom;
    //     if (bottom != up)
    //     {
    //         map<unsigned int, PLF> dismap = tnodes[bottom].edges;
    //         NodeId anc = tnodes[bottom].pnodeid;
    //         while (anc!=up)
    //         {
    //             vector<NodeId> uVertex = uVertexes(s,anc);
    //             vector<NodeId> vVetex = vVetexes(s, anc);

    //             cout << " ********************************************* " << endl;
    //             cout << "check node: " << anc << endl;
    //             // cout << "uVertexes: ";
    //             // for (auto &uid : uVertex)
    //             // {
    //             //     cout << uid << ",";
    //             // }
    //             // cout << endl;
    //             // cout << "vVetexes: ";
    //             // for (auto &vid : vVetex)
    //             // {
    //             //     cout << vid << ",";
    //             // }
    //             // cout << endl;            

    //             for (auto &uid : uVertex)
    //             {
    //                 PLF plfvu;
    //                 dismap[anc].compound(tnodes[anc].edges[uid],plfvu,anc);
    //                 dismap[uid] = plfvu;
    //                 cout << " uVertex: " << uid << endl;
    //                 cout << "source->" << anc << " compound to " << anc << "->" << uid << endl;
    //                 cout << dismap[anc] << " compound to " << tnodes[anc].edges[uid] << endl;
    //                 cout << "compouned: source->"<< uid <<" " << plfvu << endl;
    //             }
    //             cout << " U####################V " << endl;
    //             for (auto &vid : vVetex)
    //             {
    //                 PLF plfvu;
    //                 dismap[anc].compound(tnodes[anc].edges[vid],plfvu,anc);
                    
    //                 cout << " vVertex: " << vid << endl;
    //                 cout << "source->" << anc << " compound to " << anc << "->" << vid << endl;
    //                 cout << dismap[anc] << " compound to " << tnodes[anc].edges[vid] << endl;
    //                 cout << "compouned: source->"<< vid <<" " << plfvu << endl;
    //                 cout << "origin: source->"<< vid << dismap[vid] << endl;

    //                 dismap[vid].minimize(plfvu);
    //                 cout << dismap[vid] << endl;
    //             }

    //             s = anc;
    //             anc = tnodes[anc].pnodeid;
    //         }

    //         for (map<unsigned int, PLF>::iterator it = dismap.begin(); it!=dismap.end(); ++it)
    //         {
    //             results[it->first] = it->second;
    //         }
    //     }
    //     else
    //     {
    //         for (map<unsigned int, PLF>::iterator it = tnodes[up].edges.begin(); it!=tnodes[up].edges.end(); ++it)
    //         {
    //             results[it->first] = it->second;
    //         } 
    //     }

    //     return results;    
    // }


    //     map<NodeId, PLF> Path_computing_up2bottom(NodeId bottom, NodeId up)
    // {
    //     map<unsigned int, PLF> results;

    //     NodeId s = bottom;
    //     if (bottom != up)
    //     {
    //         map<unsigned int, PLF> dismap;
    //         //cout << "tnodes[bottom].rev_edges: " << tnodes[bottom].rev_edges.size() << endl;
    //         for (auto &p:tnodes[bottom].edges)
    //         {
    //             dismap[p.first] = tnodes[p.first].rev_edges[bottom];
    //             //dismap[p.first] = tnodes[bottom].rev_edges[p.first];
    //             //cout << "p.first: " << p.first<< endl;
    //             //cout << "tnodes[bottom].rev_edges: " << tnodes[p.first].rev_edges.size() << endl;
    //         }
    //         NodeId anc = tnodes[bottom].pnodeid;

    //         while (anc != up)
    //         {
    //             vector<NodeId> uVertex = uVertexes(s,anc);
    //             vector<NodeId> vVetex = vVetexes(s, anc); 
    //             for (auto &uid : uVertex)
    //             {
    //                 PLF plfvu;
    //                 //cout << "PLF plfvu;" << endl;

    //                 //tnodes[anc].rev_edges[uid].compound(dismap[anc],plfvu,anc);
    //                 tnodes[uid].rev_edges[anc].compound(dismap[anc],plfvu,anc);
    //                 //cout << "tnodes[anc].rev_edges[uid].compound" << endl;
    //                 dismap[uid] = plfvu;
    //             }

    //             s = anc;
    //             anc = tnodes[anc].pnodeid;               
    //         }
            
    //         for (map<unsigned int, PLF>::iterator it = dismap.begin(); it!=dismap.end(); ++it)
    //         {
    //             results[it->first] = it->second;
    //         }
    //     }
    //     else
    //     {
    //         //for (map<unsigned int, PLF>::iterator it = tnodes[up].rev_edges.begin(); it!=tnodes[up].rev_edges.end(); ++it)
    //         for (map<unsigned int, PLF>::iterator it = tnodes[up].edges.begin(); it!=tnodes[up].edges.end(); ++it)
    //         {
    //             results[it->first] = tnodes[bottom].rev_edges[it->first];
    //         } 
    //     }

    //     return results;   
  
         
    // }
