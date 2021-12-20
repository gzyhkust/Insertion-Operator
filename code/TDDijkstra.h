#ifndef TDGT_TDDIJKSTRA_H
#define TDGT_TDDIJKSTRA_H

#include <vector>
#include <queue>
#include "Constants.h"
#include "TGraph.h"

pair<int, double> tddijkstra(int s, int d, double td, TGraph &TG, std::vector<int> &prevexs) {
    auto cmp = [](std::pair<int, double> left, std::pair<int, double> right) {
        return left.second > right.second;
    };
    prevexs = std::vector<int>(TG.n, -1);
    std::vector<double> arrs(TG.n, INT_MAX);
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, decltype(cmp)>
            vexs_pq(cmp);
    vexs_pq.emplace(s, td);
    arrs[s] = td;
    int cnt_visited = 0;
    while (!vexs_pq.empty()) {
        auto vex = vexs_pq.top();
        vexs_pq.pop();
        if (vex.first == d)//done
            return make_pair(cnt_visited, arrs[d]-td);
        if (neq(vex.second, arrs[vex.first])) //only settle for entry with shortest travel time
            continue;
        cnt_visited++;
        for (int j = TG.head[vex.first]; j != -1; j = TG.next[j]) {
            double arr_new = TG.weights[j].dpt2arr(vex.second);
            assert(lt(arr_new, INT_MAX));
            if (lt(arr_new, arrs[TG.adjv[j]])) { // relax adj vertices
                arrs[TG.adjv[j]] = arr_new;
                prevexs[TG.adjv[j]] = vex.first;
                vexs_pq.emplace(TG.adjv[j], arr_new);
            }
        }
    }
    return make_pair(-1,-1);
}
void recover_path_from_prevs(std::ostream &out, int s, int d, std::vector<int> &prevs) {
    assert(prevs.size());
    std::stack<int> path;
    int vex = d;
    while (vex != s) {
        path.push(vex);
        vex = prevs[vex];
    }
    path.push(vex);
    while (!path.empty()) {
        out << path.top() << " ";
        path.pop();
    }
    out << std::endl;
}

#endif //TDGT_TDDIJKSTRA_H
