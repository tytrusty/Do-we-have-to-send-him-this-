#ifndef CLUSTER_H
#define CLUSTER_H

#include <igl/adjacency_list.h>
#include <vector>
#include <random>
#include <unordered_set>

using namespace std;
using namespace Eigen;

class Cluster
{
private:
public:
    Cluster(const MatrixXi& F, int clusters, int nvertices) : clusters_(clusters), 
        visited(nvertices, false), cluster_map(clusters, unordered_set<int>()),
        cluster_map_true(clusters, unordered_set<int>()), queues(clusters, list<int>())
    {
        igl::adjacency_list(F, vertex_to_vertices);
        std::default_random_engine gen;
        std::uniform_int_distribution<int> dis(0, nvertices-1);
        for (int i = 0; i < clusters; i++) {
            queues[i].push_back(dis(gen));
        }
    }
    vector<unsigned char> visited;
    vector<list<int>> queues;
    vector<unordered_set<int>> cluster_map;
    vector<unordered_set<int>> cluster_map_true;
    vector<vector<int>> vertex_to_vertices;
    int clusters_;

    bool queues_not_empty()
    {
        for (list<int> q : queues) {
            if (!q.empty()) {
                return true;
            }
        }
        return false;

    }

    void BFS()
    {
        std::cout << "BEGIN BFS" << std::endl;
        while (queues_not_empty()) {
            for (int i = 0; i < queues.size(); i++) {
                if (queues[i].empty()) continue;
                int v = queues[i].front(); 
                queues[i].pop_front();
                for (int adj : vertex_to_vertices[v]) {
                    if (!visited[adj]) {
                        visited[adj] = true;
                        queues[i].push_back(adj);
                        cluster_map_true[i].insert(adj);
                    } else { // If visited, still need adjacent 
                        cluster_map[i].insert(adj);
                    }

                }
            }
        }
        std::cout << "END BFS" << std::endl;
    }
};
#endif // CLUSTER_H
