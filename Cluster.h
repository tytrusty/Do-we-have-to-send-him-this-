#ifndef CLUSTER_H
#define CLUSTER_H

#include <igl/adjacency_list.h>
#include <vector>
#include <random>

using namespace std;
using namespace Eigen;

class Cluster
{
private:
public:
    Cluster(const MatrixXi& F, int clusters) : clusters_(clusters), 
        visited(F.rows(), false), cluster_map(clusters, vector<int>()),
        cluster_map_true(clusters, vector<int>()), queues(clusters, list<int>())
    {
        igl::adjacency_list(F, A);
        std::default_random_engine gen;
        std::uniform_int_distribution<int> dis(0, F.rows()-1);
        for (int i = 0; i < clusters; i++) {
            queues[i].push_back(dis(gen));
        }
    }
    vector<unsigned char> visited;
    vector<list<int>> queues;
    vector<vector<int>> cluster_map;
    vector<vector<int>> cluster_map_true;
    vector<vector<int>> A;
    int clusters_;

    void BFS(int initial)
    {
        //while (!q.empty())
        //{
        //    s = q.front();
        //    cout<<s<<" ";
        //    q.pop_front();
        //    for(auto i : Edge[s])
        //    {
        //        if(!visited[i])
        //        {
        //            visited[i] = true;
        //            q.push_back(i);
        //        }
        //    }
		//}
    }
};
#endif // CLUSTER_H
