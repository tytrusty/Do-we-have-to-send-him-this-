#ifndef CLUSTER_H
#define CLUSTER_H

#include <igl/adjacency_list.h>
#include <vector>

using namespace std;
using namespace Eigen;

class Cluster
{


private:
public:
    Cluster(const MatrixXi& F, int clusters) : clusters_(clusters), 
        visited(F.rows(), false), cluster_map(clusters, vector<int>()),
        cluster_map_true(clusters, vector<int>())
    {
        igl::adjacency_list(F, A);
        visited.resize(F.rows());
    }
    vector<unsigned char> visited;
    vector<vector<int>> cluster_map;
    vector<vector<int>> cluster_map_true;
    vector<vector<int>> A;
    int clusters_;

    void BFS(int initial)
    {
        //vector<list<int>> queues;
        //visited[s] = true;
        //q.push_back(s);
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
