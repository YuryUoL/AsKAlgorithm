#ifndef PATHTOGRAPH_H
#define PATHTOGRAPH_H

#include "Graph/Graph.h"

class PathToGraph
{
public:
    PathToGraph() {}
    static void Convert(std::vector<std::vector<int>> & paths, arma::mat & cloud, MyGraphType & G)
    {
        std::vector<int> CloudToGraph(cloud.n_cols);
        std::vector<bool> duplicateRemover(cloud.n_cols);
        for (int i = 0; i < paths.size(); i++)
        {
            for (int j = 0; j < paths[i].size(); j++)
            {
                if (duplicateRemover[i] == false)
                {
                    int indx = Graph::add_vertex(G, cloud.col(paths[i][j]));
                    duplicateRemover[i] == true;
                    CloudToGraph[paths[i][j]] = indx;
                }

            }
        }
        for (int i = 0; i < paths.size(); i++)
        {
            for (int j = 0; j < paths[i].size() - 1; j++)
            {
                int first =paths[i][j];
                int second = paths[i][j+1];
                Graph::add_edge(G, CloudToGraph[first], CloudToGraph[second]);
            }
        }
    }


protected:

private:
};

#endif // PATHTOGRAPH_H
