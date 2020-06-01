#ifndef MATTOGRAPH_H
#define MATTOGRAPH_H

#include "Graph/Graph.h"
#include <mlpack/core.hpp>

class MatToGraph
{
public:
    MatToGraph() {}
    static void Evaluate(arma::mat & originalData, arma::mat & edges, MyGraphType & G)
    {
        int vSize = originalData.n_cols;
        for (int i = 0; i < vSize; i++)
        {
            Graph::add_vertex(G,originalData.col(i));
        }

        int edgeSize = edges.n_cols;
      //  std::cout << "Number of edges in original mst: " << edgeSize << std::endl;
        for (int i = 0; i < edgeSize; i++)
        {
            int one = (int) edges(0,i);
            int two = (int) edges(1,i);
            Graph::add_edge(G,one,two);
        }
      //  std::cout << "Graph converted: " << std::endl;
     //   std::cout << "Number of vertices: " << boost::num_vertices(G) << std::endl;
      //  std::cout << "Number of edges: " << boost::num_edges(G) << std::endl;
    }

    static void Reverse(arma::mat & data, MyGraphType & G)
    {
        data.resize(G[0].p.n_elem, boost::num_vertices(G));
        for (int i = 0; i < boost::num_vertices(G); i++)
        {
            data.col(i) = G[i].p;
        }


    }
    // static void MatToMyGraphType(arma::mat & originalData, arma::mat & edges, MyGraphType & G);
    virtual ~MatToGraph();

protected:

private:
};

#endif // MATTOGRAPH_H
