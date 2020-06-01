#include "Graph.h"


Graph::Graph()
{
}
//  Point p;
//  double potentialvalue;
//  double maximumlength;
//  bool branchedcandidate;


vertex_descriptor Graph::add_vertex(MyGraphType &G, arma::vec qq)
{
    vertex_descriptor v = boost::add_vertex(G);
    G[v].p = qq;
    return v;
}



edge_descriptor Graph::add_edge(MyGraphType &G,
                                vertex_descriptor v1,
                                vertex_descriptor v2)
{
    double d = mlpack::metric::EuclideanDistance::Evaluate(G[v1].p, G[v2].p);
    edge_descriptor e = boost::add_edge(v1, v2, G).first;
    G[e].segment.resize(G[v1].p.n_rows,2);
    G[e].segment.col(0) = G[v1].p;
    G[e].segment.col(1) = G[v2].p;
    boost::property_map<MyGraphType, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, G);
    weightmap[e] = d;
    G[e].distance = d;
    //!
    G[e].indices = std::vector<int>();
    return e;
}

edge_descriptor Graph::add_custom_edge(MyGraphType &G,
                                vertex_descriptor v1,
                                vertex_descriptor v2, double d)
{
    edge_descriptor e = boost::add_edge(v1, v2, G).first;
    G[e].segment.resize(G[v1].p.n_rows,2);
    G[e].segment.col(0) = G[v1].p;
    G[e].segment.col(1) = G[v2].p;
    boost::property_map<MyGraphType, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, G);
    weightmap[e] = d;
    G[e].distance = d;
    //!
    G[e].indices = std::vector<int>();
    return e;
}



