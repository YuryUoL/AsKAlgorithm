#ifndef GRAPH_H
#define GRAPH_H

#include <mlpack/core.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/program_options.hpp>
#include <boost/graph/graphviz.hpp>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/connected_components.hpp>

struct VertexData;
struct EdgeData;

//! This is graph with no extra information (can be used to compute connected components for example)
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> AbstractGraphType;

//! This is the customized graph type for R^3
typedef boost::adjacency_list<boost::hash_setS, boost::vecS,
        boost::undirectedS,
        VertexData,
        boost::property<boost::edge_weight_t, double, EdgeData>
        > MyGraphType;

typedef typename boost::graph_traits<MyGraphType>::vertex_descriptor vertex_descriptor;
typedef typename boost::graph_traits<MyGraphType>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<MyGraphType>::vertex_iterator vertex_iter;
typedef boost::graph_traits<MyGraphType>::edge_iterator edge_iter;


struct EdgeData
{
    double distance;
    std::vector<int> indices;
    arma::mat segment;
};
struct VertexData
{
    arma::vec p;
    //vertex_descriptor index;
    //vertex_descriptor correspondance;
    //double potentialvalue;
    //double maximumlength;
    //bool branchedcandidate;
    //bool boundary;
    //int interval;
    //vertex_descriptor prev;
    //vertex_descriptor prev2;
    //std::vector<std::pair<vertex_descriptor,MyGraphType*>> connectors;

};

class Graph
{
public:
    Graph();
    static vertex_descriptor add_vertex(MyGraphType &G, arma::vec p);
    static edge_descriptor add_edge(MyGraphType &G,
                         vertex_descriptor v1,
                        vertex_descriptor v2);

    static edge_descriptor add_custom_edge(MyGraphType &G,
                                vertex_descriptor v1,
                                vertex_descriptor v2, double d);


protected:

private:
};






#endif // GRAPH_H
