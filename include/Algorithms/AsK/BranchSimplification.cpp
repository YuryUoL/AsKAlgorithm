#include "BranchSimplification.h"


void BranchSimplification::ConvertToCustomGraph(std::vector<std::vector<int>> & paths, arma::mat & cloud, MyGraphType & G, double thresh, std::vector<int> & CloudToGraph)
{
//! First add all the indices
    std::vector<bool> duplicateRemover(cloud.n_cols,0);

    for (int i = 0; i < paths.size(); i++)
    {

        int first = paths[i][0];
        int second = paths[i][paths[i].size()-1];
        if (duplicateRemover[first] == false)
        {
            int indx = Graph::add_vertex(G, cloud.col(first));
            duplicateRemover[first] = true;
            CloudToGraph[first] = indx;
        }
        if (duplicateRemover[second] == false)
        {
            int indx = Graph::add_vertex(G, cloud.col(second));
            duplicateRemover[second] = true;
            CloudToGraph[second] = indx;
        }
    }
    for (int i = 0; i < paths.size(); i++)
    {
        double d = 0;
        std::vector<int> tmpContainer;
        int Ff = paths[i][0];
        int Ss = paths[i][paths[i].size()-1];
        for (int j = 0; j < paths[i].size() - 1; j++)
        {
            int first =paths[i][j];
            int second = paths[i][j+1];
            d = d + mlpack::metric::EuclideanDistance::Evaluate(cloud.col(first),cloud.col(second));
            tmpContainer.push_back(first);
            //Graph::add_edge(G, CloudToGraph[first], CloudToGraph[second]);
        }
        tmpContainer.push_back(paths[i][paths[i].size()-1]);
        if (d < thresh )
        {
            edge_descriptor e = Graph::add_custom_edge(G,Ff, Ss, d);
            G[e].indices = tmpContainer;
        }
    }
}

void BranchSimplification::SimplifyIt(std::vector<std::vector<int>> & in, MyGraphType & Final, double threshold, arma::mat & cloud)
{


    std::vector<int> CloudToGraph(cloud.n_cols);

    MyGraphType Simplified;

    BranchSimplification::ConvertToCustomGraph(in, cloud, Simplified,threshold, CloudToGraph);

    auto allVertices = boost::vertices(Simplified);


    std::vector<int> componentMap(num_vertices(Simplified));
    int num = boost::connected_components(Simplified, componentMap.data());


    std::vector<bool> duplicateCatcher(boost::num_vertices(Simplified));
    std::vector<arma::vec> barycenterSums(num, arma::vec(3));
    for (int i = 0; i < num; i++)
    {
        barycenterSums[i].fill(0);
    }
    std::vector<double> sumComputation(num,0);

    auto edges = boost::edges(Simplified);

    for (int i = 0; i < in.size(); i++)
    {
        int first = in[i][0];
        int second = in[i][in[i].size()-1];
       if (duplicateCatcher[first] == false)
        {
            int vertex = CloudToGraph[first];
            barycenterSums[componentMap[vertex]] = barycenterSums[componentMap[vertex]] + cloud.col(first);
            sumComputation[componentMap[vertex]] = sumComputation[componentMap[vertex]] + (double) 1;
            duplicateCatcher[first] = true;
        }
        if (duplicateCatcher[second] == false)
        {
            int vertex = CloudToGraph[second];
            barycenterSums[componentMap[vertex]] = barycenterSums[componentMap[vertex]] + cloud.col(second);
            sumComputation[componentMap[vertex]] = sumComputation[componentMap[vertex]] + (double) 1;
            duplicateCatcher[second] = true;
        }

    }

    for (auto edgeit = edges.first ; edgeit != edges.second ; edgeit++)
    {
        int vertex = boost::source(*edgeit,Simplified);
        for(int i = 0; i < Simplified[*edgeit].indices.size(); i++)
        {
            int indx = Simplified[*edgeit].indices[i];
            if (duplicateCatcher[indx] == false)
            {
                barycenterSums[componentMap[vertex]] = barycenterSums[componentMap[vertex]] + cloud.col(indx);
                sumComputation[componentMap[vertex]] = sumComputation[componentMap[vertex]] + (double) 1;
                duplicateCatcher[indx] = true;
            }
        }
    }
    for (int i = 0; i < num; i++)
    {
        barycenterSums[i] = barycenterSums[i] / sumComputation[i];
    }


    std::vector<bool> newDuplicateCatcher(cloud.n_cols);


    for (int i = 0; i < num; i++)
    {
        int jjj = Graph::add_vertex(Final,barycenterSums[i]);
    }


    for (int i = 0; i < in.size(); i++)
    {
        int lastindx = in[i].size()-1;
        if (componentMap[CloudToGraph[in[i][0]]] != componentMap[CloudToGraph[in[i][lastindx]]])
        {
            if(in[i].size() == 2)
            {
                Graph::add_edge(Final, componentMap[CloudToGraph[in[i][0]]],componentMap[CloudToGraph[in[i][lastindx]]]);
            }
            int prev = -1;
            int prev2 = -1;
            for (int j = 0;  j < in[i].size() - 1; j++)
            {
                if (j == 0)
                {
                    prev = Graph::add_vertex(Final, cloud.col(in[i][1]));
                    Graph::add_edge(Final, prev, componentMap[CloudToGraph[in[i][0]]]);
                }
                else if (j == in[i].size() - 2)
                {
                    Graph::add_edge(Final, prev, componentMap[CloudToGraph[in[i][lastindx]]]);



                }
                else
                {

                    prev2 = Graph::add_vertex(Final, cloud.col(in[i][j+1]));
                    Graph::add_edge(Final, prev, prev2);
                    prev = prev2;

                }

            }
        }

    }


//auto edgePair = boost::edges(Simplified)



}



