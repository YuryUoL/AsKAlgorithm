#ifndef POSTPROCESSOR_H
#define POSTPROCESSOR_H


#include "Graph/Graph.h"
#include <mlpack/core.hpp>
#include "Algorithms/Convertors/npy.hpp"

class PostProcessor
{
public:
    PostProcessor() {}

    static void CloudToNpy(std::string cloudout, arma::mat & cloud)
    {
        std::vector<std::vector<double>> cloudV(cloud.n_cols, std::vector<double>(cloud.n_rows));
        for (int i = 0 ; i < cloud.n_cols; i++)
        {
            for (int j = 0; j < cloud.n_rows; j++)
            {
                cloudV[i][j] = cloud(j,i);
            }
        }
        std::vector<double> cloudVinv;
        PostProcessor::Convertor(cloudV, cloudVinv);
        const long unsigned leshape [] = {cloud.n_cols,cloud.n_rows};
        npy::SaveArrayAsNumpy(cloudout, false, 2, leshape, cloudVinv);


    }

    static void GraphToNpy(std::string abstractGraph, std::string ConcreteGraph, std::string allocationOut, std::string pointIndicesOut, MyGraphType & G, std::vector<int> & allocation, std::vector<int> & origAlloc)
    {

        std::vector<std::vector<int>> edgeStructure(boost::num_edges(G), std::vector<int>(2));
        std::vector<edge_descriptor> edges(boost::num_edges(G));
        int dim = G[0].p.n_elem;
        std::vector<std::vector<double>> cloudV(boost::num_vertices(G), std::vector<double>(dim));

        auto edgemass = boost::edges(G);
        int countt = 0;
        int maxDegree = 0;
        for (auto eit = edgemass.first;  eit != edgemass.second; eit++)
        {
            edges[countt] = *eit;
            edgeStructure[countt][0] = boost::source(*eit,G);
            edgeStructure[countt][1] = boost::target(*eit,G);
            G[*eit].order = countt;
            countt++;
        }
        for (int i = 0; i < boost::num_vertices(G); i++)
        {
            int opq = boost::out_degree(i,G);
            maxDegree = std::max(maxDegree, opq);
            for (int j = 0; j < dim; j++)
            {
                cloudV[i][j] = G[i].p(j);
            }
        }


        std::vector<std::vector<int>> allocationArray(allocation.size(), std::vector<int>(maxDegree));
        std::vector<std::vector<int>> pointAllocationArray(boost::num_vertices(G), std::vector<int>());

        for (int i = 0; i < allocation.size(); i++)
        {
            pointAllocationArray[origAlloc[i]].push_back(i);
            int kk = boost::degree(allocation[i],G);
            int j = 0;
            auto bbb = boost::out_edges(allocation[i],G);
            for (auto it = bbb.first; it != bbb.second; it++)
            {
                allocationArray[i][j] = G[*it].order;
                j++;
            }
            for (int abc = j; abc < maxDegree; abc++)
            {
                allocationArray[i][abc] = allocationArray[i][0];
            }

        }



        int maxSize = 0;
        for (int i = 0; i < boost::num_vertices(G); i++)
        {
            maxSize = std::max((int) pointAllocationArray[i].size(), maxSize);
        }


        for (int i = 0; i < boost::num_vertices(G); i++)
        {
            for (int j =  pointAllocationArray[i].size(); j < maxSize; j++ )
            {
                int kk = pointAllocationArray[i][0];
                (pointAllocationArray[i]).push_back(kk);
            }

        }



        std::vector<int> edgeStructureReal;
        PostProcessor::Convertor(edgeStructure, edgeStructureReal);

        const long unsigned leshape [] = {boost::num_edges(G),2};
        npy::SaveArrayAsNumpy(abstractGraph, false, 2, leshape, edgeStructureReal);

        std::vector<double> cloudVReal;
        PostProcessor::Convertor(cloudV, cloudVReal);

        const long unsigned leshapeTwo [] = {boost::num_vertices(G),dim};
        npy::SaveArrayAsNumpy(ConcreteGraph, false, 2, leshapeTwo, cloudVReal);

        std::vector<int> allReal;
        PostProcessor::Convertor(allocationArray, allReal);
        const long unsigned leshapeThree [] = {allocation.size(),maxDegree};
        npy::SaveArrayAsNumpy(allocationOut, false, 2, leshapeThree, allReal);


        std::vector<int> tmpAsdAsd;
        PostProcessor::Convertor(pointAllocationArray, tmpAsdAsd);
        const long unsigned leshapeFour [] = {boost::num_vertices(G),maxSize};
        npy::SaveArrayAsNumpy(pointIndicesOut, false, 2, leshapeFour, tmpAsdAsd);



    }

    static void Convertor(std::vector<std::vector<double>> & in , std::vector<double> & out )
    {
        //out.reserve(in.size() * in[0].size());
        for (int i = 0; i < in.size(); i++)
        {
            for (int j  = 0 ; j < in[i].size(); j++)
            {
                out.push_back(in[i][j]);
            }

        }
    }
    static void Convertor(std::vector<std::vector<int>> & in , std::vector<int> & out )
    {
        //out.reserve(in.size() * in[0].size());
        for (int i = 0; i < in.size(); i++)
        {
            for (int j  = 0 ; j < in[i].size(); j++)
            {
                out.push_back(in[i][j]);
            }

        }
    }

    static void GraphToNpySimple(std::string abstractGraph, std::string ConcreteGraph, std::string allocationOut, MyGraphType & G, int numCloud)
    {
        std::vector<std::vector<int>> edgeStructure(boost::num_edges(G), std::vector<int>(2));
        int dim = G[0].p.n_elem;

        auto edgemass = boost::edges(G);

        int countt = 0;
        for (auto eit = edgemass.first;  eit != edgemass.second; eit++)
        {
            edgeStructure[countt][0] = boost::source(*eit,G);
            edgeStructure[countt][1] = boost::target(*eit,G);
            countt++;
        }

        const long unsigned leshape [] = {boost::num_edges(G),2};

        //! Convertor(std::vector<std::vector<double>> & in , std::vector<double> & out )
        std::vector<int> edgeStructureConv;
        PostProcessor::Convertor(edgeStructure , edgeStructureConv );

        npy::SaveArrayAsNumpy(abstractGraph, false, 2, leshape, edgeStructureConv);
        std::vector<std::vector<double>> cloudV(boost::num_vertices(G), std::vector<double>(dim));




        for (int i = 0; i < boost::num_vertices(G); i++)
        {
            for (int j = 0; j < dim; j++)
            {
                cloudV[i][j] = G[i].p(j);
            }
        }

        std::vector<double> cloudVConv;
        PostProcessor::Convertor(cloudV , cloudVConv );

        const long unsigned leshapeTwo [] = {boost::num_vertices(G),dim};
        npy::SaveArrayAsNumpy(ConcreteGraph, false, 2, leshapeTwo, cloudVConv);



        std::vector<std::vector<int>> allocationArray(numCloud, std::vector<int>(boost::num_edges(G)));


        for (int i = 0; i < numCloud; i++)
        {

            for (int j = 0; j < boost::num_edges(G); j++)
            {
                allocationArray[i][j] = j;
            }

        }
        std::vector<int> allocationArrayConv;
        PostProcessor::Convertor(allocationArray , allocationArrayConv );
        const long unsigned leshapeThree [] = {numCloud, boost::num_edges(G)};
        npy::SaveArrayAsNumpy(allocationOut, false, 2, leshapeThree, allocationArrayConv);



    }

    static void PostProcessSimple(std::string cloudout, std::string abstractGraph, std::string ConcreteGraph, std::string allocationOut, MyGraphType & G, arma::mat & cloud)
    {
            PostProcessor::CloudToNpy(cloudout, cloud);
            PostProcessor::GraphToNpySimple(abstractGraph, ConcreteGraph, allocationOut, G, cloud.n_cols);


    }

    //        PostProcessor::PostProcess(this->cloudout, this->abstractGraphOut, this->concreateGraphOut, this->belonghoodOut,  belonghood, G,  Cloud);

    static void PostProcess(std::string cloudout, std::string abstractGraph, std::string ConcreteGraph, std::string allocationOut, std::string pointIndicesOut, MyGraphType & G, std::vector<int> & allocation, std::vector<int> & origAlloc, arma::mat & cloud)
    {
        PostProcessor::CloudToNpy(cloudout, cloud);
        PostProcessor::GraphToNpy(abstractGraph, ConcreteGraph, allocationOut, pointIndicesOut, G, allocation, origAlloc);
   //     PostProcessor::GraphToNpy(abstractGraph, ConcreteGraph, allocationOut, G, allocation);

    }

protected:

private:
};

#endif // POSTPROCESSOR_H
