#ifndef ASKLAUNCHER_H
#define ASKLAUNCHER_H

#include <mlpack/core.hpp>
#include "Graph/Graph.h"

class AsKLauncher
{
    public:
        AsKLauncher(double a, double b, double c,std::string name, bool simpleMode, bool postProcessSwitch);

        void Run(arma::mat & data, MyGraphType & out);
        double ComputeAverageEdge(MyGraphType & G);
        void ConvertToGraph(std::vector<std::vector<int>> & paths, arma::mat & cloud, MyGraphType & G);
        void CompletePostProcessing(std::vector<int> & belonghood, std::vector<int> & originAlloc, MyGraphType & G, arma::mat & Cloud );
        void SetPostProcessing(std::string cloudout, std::string abstractGraphOut, std::string GraphOut, std::string belonghoodOut, std::string PointIndices);
        void ReIndexBelonghood(std::vector<int> & belonghood, std::vector<int> & oldToNew);
        void FixedBelonghood(std::vector<int> & belonghood, MyGraphType & out);



      //  virtual ~AsKLauncher();

    protected:

    private:
    bool postProcessSwitch;
    bool simpleMode;
    double branch_detection;
    double branchd;
    double approxError;
    double simplificationError;
    int numberOfRuns;
    std::string settings;

    //! Post processing parameters:
    //    static void PostProcess(std::string cloudout, std::string abstractGraph, std::string ConcreteGraph, std::string allocationOut , MyGraphType & G, std::vector<int> & allocation, arma::mat & cloud)
    std::string cloudout;
    std::string abstractGraphOut;
    std::string concreateGraphOut;
    std::string belonghoodOut;
    std::string PointIndices;





};

#endif // ASKLAUNCHER_H
