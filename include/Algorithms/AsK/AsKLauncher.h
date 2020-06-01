#ifndef ASKLAUNCHER_H
#define ASKLAUNCHER_H

#include <mlpack/core.hpp>
#include "Graph/Graph.h"

class AsKLauncher
{
    public:
        AsKLauncher(double a, double b, double c,std::string name);

        void Run(arma::mat & data, MyGraphType & out);
        double ComputeAverageEdge(MyGraphType & G);
        void ConvertToGraph(std::vector<std::vector<int>> & paths, arma::mat & cloud, MyGraphType & G);

      //  virtual ~AsKLauncher();

    protected:

    private:
    double branch_detection;
    double branchd;
    double approxError;
    double simplificationError;
    int numberOfRuns;
    std::string settings;
};

#endif // ASKLAUNCHER_H
