#ifndef BRANCHSIMPLIFICATION_H
#define BRANCHSIMPLIFICATION_H
#include "Graph/Graph.h"


class BranchSimplification
{
    public:
        BranchSimplification();
        static void ConvertToCustomGraph(std::vector<std::vector<int>> & paths, arma::mat & cloud, MyGraphType & G, double thresh, std::vector<int> & CloudToGraph);
        static void SimplifyIt(std::vector<std::vector<int>> & in, MyGraphType & Final, double threshold, arma::mat & cloud);
//void BranchSimplification::ConvertToCustomGraph(std::vector<std::vector<int>> & paths, arma::mat & cloud, MyGraphType & G, double thresh, std::vector<int> & CloudToGraph)


    protected:

    private:
};

#endif // BRANCHSIMPLIFICATION_H
