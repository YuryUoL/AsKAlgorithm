#ifndef BRANCHDETECTION_H
#define BRANCHDETECTION_H

#include "Graph/Graph.h"
#include <functional>
#include <queue>
#include <vector>

typedef std::priority_queue<std::pair<int,double>, std::vector<std::pair<int,double>>,
        std::function<bool(std::pair<int,double> & a, std::pair<int,double> & b)>> VH;

class BranchDetection
{
public:
    BranchDetection();
    static void branchIndexation(MyGraphType & G, double thresh, MyGraphType & out, std::vector<int> & oldToNew, std::vector<int> & NewToOld, std::vector<int> & next);

    //! void BranchDetection::branchIndexation(MyGraphType & G, double thresh, MyGraphType & out, std::vector<int> & OldToNew, std::vector<int> & NewToOld, std::vector<int> & next)

    static void CorersponanceMarker( std::vector<bool> & coloring, std::vector<int> & next, std::vector<int> & edgePoints);


protected:

private:
};

#endif // BRANCHDETECTION_H
