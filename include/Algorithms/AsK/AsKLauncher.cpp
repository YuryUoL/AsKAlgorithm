#include "Algorithms/AsK/AsKLauncher.h"
#include <mlpack/methods/emst/dtb.hpp>
#include "Algorithms/Convertors/MatToGraph.h"
#include "Algorithms/AsK/BranchDetection.h"
#include "Algorithms/AsK/CorrectStraightening.h"
#include "Algorithms/AsK/BranchSimplification.h"


AsKLauncher::AsKLauncher(double a, double b, double c, std::string name = "AsK"):
    branch_detection(a),approxError(b),simplificationError(c),numberOfRuns(0)
{
}
void AsKLauncher::Run(arma::mat & data, MyGraphType & out)
{
    //! DebugPrinting

    mlpack::emst::DualTreeBoruvka<> MSTOP(data);
    arma::mat results;
    MSTOP.ComputeMST(results);

    MyGraphType G;

    MatToGraph::Evaluate(data,results,G);

    MyGraphType Simplified;

    double branchingParameter = this->branch_detection * this->ComputeAverageEdge((G));

    std::vector<int> oldToNew;
    std::vector<int> newToOld;
    std::vector<int> correspondance;

    BranchDetection::branchIndexation(G, branchingParameter, Simplified, oldToNew,  newToOld, correspondance);
    // out = Simplified;

    //BranchDetection::PrintStupidStuff();
    //! static void branchIndexation(MyGraphType & G, double thresh, MyGraphType & out, std::vector<double> & oldToNew, std::vector<double> & NewToOld);
    std::vector<std::vector<int>> SimplifiedPaths;
    //!  static double ComputeStraightening(MyGraphType & G, arma::mat & cloud, std::vector<std::vector<int>> & out, double e, std::vector<int> & CloudToGraph, std::vector<int> & GraphToCloud)
    //OutputPrinter::DebugGraph("/home/yury/LocalTests4/Debug/zero.txt" , Simplified);

    double distanceThreshttt = CorrectStraightening::ComputeStraightening(Simplified, data, SimplifiedPaths, this->approxError, oldToNew, newToOld, correspondance);


    BranchSimplification::SimplifyIt(SimplifiedPaths, out, this->simplificationError * distanceThreshttt, data);
//
//    //! Temporalily:
//    std::cout << "Final paths of the form: " << SimplifiedPaths.size() << std::endl;
//    for (int i = 0; i < SimplifiedPaths.size(); i++)
//    {
//        std::cout << "Another size: " << i << ": " << SimplifiedPaths[i].size() << std::endl;
//
//    }
//! Temporal
   //this->ConvertToGraph(SimplifiedPaths, data, out);
//    //BranchSimplification::
//
//
//    std::string actualSettings = this->settings;

}

void AsKLauncher::ConvertToGraph(std::vector<std::vector<int>> & paths, arma::mat & cloud, MyGraphType & G)
{
//! First add all the indices
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



double AsKLauncher::ComputeAverageEdge(MyGraphType & G)
{
    auto vpair = boost::edges(G);
    double sum = 0;
    for (auto it = vpair.first ; it != vpair.second ; it++)
    {

        sum = sum + G[*it].distance;

    }
    return sum / (double) boost::num_edges(G);
}

//void AsKLauncher::Convert()
//{




//}

//void ClassicUIHandle::CorrectCalculationToGraphIntermediate(MyGraphType & FG, std::list<std::list<Point>> & optipath,
//std::list<std::list<Point>> & branchsimplified, std::list<Point> & cloud, double & mstlength, std::string settings)
//{
//MyGraphType tree = Computation::computeMST(cloud);
//MyGraphType optiout;
//mstlength = Computation::AverageEdgelength(tree);
//double param1 = ClassicUIHandle::branchDetect*mstlength;
//BranchDetection::SimplifyIt(tree,optiout,param1,"",settings);
//double ddd = StraighteningMethods::ClassicStraightening(optiout, cloud, optipath, this->Straightening );
//std::cout << "The straightening parameter: " << ddd << std::endl;
//double valuev = ddd*ClassicUIHandle::branchCollapse;
//BranchSimplification::SimplifyIt(optipath, branchsimplified,valuev);
//BranchSimplification::PathToGraphProper(FG, branchsimplified);
//}


//class AsKLauncher : public AbstractAlgorithm
//{
//    public:
//        AsKLauncher(double a, double b, double c);
//        void Run(std::list<Point> & cloudlist, MyGraphType & out) override;
//
//      //  virtual ~AsKLauncher();
//
//    protected:
//
//    private:
//    double branch_detection;
//    double approxError;
//    double simplificationError;
//};
