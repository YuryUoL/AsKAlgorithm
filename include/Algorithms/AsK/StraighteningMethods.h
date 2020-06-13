//#ifndef STRAIGHTENINGMETHODS_H
//#define STRAIGHTENINGMETHODS_H
#include "Graph/Graph.h"
//
//typedef std::pair<Point,std::list<Point>> PointPlus;
//
//class StraighteningMethods
//{
//public:
//    StraighteningMethods();
//    static bool OptimizeSingle(std::vector<PointPlus> & path,  std::list<Point> & thepath, double e);
//    static void Optimize(std::list<std::list<Point>> & out, std::vector<std::vector<PointPlus>> & paths, double e);
//    static void Allocate(std::vector<std::vector<PointPlus>> & paths, MyGraphType & G, std::list<Point> & cloud);
//    static void GraphToPaths(std::vector<std::vector<PointPlus>> & paths, MyGraphType & G);
//    static double ClassicStraightening(MyGraphType & G, std::list<Point> & cloud, std::list<std::list<Point>> & out, double e);
//
//    //! New straightening methods:
//    static void ConvertGraphToArmaMat(arma::mat & pathpoints, MyGraphType & G);
//    static void DualTreeAllocation(MyGraphType & G, std::list<Point> & cloud, arma::Mat<size_t> & theNeighbors);
//    static void DualTreeAllocatorPostProcessor(MyGraphType & G, arma::Mat<size_t> & theNeighbors, std::list<Point> & clouds, std::vector<std::list<Point>> & pAllocation);
//    static void addToPathNew(std::vector<PointPlus> & path, MyGraphType & G, vertex_descriptor v,std::vector<std::list<Point>> & pAllocation);
//    static void GraphToPathsNew(std::vector<std::vector<PointPlus>> & paths, MyGraphType & G, std::vector<std::list<Point>> & pAllocation);
//    static double UpgradedStraightening(MyGraphType & G, std::list<Point> & cloud, std::list<std::list<Point>> & out, double e,arma::Mat<size_t> & theNeighbors);
//
//
//
//protected:
//
//private:
//};
//
//#endif // STRAIGHTENINGMETHODS_H
