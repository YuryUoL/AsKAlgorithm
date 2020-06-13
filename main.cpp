#include <iostream>
#include "Algorithms/Convertors/OutputPrinter.h"
#include "Algorithms/Convertors/InputReader.h"
//#include "Algorithms/Amst/AMSTComputator.h"
//#include "Algorithms/Kde/Functions/ExponentialDensity.h"
#include "Algorithms/AsK/AsKLauncher.h"
#include <boost/filesystem.hpp>

using namespace std;
//
//
//
//void TestNewAmst(std::string inputpath, std::string outputpath )
//{
//    arma::mat cloud;
//    MyGraphType G;
//    InputReader::XYZtoMAT(inputpath,cloud);
//    AMSTComputator<mlpack::metric::EuclideanDistance, arma::mat, ExponentialDensity,mlpack::tree::KDTree> comp(0.01,2.5);
//    cloud = cloud.t();
//    comp.Run(cloud, G);
//    OutputPrinter::GraphToVtk(outputpath,G);
//
//}

void TestAsK(std::string folder, std::string inputpath, double param = 20.0)
{
    std::string outputpath = folder + "InitialOutput.vtk";
    arma::mat cloud;
    MyGraphType G;
    InputReader::XYZtoMAT(inputpath,cloud);
    AsKLauncher launcher(param, 1.0, 1.1, "AsK",false, false);
//!    launcher.SetPostProcessing(folder + "cloud.npy", folder + "abstractGraph.npy", folder + "Edges.npy", folder + "Belonghood.npy");

    cloud = cloud.t();

    launcher.Run(cloud, G);
    std::cout << "Launcher finnished running" << std::endl;
    OutputPrinter::GraphToVtk(outputpath,G);

}

void TestAsKNew(std::string inputpath, std::string path, int number, double param, double str)
{
    std::string outputpath = path + "InitialVTK/output" + std::to_string(number) + ".vtk";
    arma::mat cloud;
    MyGraphType G;
    InputReader::XYZtoMAT(inputpath,cloud);
    AsKLauncher launcher(param, str, 1.1, "AsK",false, true);
 //   launcher.SetPostProcessing(folder + "cloud.npy", folder + "abstractGraph.npy", folder + "Edges.npy", folder + "Belonghood.npy");
    cloud = cloud.t();
    std::string folder = path;
    launcher.SetPostProcessing(folder + "Points/cloud" + std::to_string(number) + ".npy" ,
     folder + "Graph/Graph" + std::to_string(number) + ".npy",
     folder + "Edges/Edge " + std::to_string(number) + ".npy",
     folder + "Segmentindices/si"+ std::to_string(number)  + ".npy",
     folder + "Pointindices/pi" + std::to_string(number) + ".npy");

    launcher.Run(cloud, G);

    std::cout << "outputpath " << outputpath << std::endl;
    OutputPrinter::GraphToVtk(outputpath,G);

}


void MassiveTestAsk(std::string pathName, std::string outputpathName)
{
    std::vector<boost::filesystem::directory_entry> directories;
    boost::filesystem::path p(pathName);
    std::copy(boost::filesystem::directory_iterator(p), boost::filesystem::directory_iterator(), back_inserter(directories));
    std::sort(directories.begin(), directories.end());
    for (int i = 0; i < directories.size(); i++)
    {

        std::cout << "case " << i << std::endl;
        TestAsKNew(directories[i].path().string(), outputpathName, i , 20.0, 1.05);

     //   std::cout << directories[i].path().string() << std::endl;
    }

}

int main()
{
std::cout << "Test " << std::endl;

}

