#include <iostream>
#include "Algorithms/Convertors/OutputPrinter.h"
#include "Algorithms/Convertors/InputReader.h"
#include "Algorithms/AsK/AsKLauncher.h"

using namespace std;


void TestAsK(std::string outputpath, std::string inputpath, double param = 20.0)
{
arma::mat cloud;
MyGraphType G;
InputReader::XYZtoMAT(inputpath,cloud);
AsKLauncher launcher(param, 1.1, 1.1, "AsK");
cloud = cloud.t();
launcher.Run(cloud, G);
OutputPrinter::GraphToVtk(outputpath,G);

}

int main()
{

    std::string cloudpath = "";
    std::string outputpath = "";
    TestAsK(outputpath,cloudpath , 30.0);


}
