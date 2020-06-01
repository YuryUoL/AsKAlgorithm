#include "InputReader.h"
#include <fstream>
#include "Graph/Graph.h"

InputReader::InputReader()
{
    //ctor
}

void InputReader::XYZtoMAT(std::string fl, arma::mat & data)
{
    // We make sure that the data-matrix is empy
    data.reset();
    std::ifstream infile(fl.c_str());
    double x;
    double y;
    double z;
    std::string w;
    int k;
    infile >> k;
    data.set_size(k,3);
    getline(infile,w);
    getline(infile,w);
    //  std::vector<std::vector<double>> toArma;
    int i = 0;
    while (infile >> w >> x >> y >> z)
    {
        data(i,0) = x;
        data(i,1) = y;
        data(i,2) = z;
        i++;
    }
}

//InputReader::~InputReader()
//{
    //dtor
//}
