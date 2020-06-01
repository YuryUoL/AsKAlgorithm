#include "OutputPrinter.h"


OutputPrinter::OutputPrinter()
{
    //ctor
}

void OutputPrinter::PrintThreshedOnly(std::string path, MyGraphType & G, std::vector<double> & threshes, double thresh, std::vector<bool> & branchedCandidate)
{
    std::ofstream mystream;
    mystream.open(path);
    mystream << "x, y, z" << std::endl;
    auto vertices = boost::vertices(G);
    for (auto it = vertices.first; it != vertices.second; it++)
    {
        if ((branchedCandidate[*it] == true))
        {
            for (int j = 0; j < 3; j++)
            {
                mystream << G[(*it)].p(j) << ", ";
            }
              mystream << std::endl;

        }
        //  double k = (double) threshes[*it];
        // mystream << k << std::endl;



    }




    mystream.close();

}

void OutputPrinter::ThreshColored(std::string path, MyGraphType & G, std::vector<double> & threshes)
{
    std::ofstream mystream;
    mystream.open(path);
    mystream << "x, y, z, t" << std::endl;
    auto vertices = boost::vertices(G);
    for (auto it = vertices.first; it != vertices.second; it++)
    {
        for (int j = 0; j < 3; j++)
        {
            mystream << G[(*it)].p(j) << ", ";
        }
        double k = (double) threshes[*it];
        mystream << k << std::endl;



    }
    mystream.close();



}

void OutputPrinter::InfoToCsv(std::string path, std::vector<int> & indices, arma::mat & cloud)
{
    std::ofstream mystream;
    mystream.open(path);
    mystream << "x, y, z, t" << std::endl;
    for (int i = 0 ; i < indices.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            mystream << cloud(j,indices[i]);

            mystream << ", ";

        }
        double k = (double) i / (double) indices.size();
        mystream << k << std::endl;



    }
    mystream.close();




}

void OutputPrinter::GraphToVtk(std::string path, MyGraphType & G)
{
    std::ofstream mystream;
    mystream.open(path);
    mystream << "# vtk DataFile Version 1.0\n";
    mystream << "3D triangulation data\n";
    mystream << "ASCII\n";
    mystream << std::endl;
    mystream << "DATASET POLYDATA\n";

    mystream << "POINTS " << num_vertices(G) << " float\n";
    for(int i=0; i<num_vertices(G); i++)
    {
        for (int j = 0 ; j < G[i].p.n_cols; j++)
        {
            for(int k = 0; k < G[i].p.n_rows; k++)
            {
                mystream << G[i].p(k,j) << " ";

            }
            mystream << std::endl;
        }
        //mystream << G[i].p << std::endl;
    }

    mystream << "LINES " << (num_edges(G)) << " " << (num_edges(G))*3 << std::endl;

    auto epair = edges(G);
    for(auto iter=epair.first; iter!=epair.second; iter++)
    {
        mystream  << "2 " << source(*iter, G) << " " << target(*iter, G) <<std::endl;
    }
    mystream.close();
}

void OutputPrinter::DebugGraph(std::string path, MyGraphType &Simplified)
{
    std::ofstream mystream;
    mystream.open(path);
    auto vpairxD = boost::vertices(Simplified);
    for (auto it = vpairxD.first; it != vpairxD.second; it++)
    {
        mystream << *it << " : degree : " << boost::degree(*it,Simplified) << std::endl;
    }
    mystream.close();

}

void OutputPrinter::DebugArray(std::string path, std::vector<int> & correspondance, std::vector<bool> & coloring, MyGraphType & G)
{
    std::ofstream mystream;
    mystream.open(path);
    for (int i = 0; i < correspondance.size(); i++)
    {
        mystream << " i " << i << " : value : " << correspondance[i] << " : coloring : " << coloring[i] << " degree: " << boost::degree(i,G)<< std::endl;
    }

    mystream.close();
}


//OutputPrinter::~OutputPrinter()
//{
//dtor
//}
