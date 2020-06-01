#ifndef OUTPUTPRINTER_H
#define OUTPUTPRINTER_H

#include "Graph/Graph.h"


class OutputPrinter
{
    public:
        OutputPrinter();
        static void ThreshColored(std::string path, MyGraphType & G, std::vector<double> & threshes);
        static void GraphToVtk(std::string path, MyGraphType &G);
        static void DebugGraph(std::string path, MyGraphType &G);
        static void DebugArray(std::string path, std::vector<int> & correspondance, std::vector<bool> & coloring,  MyGraphType & G);
        static void InfoToCsv(std::string path, std::vector<int> & indices, arma::mat & cloud);
        static void PrintThreshedOnly(std::string path, MyGraphType & G, std::vector<double> & threshes, double thresh, std::vector<bool> & branchedCandidate);



        //virtual ~OutputPrinter();

    protected:

    private:
};

#endif // OUTPUTPRINTER_H
