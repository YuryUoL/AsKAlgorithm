#include "BranchDetection.h"
#include "Algorithms/Convertors/OutputPrinter.h"

BranchDetection::BranchDetection()
{

}



bool CompareD(std::pair<int,double> & a, std::pair<int, double> & b)
{
    return a.second > b.second;
}



edge_descriptor getOnlyEdge(vertex_descriptor v, MyGraphType & G)
{

    auto containerEdgeN = out_edges(v, G);
    auto pointerEdgeN = containerEdgeN.first;
    edge_descriptor edgebetween = *pointerEdgeN;
    return edgebetween;


}

void BranchDetection::CorersponanceMarker( std::vector<bool> & coloring, std::vector<int> & next, std::vector<int> & edgePoints)
{
    for (int i = 0; i < edgePoints.size(); i++)
    {
        int j = edgePoints[i];
        std::deque<int> elems;
        while (coloring[j] == false )
        {
            elems.push_back(j);
            j = next[j];
        }
        for (auto it : elems)
        {
            next[it] = j;
        }

    }

}

void BranchDetection::branchIndexation(MyGraphType & G, double thresh, MyGraphType & out, std::vector<int> & OldToNew, std::vector<int> & NewToOld, std::vector<int> & next)
{
    next.resize(boost::num_vertices(G));
    for (int i = 0; i < next.size(); i++)
    {
        next[i] = i;
    }

//! STEP ONE:
    MyGraphType tmp = G;
    std::vector<int> DegreeZeroVertices(boost::num_vertices(G));
    std::vector<bool> branchedCandidate(boost::num_vertices(G));
    std::vector<bool> boundary(boost::num_vertices(G));
    std::vector<double> PotentialValue(boost::num_vertices(G));
    std::vector<double> maximumLength(boost::num_vertices(G));
    std::vector<int> prev(boost::num_vertices(G));
    std::vector<int> prev2(boost::num_vertices(G));
    std::vector<std::pair<int,double>> vertexContainer;
    std::vector<std::vector<std::pair<int,double>>> Values(boost::num_vertices(G), std::vector<std::pair<int,double>>());
    int specialCase = -1;
    //std::vector<int> next(boost::num_vertices(G));
    auto vpair = vertices(G);
    for (auto iter = vpair.first; iter != vpair.second; iter++)
    {
        int sizedges = out_degree(*iter,G);
        if (sizedges == 1)
        {
            auto container = out_edges(*iter, G);
            auto pointer = container.first;
            edge_descriptor e = *pointer;
            double d = G[e].distance;
            PotentialValue[*iter] = d;
            boundary[*iter] = true;
            vertexContainer.push_back(std::make_pair(*iter,d));
        }
        else if (sizedges >= 3)
        {
            branchedCandidate[*iter] = true;
        }
    }
    VH heap = VH(CompareD,vertexContainer);


    while(!heap.empty())
    {
        std::pair<int,double> cur = heap.top();
        heap.pop();
        vertex_descriptor it = cur.first;

        if (boost::out_degree(it,tmp) == 0)
        {
            continue;
        }

        edge_descriptor edgebetween = getOnlyEdge(it, tmp);


        auto container = boost::adjacent_vertices(it, tmp);
        auto pointer = container.first;

        vertex_descriptor neighbor = *pointer;

        boost::remove_edge(edgebetween, tmp);

        int degree = out_degree(neighbor,tmp);

        //double pot = PotentialValue[it];
        next[it] = neighbor;

      //  if (degree != 0)
       // {
        Values[neighbor].push_back(std::make_pair(it,cur.second));

          //  kwriter << "Reached degree 0 !!!!!! " << " | Culprit: " << it << " Its neighbor: " <<  neighbor << std::endl;

        //}
        if (degree == 1)
        {
            edge_descriptor e = getOnlyEdge(neighbor, tmp);

            double itsdistance = tmp[e].distance;
            double pot = cur.second;
            double totaldistance = itsdistance + pot;

          //  kwriter << "Adding to heap " << " Guy: " << it << " Its neighbor: " <<  neighbor << " | Totaldistance: " << totaldistance << std::endl;
            heap.push(std::make_pair(neighbor, totaldistance));
        }
    }
//! Step two:
    //static void PrintThreshedOnly(std::string path, MyGraphType & G, std::vector<double> & threshes, double & thresh);
    std::deque<int> indexContainer;
    std::vector<bool> coloring(boost::num_vertices(G));
    std::vector<bool> accepted(boost::num_vertices(G));
    for (int i = 0; i < boost::num_vertices(G); i++)
    {

        if (Values[i].size() > 1)
        {
            if(Values[i][Values[i].size()-2].second > thresh)
                {
                    accepted[i] = true;
                    coloring[i] = true;
                    indexContainer.push_back(i);

                }
            }

    }
  //  OutputPrinter::PrintThreshedOnly("/home/yury/LocalTests4/ChristmasDebug/NewThreshes/threshZero.csv", G, PotentialValue, 0.0, coloring);
    while(!indexContainer.empty())
    {
        int iter = indexContainer.front();
        indexContainer.pop_front();
        bool activator = true;
        for (int j = Values[iter].size() - 1; j > -1; j--)
        {
            double value = Values[iter][j].second;
            if (value < thresh)
            {
                break;
            }
            else
            {
                activator = false;

                int goodguy =  Values[iter][j].first;
                Graph::add_edge(tmp, iter,goodguy);
                if ((coloring[goodguy] != true)&&(boost::degree(goodguy,G) > 1))
                {
                    indexContainer.push_back(Values[iter][j].first);
                }

            }
        }
        if (activator)
        {
            int guest = Values[iter][Values[iter].size()-1].first;
            Graph::add_edge(tmp, iter, guest);
            if ((coloring[guest] != true)&&(boost::degree(guest,G) > 1))
            {
                indexContainer.push_back(guest);
            }
        }


    }
  //  OutputPrinter::GraphToVtk("/home/yury/LocalTests4/ChristmasDebug/intermediateExperimental.vtk", tmp);


  //  OutputPrinter::PrintThreshedOnly("/home/yury/LocalTests4/ChristmasDebug/NewThreshes/threshFirst.csv", G, PotentialValue, 0.0, coloring);

//! Step three:
//!  std::vector<int> & CloudToGraph
    OldToNew.resize(boost::num_vertices(G));
    std::fill(OldToNew.begin(), OldToNew.end(), -1);
    //OldToNew.fill(-1);

    int numColoring = 0;
    for (int i = 0; i < coloring.size(); i++)
    {
        if (boost::degree(i,tmp) != 0)
        {
            coloring[i] = true;
        }
        if (coloring[i] == true)
        {
            numColoring++;
        }

    }

   // OutputPrinter::PrintThreshedOnly("/home/yury/LocalTests4/ChristmasDebug/NewThreshes/threshSecond.csv", G, PotentialValue, 0.0, coloring);

    NewToOld.resize(numColoring);


    for (int i = 0; i < boost::num_vertices(G); i++)
    {
        if ((boost::degree(i,tmp) > 0))
        {
            int neww = Graph::add_vertex(out,G[i].p);
            NewToOld[neww] = i;
            OldToNew[i] = neww;
        }
    }


    auto edges = boost::edges(tmp);
    for (auto edgeit = edges.first; edgeit != edges.second; edgeit++)
    {
        int xDxD = boost::source(*edgeit, tmp);
        int firstV = OldToNew[xDxD];
        int secondV = OldToNew[boost::target(*edgeit,tmp)];


        Graph::add_edge(out,firstV,secondV);
    }


    //! void BranchDetection::CorersponanceMarker( std::vector<bool> & coloring, std::vector<int> & next, std::vector<int> & edgePoints)
    std::vector<int> edgePoints;
    for (int i = 0; i < boost::num_vertices(G); i++)
    {
        if (boost::degree(i, G) == 1)
        {
            edgePoints.push_back(i);
        }

    }
    //  std::cout << "Before correspondanceMarker" << std::endl;

    for (int i = 0; i < coloring.size(); i++)
    {
        if (coloring[i] == true)
        {
            next[i] = i;

        }

    }

    BranchDetection::CorersponanceMarker( coloring, next, edgePoints);


    //OutputPrinter::DebugArray("/home/yury/LocalTests4/MicelleOutputsWorm/nextStructure.txt", next,  coloring, G);

    //! void BranchDetection::CorersponanceMarker( std::vector<bool> & coloring, std::vector<int> & next, std::vector<int> & edgePoints)


//   std::cout << "After correspondanceMarker" << std::endl;
    //  for (int i = 0; i < next.size(); i++)
    //  {
    //     std::cout << "i: " << i << " : " << next[i] << std::endl;
    //  }


}


