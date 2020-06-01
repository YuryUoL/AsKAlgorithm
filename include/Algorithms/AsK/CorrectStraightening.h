#ifndef CORRECTSTRAIGHTENING_H
#define CORRECTSTRAIGHTENING_H


#include "Graph/Graph.h"
#include <mlpack/core.hpp>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include "Algorithms/Distances/PointSegmentDistance.h"
#include "Algorithms/Convertors/OutputPrinter.h"


class CorrectStraightening
{
public:
    CorrectStraightening() {}
    //! First we do the allocation



    static int GetOnlyCorrectNeighbor(MyGraphType & G, int kk, std::vector<bool> & coloring)
    {
        auto vpair = boost::adjacent_vertices(kk,G);

        for (auto it = vpair.first; it != vpair.second ; it++)
        {
            if (coloring[(*it)] == false)
            {
                return (*it);
            }

        }
        //return (*(iter.first));


    }

    static void PathCreator(MyGraphType & G, arma::mat & cloud, std::vector<std::vector<int>> & resultingPaths, std::vector<int> & CloudToGraph, std::vector<int> & GraphToCloud, double & MaxDistance, std::vector<int> & correspondance)
    {

        CorrectStraightening::PointsToSegments( G, cloud, CloudToGraph, MaxDistance,correspondance);

        std::vector<std::vector<arma::mat>> edges;
        std::vector<std::vector<int>> paths;
        std::vector<int> GraphToMat(boost::num_vertices(G));
        std::deque<int> la;
        std::vector<bool> coloring(boost::num_vertices(G));

        for (int i = 0; i < coloring.size(); i++)
        {
            coloring[i] = false;

        }
        for (int i = 0; i < boost::num_vertices(G); i++)
        {
            if ((boost::degree(i,G) != 2)&&(boost::degree(i,G) != 0 ))
            {
                la.push_back(i);
                break;
            }
        }
        while(!la.empty())
        {
            int iter = la.front();
            la.pop_front();
            //paths.push_back(std::vector<int>{iter})
            auto vpair = boost::adjacent_vertices(iter,G);
            coloring[iter] = true;
            for (auto it = vpair.first; it != vpair.second; it++)
            {
                if (coloring[*it] == 1)
                {
                    continue;
                }
                paths.push_back(std::vector<int> {iter});
                int stupidcount = 0;
                int flyingSheep = *it;
                coloring[flyingSheep] = true;
                while (boost::degree(flyingSheep, G) == 2)
                {
                    if (stupidcount < 100)
                    {
                        stupidcount++;
                    }
                    //!  std::cout << "In loop " << flyingSheep << "degree: " <<  boost::degree(flyingSheep, G) << std::endl;
                    paths[paths.size()-1].push_back(flyingSheep);
                    //! flyingSheep --- //
                    coloring[flyingSheep] = true;
                    flyingSheep = GetOnlyCorrectNeighbor(G, flyingSheep, coloring);
                }
                paths[paths.size()-1].push_back(flyingSheep);
                if (boost::degree(flyingSheep, G) != 1)
                {
                    la.push_back(flyingSheep);
                }

            }
        }
        //resultingPaths.resize(paths.size());



        for(int i = 0; i < paths.size(); i++)
        {
            resultingPaths.push_back(std::vector<int>());
            resultingPaths[resultingPaths.size()-1].push_back(GraphToCloud[paths[i][0]]);

            for (int j = 0 ; j < paths[i].size() - 1; j++)
            {
                arma::mat segment(cloud.n_rows, 2);
                segment.col(0) = G[paths[i][j]].p;
                segment.col(1) = G[paths[i][j+1]].p;
                auto eitt = boost::edge(paths[i][j], paths[i][j+1], G).first;
                std::sort(G[eitt].indices.begin(), G[eitt].indices.end(), [&segment, &cloud](int a, int b)
                {
                    double t = arma::dot(cloud.col(a)-segment.col(0),segment.col(1)-segment.col(0)) / (arma::dot(segment.col(1)-segment.col(0),segment.col(1)-segment.col(0)));
                    double t2 = arma::dot(cloud.col(b)-segment.col(0),segment.col(1)-segment.col(0)) / (arma::dot(segment.col(1)-segment.col(0),segment.col(1)-segment.col(0)));
                    return t < t2;
                });

                for (int xx = 0; xx < G[eitt].indices.size(); xx++)
                {
                    resultingPaths[resultingPaths.size()-1].push_back((G[eitt].indices)[xx]);
                }

            }
            int tmpSize = paths[i].size();
            resultingPaths[resultingPaths.size()-1].push_back(GraphToCloud[paths[i][tmpSize-1]]);


        }
    }

    static void PointsToSegments(MyGraphType & G, arma::mat & cloud, std::vector<int> & CloudToGraph, double & maxDistance, std::vector<int> & correspondance)
    {

        arma::mat referenceSet;
        MatToGraph::Reverse(referenceSet, G);
        mlpack::neighbor::KNN a(referenceSet);
        arma::Mat<size_t> resultingNeighbors;
        arma::mat resultingDistances;
        a.Search(cloud, 1, resultingNeighbors, resultingDistances);


        for (int i = 0; i < correspondance.size(); i++)
        {
            int iter = resultingNeighbors(0,i);
          //  int planetEarth = CloudToGraph[iter];
            int planetEarth = iter;
            //  int iter = correspondance[i];


            if (CloudToGraph[i] != -1)
            {
                continue;
                //if (boost::out_degree(CloudToGraph[i], G) != 2)
                //{
                // marked[i] = true;
                //}
            }
            //std::cout << "i: " << i << " | correspondance[i]: " << iter << " | CloudToGraph[y]:" << planetEarth << std::endl;
            //std::cout << "i: " << i << " | CloudToGraph[i]: " << CloudToGraph[i] <<  std::endl;


            auto edges = boost::out_edges(planetEarth,G);
            double minDz = std::numeric_limits<double>::max();
            edge_descriptor found;
            int foundsource = -1;
            int foundtarget =  -1;
            for (auto eit = edges.first ; eit != edges.second ; eit++)
            {
                double ddd = sqrt(PointSegmentDistance::Calculate(G[(*eit)].segment, cloud.col(i)));
                if (ddd < minDz)
                {
                    minDz = ddd;
                    found = *eit;
                    foundsource = boost::source(*eit,G);
                    foundtarget = boost::target(*eit,G);
                }
            }

            G[found].indices.push_back(i);
            maxDistance = std::max(minDz, maxDistance);

        }

        auto vpair = boost::edges(G);
        int sum = 0;
        int numm = 0;
        for (auto it = vpair.first ; it != vpair.second ; it++)
        {
            numm++;
            sum = sum + G[*it].indices.size();
        }

    }


    static double maximalDistance(int i, int j, std::vector<int> & points, arma::mat & cloud)
    {
        if (j > points.size() - 1)
        {
            return std::numeric_limits<double>::max();
            //return 0;
        }
        arma::mat iterSegment(cloud.n_rows, 2);
        //Segment s = Segment(points[i],points[j]);
        iterSegment.col(0) = cloud.col(points[i]);
        iterSegment.col(1) = cloud.col(points[j]);

        double maxD = 0;
        for (int k = i; k <= j; k++)
        {
            maxD = std::max(maxD, sqrt(PointSegmentDistance::Calculate(iterSegment,cloud.col(points[k]))));
        }
        return maxD;
    }


    static bool RunTheAlgorithm(arma::mat & cloud, std::vector<int> & points, std::vector<int> & thepath, double e)
    {
        thepath.push_back(points[0]);
        int sizepath = points.size();
        int a = 0;
        int b = 0;
        while (b < points.size() - 1)
        {
            int l = 0;
            //  l = l + 1; //! Is this required? Most likely now
            while((CorrectStraightening::maximalDistance(b, b + pow(2,(l+1)), points, cloud) <= e)
                    &&(b + pow(2,(l+1)) < points.size()))
            {
                l++;
            }
            int low = pow(2,l);
            int sizesizesize = points.size();
            int powerpowerpower = pow(2,l+1);
            int high = std::min(powerpowerpower,sizesizesize - b);
            while(low < high - 1)
            {
                int mid = (low + high)/2;
                int sumdxd = b + mid;
                double maxd = CorrectStraightening::maximalDistance(b, b + mid, points,cloud);
                if (maxd <= e)
                {
                    low = mid;
                }
                else
                {
                    high = mid;
                }
            }
            b = b + low;
            if (b < points.size())
            {
                thepath.push_back(points[b]);
            }
            else
            {
                thepath.push_back(points[points.size()-1]);
            }
            a = a + 1;
        }
        return true;
    }


    static double ComputeStraightening(MyGraphType & Simplified, arma::mat & cloud, std::vector<std::vector<int>> & out, double e, std::vector<int> & CloudToGraph, std::vector<int> & GraphToCloud, std::vector<int> & correspondance)
    {
        // OutputPrinter::DebugGraph("/home/yury/LocalTests4/MicelleOutputsWorm/graphDebug.txt", Simplified);

        double MaxDistance = 0;
        std::vector<std::vector<int>> tmp;
        CorrectStraightening::PathCreator(Simplified, cloud, tmp, CloudToGraph, GraphToCloud, MaxDistance, correspondance);
        //OutputPrinter::InfoToCsv(std::string path, std::vector<int> & indices, arma::mat & cloud)


        out.resize(tmp.size());
        for (int i = 0; i < tmp.size(); i++)
        {
            CorrectStraightening::RunTheAlgorithm(cloud, tmp[i], out[i], e * MaxDistance);
        }
        return MaxDistance;


    }

protected:

private:
};

#endif // CORRECTSTRAIGHTENING_H
