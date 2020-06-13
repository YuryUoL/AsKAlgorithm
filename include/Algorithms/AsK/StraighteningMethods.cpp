#include "StraighteningMethods.h"
#include "Computation.h"
#include "GeneralConvertor.h"

StraighteningMethods::StraighteningMethods()
{
    //ctor
}

double maximalDistance(int i, int j, Segment l, std::vector<PointPlus> & path )
{
    double maxdistance = 0;
    for (int k = i; k <= j; k++)
    {
        for(auto it = path[k].second.begin(); it != path[k].second.end(); it++)
        {
            double distance = sqrt(CGAL::squared_distance((*it), l));
            maxdistance = std::max(distance, maxdistance);
        }
    }
    return maxdistance;
}

double maximalDistance(int i, int j, std::vector<PointPlus> & path )
{
    int bumbum = path.size()-1;
    int jp = std::min(j,bumbum);
    return maximalDistance(i,jp, Segment(path[i].first, path[j].first), path);
}


bool StraighteningMethods::OptimizeSingle(std::vector<PointPlus> & path,  std::list<Point> & thepath, double e)
{
    int sizepath = path.size();
    int a = 0;
    int b = 0;
    while (b < path.size() - 1)
    {
        int l = 0;
        bool firsttest = maximalDistance(b, b+1, path) <= e;
        if (firsttest == false)
        {
            std::cout << "STRAIGHTENING UNSUCCESFUL!!!!!!!!!!!!!!!!" << std::endl;
           // throw std::exception("You win, but you lose, game over");
            return false;
        }
        l = l + 1; //! Is this required?
        while((maximalDistance(b, b + pow(2,(l+1)), path) <= e)
                &&(b + pow(2,(l+1)) < path.size()))
        {
            l++;
        }
        int low = pow(2,l);
        int sizesizesize = path.size();
        int powerpowerpower = pow(2,l+1);
        int high = std::min(powerpowerpower,sizesizesize - b);
        while(low < high - 1)
        {
            int mid = (low + high)/2;
            int sumdxd = b + mid;
            double maxd = maximalDistance(b, b + mid, path);
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
        if (b < path.size())
        {
            thepath.push_back(path[b].first);
        }
        else
        {
            thepath.push_back(path[path.size()-1].first);
        }
        a = a + 1;
    }
    return true;
}

void StraighteningMethods::Optimize(std::list<std::list<Point>> & out, std::vector<std::vector<PointPlus>> & paths, double e)
{
    for (int i = 0; i < paths.size(); i++)
    {
        std::list<Point> thepath;
        thepath.push_back(paths[i][0].first);
        StraighteningMethods::OptimizeSingle(paths[i],thepath,e);
        out.push_back(thepath);
    }

}

void addToPath(std::vector<PointPlus> & path, MyGraphType & G, vertex_descriptor v)
{
    path.push_back(std::make_pair(G[v].p, std::list<Point>()));
}

void StraighteningMethods::GraphToPaths(std::vector<std::vector<PointPlus>> & paths, MyGraphType & G)
{
// Find vertex of degree non two.
    std::deque<std::pair<vertex_descriptor,vertex_descriptor>> ddd;
    auto vpair = boost::vertices(G);
    for (auto it = vpair.first; it != vpair.second ; it++)
    {
        if (boost::degree(*it, G) != 2)
        {
            ddd.push_back(std::make_pair(*it,*it));
            break;
        }
    }
    while (!ddd.empty())
    {
        vertex_descriptor a = ddd.front().first;
        vertex_descriptor wrong = ddd.front().second;

        //     std::cout << a << std::endl;
        ddd.pop_front();
        auto epair = boost::adjacent_vertices(a,G);
        for (auto it = epair.first ; it != epair.second; it++)
        {
            //  std::cout << *it << std::endl;
            if (*it == wrong)
            {
                continue;
            }
            std::vector<PointPlus> path;
            addToPath(path,G,a);
            vertex_descriptor ot = *it;
            vertex_descriptor prev = a;
            while (boost::degree(ot,G) == 2)
            {
                //    std::cout << "InnerLoop: " << ot << std::endl;
                auto innerpair = boost::adjacent_vertices(ot,G);
                for (auto bt = innerpair.first ; bt != innerpair.second; bt++)
                {
                    if (*bt != prev)
                    {
                        prev = ot;
                        ot = *bt;
                        break;
                    }
                }
                addToPath(path,G,ot);

            }
            addToPath(path,G,ot);
            paths.push_back(path);
            if (boost::degree(ot,G) != 1)
            {
                ddd.push_back(std::make_pair(ot, prev));
            }
        }
    }

}

void fairsplit(std::vector<PointPlus> & path, int i, int j)
{
    std::list<Point> container = path[i].second;
    path[i].second.clear();

    for (auto it = container.begin(); it != container.end(); it++)
    {
        double di = sqrt(CGAL::squared_distance(path[i].first, *it));
        double dj = sqrt(CGAL::squared_distance(path[j].first, *it));
        if (di >= dj)
        {
            path[j].second.push_back(*it);
        }
        else
        {
            path[i].second.push_back(*it);
        }
    }

}


void StraighteningMethods::Allocate(std::vector<std::vector<PointPlus>> & paths, MyGraphType & G, std::list<Point> & cloud)
{

    Triangulation T;
    std::map<Point, std::pair<int, int>> generalmap;
    for (int i = 0; i < paths.size(); i++)
    {
        for(int j = 1; j < paths[i].size()-1; j++)
        {
            Point pt = paths[i][j].first;
            generalmap[pt] = std::make_pair(i,j);
            T.insert(pt);
        }
    }
    for (auto it = cloud.begin(); it != cloud.end(); it++)
    {
        Point result = (T.nearest_vertex(*it))->point();
        auto pair = generalmap[result];
        paths[pair.first][pair.second].second.push_back(*it);

    }
    for (int i = 0; i < paths.size(); i++)
    {
        if (paths[i].size() != 1)
        {
            fairsplit(paths[i], 1,0 );
            fairsplit(paths[i], paths[i].size()-2, paths[i].size()-1);
        }

    }
}

void PathGrabber(std::vector<std::vector<PointPlus>> & pp, std::list<std::list<Point>> & bt)
{

    for (auto it = pp.begin(); it != pp.end(); it++)
    {
        std::list<Point> aatt;
        for (auto jt = (*it).begin(); jt != (*it).end(); jt++)
        {
            Point q = (*jt).first;
            aatt.push_back(q);
        }
        bt.push_back(aatt);

    }

}

double ComputeOriginalDistance(std::list<Point> & cloud, std::vector<std::vector<PointPlus>> & pp )
{
    std::list<std::list<Point>> bt;
    PathGrabber(pp, bt);
    return Computation::AABBDistance(bt,cloud);

}
//  data.set_size(3,k);
//          data(0,j) = it->x();

//void StraighteningMethods::NumBranchedPoints(MyGraphType & G)
//{
//
//
//
//}

void StraighteningMethods::ConvertGraphToArmaMat(arma::mat & pathpoints, MyGraphType & G)
{
    size_t numVertices = boost::num_vertices(G);
    pathpoints.set_size(3,numVertices);
    auto vpair = boost::vertices(G);
    for (auto it = vpair.first; it != vpair.second; it++)
    {
        pathpoints(0,*it) = G[*it].p.x();
        pathpoints(1,*it) = G[*it].p.y();
        pathpoints(2,*it) = G[*it].p.z();
    }
//    int globalcounter = 0;
//    for (int i = 0 ; i < paths.size(); i++)
//    {
//        for (int j = 0; j < paths[i].size(); j++)
//        {
//            pathpoints()
//            globalcounter++;
//        }
//    }
}

void StraighteningMethods::DualTreeAllocation(MyGraphType & G, std::list<Point> & cloud, arma::Mat<size_t> & theNeighbors)
{
//! This method provides us with linear time allocation for the cloud:

    arma::mat referenceSet;
    arma::mat querySet;
    StraighteningMethods::ConvertGraphToArmaMat(referenceSet,G);
    querySet.set_size(3,cloud.size());
    GeneralConvertor::ListToMatTransposed(cloud,querySet);
    //  arma::mat crf = referenceSet.t();
//    arma::mat qrf = querySet.t();
    DualTreeComputation::NearestNeighborForTwoKDTrees(theNeighbors,referenceSet,querySet);
    //  DualTreeComputation::NearestNeighborForTwoKDTrees(theNeighbors,crf,qrf);
//DualTreeComputation::NearestNeighborForTwoKcurDistanaceDTrees(arma::Mat<size_t> & results,arma::mat & referenceData, arma::mat & queryData)
}       // minvalue = minvalueFinder(std::get<0>(data), scale, settings);


void StraighteningMethods::DualTreeAllocatorPostProcessor(MyGraphType & G, arma::Mat<size_t> & theNeighbors, std::list<Point> & clouds, std::vector<std::list<Point>> & pAllocation)
{
    auto vpair = boost::vertices(G);
    std::list<vertex_descriptor> theBadBoys;

    //! Allocate stuff properly:
    int indx = 0;
    for (auto clit = clouds.begin(); clit != clouds.end(); clit++)
    {
        size_t corIndx = theNeighbors(0,indx);
        pAllocation[corIndx].push_back(*clit);
        indx++;
    }

    //! Find the bad boys:
    for (auto it = vpair.first; it != vpair.second ; it++)
    {
        if (boost::degree(*it,G) > 2)
        {
            theBadBoys.push_back(*it);
        }
    }
    //! Process the bad boys
    for (auto vit = theBadBoys.begin(); vit != theBadBoys.end(); vit++)
    {
        for (auto clit = pAllocation[*vit].begin(); clit != pAllocation[*vit].end(); clit++)
        {
            auto bpair = boost::adjacent_vertices(*vit,G);
            vertex_descriptor kw;
            double curDistance = std::numeric_limits<double>::max();
            for (auto bit = bpair.first; bit != bpair.second; bit++)
            {
                double tmpdistance = sqrt(CGAL::squared_distance(*clit, G[*vit].p));
                if (tmpdistance < curDistance)
                {
                    curDistance = tmpdistance;
                    kw = *bit;
                }
            }
            pAllocation[kw].push_back(*clit);
        }
        pAllocation[*vit].clear();
    }

}

void StraighteningMethods::addToPathNew(std::vector<PointPlus> & path, MyGraphType & G, vertex_descriptor v,std::vector<std::list<Point>> & pAllocation)
{
    path.push_back(std::make_pair(G[v].p, pAllocation[v]));
}

void StraighteningMethods::GraphToPathsNew(std::vector<std::vector<PointPlus>> & paths, MyGraphType & G, std::vector<std::list<Point>> & pAllocation)
{
// Find vertex of degree non two.
    std::deque<std::pair<vertex_descriptor,vertex_descriptor>> ddd;
    auto vpair = boost::vertices(G);
    for (auto it = vpair.first; it != vpair.second ; it++)
    {
        if (boost::degree(*it, G) != 2)
        {
            ddd.push_back(std::make_pair(*it,*it));
            break;
        }
    }
    while (!ddd.empty())
    {
        vertex_descriptor a = ddd.front().first;
        vertex_descriptor wrong = ddd.front().second;

        //     std::cout << a << std::endl;
        ddd.pop_front();
        auto epair = boost::adjacent_vertices(a,G);
        for (auto it = epair.first ; it != epair.second; it++)
        {
            //  std::cout << *it << std::endl;
            if (*it == wrong)
            {
                continue;
            }
            std::vector<PointPlus> path;
            StraighteningMethods::addToPathNew(path,G,a,pAllocation);
            vertex_descriptor ot = *it;
            vertex_descriptor prev = a;
            while (boost::degree(ot,G) == 2)
            {
                //    std::cout << "InnerLoop: " << ot << std::endl;
                auto innerpair = boost::adjacent_vertices(ot,G);
                for (auto bt = innerpair.first ; bt != innerpair.second; bt++)
                {
                    if (*bt != prev)
                    {
                        prev = ot;
                        ot = *bt;
                        break;
                    }
                }
                StraighteningMethods::addToPathNew(path,G,ot,pAllocation);

            }
            StraighteningMethods::addToPathNew(path,G,ot, pAllocation);
            paths.push_back(path);
            if (boost::degree(ot,G) != 1)
            {
                ddd.push_back(std::make_pair(ot, prev));
            }
        }
    }

}


double StraighteningMethods::UpgradedStraightening(MyGraphType & G, std::list<Point> & cloud, std::list<std::list<Point>> & out, double e, arma::Mat<size_t> & theNeighbors)
{
    std::vector<std::vector<PointPlus>> paths;
    std::vector<std::pair<bool,bool>> status;
    std::vector<std::list<Point>> pAllocation(boost::num_vertices(G), std::list<Point>());
    StraighteningMethods::DualTreeAllocation(G, cloud,  theNeighbors);
    StraighteningMethods::DualTreeAllocatorPostProcessor(G, theNeighbors, cloud,pAllocation);
    StraighteningMethods::GraphToPathsNew(paths, G, pAllocation);

//! Repair Allocation // Actually this wont be done, because we are running out of time, it would be a complete nonsense.

    for (int i = 0; i < paths.size(); i++)
    {
        if (paths[i].size() != 1)
        {
            fairsplit(paths[i],1,0);
            fairsplit(paths[i], paths[i].size()-2, paths[i].size()-1);
        }
    }
    double ddd = ComputeOriginalDistance(cloud,paths);
  //  std::cout << "Parameter ddd: " << ddd << std::endl;
    StraighteningMethods::Optimize(out, paths, ddd * e);
//    std::cout << "Parameter ddd*e: " << ddd * e << std::endl;

//    int ctr= 0;
//    for (auto it = out.begin(); it != out.end(); it++)
//    {
//        std::cout << "(Upgraded straightening) Path: " << ctr << " size: " << (*it).size() << std::endl;
//        ctr++;
//    }

    return ddd;



}

double StraighteningMethods::ClassicStraightening(MyGraphType & G, std::list<Point> & cloud, std::list<std::list<Point>> & out, double e)
{

    std::vector<std::vector<PointPlus>> paths;
    StraighteningMethods::GraphToPaths(paths,  G);



    //! This one might require modifications:
    double ddd = ComputeOriginalDistance(cloud,paths);

//   std::cout << "Original error: " << ddd << std::endl;
    // std::cout << "e factor: "<< e << std::endl;
    StraighteningMethods::Allocate(paths, G, cloud);
    StraighteningMethods::Optimize(out, paths, ddd * e);

    return ddd;

}


