#ifndef POINTSEGMENTDISTANCE_H
#define POINTSEGMENTDISTANCE_H

#include <mlpack/core.hpp>


class PointSegmentDistance
{
public:
    PointSegmentDistance()
    {

    }
    static double Calculate(const arma::mat & segment, const arma::vec & p)
    {
        double t = arma::dot(p-segment.col(0),segment.col(1)-segment.col(0)) / (arma::dot(segment.col(1)-segment.col(0),segment.col(1)-segment.col(0)));
        double u = std::min(std::max(t,(double) 0), (double) 1);
        return arma::norm(segment.col(0)+u*(segment.col(1)-segment.col(0))-p);
        //  return 0;

    }
    static double Calculate(arma::vec a, arma::vec b, arma::vec & p)
    {
        arma::mat tmp(p.n_rows, 2);
        tmp.col(0) = a;
        tmp.col(1) = b;
        return Calculate(tmp, p);
        //  return 0;

    }
protected:

private:
};

#endif // POINTSEGMENTDISTANCE_H
