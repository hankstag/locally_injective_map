#ifndef POINT_IN_POLY_H
#define POINT_IN_POLY_H

#include <Eigen/Core>
#include "segments_intersect.h"

template <typename Scalar>
bool point_in_poly(
    const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& P,
    const Eigen::Matrix<Scalar,1,2>& q
);

#endif