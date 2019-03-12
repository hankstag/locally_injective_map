#ifndef IS_SIMPLE_POLYGON_H
#define IS_SIMPLE_POLYGON_H

#include <Eigen/Core>
#include <vector>
#include <unordered_set>
#include <set>

template <typename DerivedV>
bool is_simple_polygon(const Eigen::MatrixBase<DerivedV>& P);

template bool is_simple_polygon<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&);

#endif