#ifndef IS_SIMPLE_POLYGON_H
#define IS_SIMPLE_POLYGON_H

#include <Eigen/Core>
#include <vector>
#include <unordered_set>
#include <set>

void test_is_simple_polygon();
bool is_simple_polygon(const Eigen::MatrixXd& P);

#endif