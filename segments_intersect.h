#ifndef SEGMENTS_INTERSECT_H
#define SEGMENTS_INTERSECT_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <igl/sortrows.h>
#include <vector>
#include <unordered_set>
#include <set>

// bool segment_segment_intersect(double a[2], double b[2], double c[2], double d[2], Eigen::RowVectorXd& q, bool calc=false);
// bool segment_segment_intersect(Eigen::RowVectorXd& A,Eigen::RowVectorXd& B,Eigen::RowVectorXd& C,Eigen::RowVectorXd& D,Eigen::RowVectorXd& P,bool calc=false);
// bool is_simple_polygon(const Eigen::MatrixXd& P);
// bool orientation(Eigen::MatrixXd& P);

template <typename DerivedV>
bool is_simple_polygon(const Eigen::MatrixBase<DerivedV>& P);

template <typename DerivedV>
bool segment_segment_intersect(
    const Eigen::MatrixBase<DerivedV>& A,
    const Eigen::MatrixBase<DerivedV>& B,
    const Eigen::MatrixBase<DerivedV>& C,
    const Eigen::MatrixBase<DerivedV>& D, 
    Eigen::PlainObjectBase<DerivedV>& q, 
    bool calc=false
);

template <typename DerivedV>
bool segment_segment_intersect(
    typename DerivedV::Scalar a[2], 
    typename DerivedV::Scalar b[2], 
    typename DerivedV::Scalar c[2], 
    typename DerivedV::Scalar d[2], 
    Eigen::PlainObjectBase<DerivedV>& q, 
    bool calc
);

template <typename DerivedV>
short orientation(const Eigen::MatrixBase<DerivedV>& P);

#endif