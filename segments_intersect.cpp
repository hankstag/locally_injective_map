#include "segments_intersect.h"
#include <iostream>
#include <queue>
#include <algorithm>
#include <igl/copyleft/cgal/orient2D.h>
#include <bitset> 

template <typename DerivedV>
short orientation(const Eigen::MatrixBase<DerivedV>& P){
    typedef typename DerivedV::Scalar Scalar;
    double a[2] = {double(P(0, 0)), double(P(0, 1))};
    double b[2] = {double(P(1, 0)), double(P(1, 1))};
    double c[2] = {double(P(2, 0)), double(P(2, 1))};
    return igl::copyleft::cgal::orient2D(a, b, c);
}

// assume m,n,p are colinear, check whether p is on [m,n]
template <typename Scalar>
bool on_segment(
    const Scalar m[2],
    const Scalar n[2],
    const Scalar p[2],
    const Scalar eps
){
    // if eps is 0, meaning the segments are closed 
    // (exactly equal to end points is considered interection)

    // if segment is open, meaning there is a neighborhood (-eps,eps)
    // s.t. if p landing inside the neighborhood, it's considered not on segment

    return !((p[0]-std::min(m[0],n[0])<eps) || (p[0]-std::max(m[0],n[0])>-eps));
}

template <typename Scalar>
bool segment_segment_intersect(
    const Scalar a[2], 
    const Scalar b[2], 
    const Scalar c[2], 
    const Scalar d[2],
    const Scalar eps
){

    Eigen::Matrix<Scalar, 3, 2> T1,T2,T3,T4;
    T1 << a[0],a[1],b[0],b[1],c[0],c[1];
    T2 << b[0],b[1],c[0],c[1],d[0],d[1];
    T3 << a[0],a[1],b[0],b[1],d[0],d[1];
    T4 << a[0],a[1],c[0],c[1],d[0],d[1];
    auto t1 = orientation(T1);
    auto t2 = orientation(T2);
    auto t3 = orientation(T3);
    auto t4 = orientation(T4);
    
    // colinear case        
    if((t1 == 0 && on_segment(a,b,c,eps)) ||
       (t2 == 0 && on_segment(c,d,b,eps)) ||
       (t3 == 0 && on_segment(a,b,d,eps)) ||
       (t4 == 0 && on_segment(c,d,a,eps))) return true;
    
    // ordinary case
    return (t1 != t3 && t2 != t4);
}

template <typename DerivedV>
bool segment_segment_intersect(
    const Eigen::MatrixBase<DerivedV>& A,
    const Eigen::MatrixBase<DerivedV>& B,
    const Eigen::MatrixBase<DerivedV>& C,
    const Eigen::MatrixBase<DerivedV>& D, 
    typename DerivedV::Scalar eps
){
    typedef typename DerivedV::Scalar Scalar;

    Scalar a[2] = {A(0),A(1)};
    Scalar b[2] = {B(0),B(1)};
    Scalar c[2] = {C(0),C(1)};
    Scalar d[2] = {D(0),D(1)};
    
    return segment_segment_intersect(a,b,c,d,eps);

}

void test_segment_segment_intersect(){
    Eigen::RowVector2d a(0,0);
    Eigen::RowVector2d b(0,1.1);
    Eigen::RowVector2d c(0,std::nextafter(1.1,2.0));
    Eigen::RowVector2d d(0,2.2);
    std::cout<<std::setprecision(20)<<b<<std::endl;
    std::cout<<std::setprecision(20)<<c<<std::endl;
    std::cout<<"intersect?: "<<segment_segment_intersect(a,b,c,d,0)<<std::endl;
}

template short orientation<Eigen::Matrix<double, 3, 2, 0, 3, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, 3, 2, 0, 3, 2> > const&);
template short orientation<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&);
//template bool segment_segment_intersect<Eigen::Matrix<double, 1, -1, 1, 1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> >&, bool);
//template bool segment_segment_intersect<Eigen::Matrix<double, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&, bool);