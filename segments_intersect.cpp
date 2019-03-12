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

template <typename DerivedV>
bool segment_segment_intersect(
    typename DerivedV::Scalar a[2], 
    typename DerivedV::Scalar b[2], 
    typename DerivedV::Scalar c[2], 
    typename DerivedV::Scalar d[2], 
    Eigen::PlainObjectBase<DerivedV>& q, 
    bool calc
){

    typedef typename DerivedV::Scalar Scalar;
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
    if(t1==0 && t2==0 && t3==0 && t4==0){
        if((std::max(a[0],b[0])-std::min(c[0],d[0])<1e-15 || std::min(a[0],b[0])-std::max(c[0],d[0])>-1e-15) &&
           (std::max(a[1],b[1])-std::min(c[1],d[1])<1e-15 || std::min(a[1],b[1])-std::max(c[1],d[1])>-1e-15))
            return false;
        else{
            Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> M(4,2),MS;
            Eigen::VectorXi MI;
            M<<a[0],a[1],b[0],b[1],c[0],c[1],d[0],d[1];
            igl::sortrows(M,false,MS,MI);
            q = MS.row(2);
            return true;
        }
    }

    // meet at ends case
    if(((a[0] == c[0] && a[1] == c[1]) && (b[0] != d[0] || b[1] != d[1])) ||
       ((a[0] != c[0] || a[1] != c[1]) && (b[0] == d[0] && b[1] == d[1])) || 
       ((a[0] == d[0] && a[1] == d[1]) && (b[0] != c[0] || b[1] != c[1])) ||
       ((a[0] != d[0] || a[1] != d[1]) && (b[0] == c[0] && b[1] == c[1])))
        return false;

    if(calc){
        Eigen::Matrix<Scalar,2,1> uv;
        Eigen::Matrix<Scalar,2,2> L;
        L<<b[0]-a[0],c[0]-d[0],
        b[1]-a[1],c[1]-d[1];
        Eigen::Matrix<Scalar,2,1> rhs;
        rhs<<c[0]-a[0],c[1]-a[1];
        uv = L.colPivHouseholderQr().solve(rhs);
        q.resize(2);
        q<<a[0]+(b[0]-a[0])*uv[0],a[1]+(b[1]-a[1])*uv[0];
    }
    return (t1 != t3 && t2 != t4);
}

template <typename DerivedV>
bool segment_segment_intersect(
    const Eigen::MatrixBase<DerivedV>& A,
    const Eigen::MatrixBase<DerivedV>& B,
    const Eigen::MatrixBase<DerivedV>& C,
    const Eigen::MatrixBase<DerivedV>& D, 
    Eigen::PlainObjectBase<DerivedV>& q, 
    bool calc
){
    typedef typename DerivedV::Scalar Scalar;

    Scalar a[2] = {A(0),A(1)};
    Scalar b[2] = {B(0),B(1)};
    Scalar c[2] = {C(0),C(1)};
    Scalar d[2] = {D(0),D(1)};
    
    return segment_segment_intersect(a,b,c,d,q,calc);

}


template short orientation<Eigen::Matrix<double, 3, 2, 0, 3, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, 3, 2, 0, 3, 2> > const&);
template short orientation<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&);
template bool segment_segment_intersect<Eigen::Matrix<double, 1, -1, 1, 1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> >&, bool);
template bool segment_segment_intersect<Eigen::Matrix<double, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&, bool);