#include "segments_intersect.h"
#include <igl/segment_segment_intersect.h>
#include <igl/predicates/predicates.h>

template <typename DerivedV>
short orientation(const Eigen::MatrixBase<DerivedV>& P){
    typedef typename DerivedV::Scalar Scalar;
    if(P.rows()==3){
      Eigen::RowVector2d p0,p1,p2;
      p0 << P(0,0),P(0,1);
      p1 << P(1,0),P(1,1);
      p2 << P(2,0),P(2,1);
      switch(igl::predicates::orient2d(p0,p1,p2)){
        case igl::predicates::Orientation::COLLINEAR: return 0; break;
        case igl::predicates::Orientation::NEGATIVE: return -1; break;
        case igl::predicates::Orientation::POSITIVE: return 1; break;
      }
    }else{
        double sum = 0;
        int n = P.rows();
        for(int i=0;i<P.rows();i++){
            sum += (P((i+1)%n,0)-P(i,0))*(P((i+1)%n,1)+P(i,1));
        }
        short result = 0;
        if(sum < 0)
            return 1;
        else if(sum == 0)
            return 0;
        else
            return -1;
    }
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

// assume m,n,p are colinear, check whether p is on [m,n]
// template <typename Scalar>
// bool on_segment(
//     const Scalar m[2],
//     const Scalar n[2],
//     const Scalar p[2],
//     const Scalar eps
// ){
//     // if eps is 0, meaning the segments are closed 
//     // (exactly equal to end points is considered intersection)

//     // if segment is open, meaning there is a neighborhood (-eps,eps)
//     // s.t. if p is within the neighborhood, it's considered on segment
//     return ((p[0] >= std::min(m[0],n[0])-eps) &&
//             (p[0] <= std::max(m[0],n[0])+eps) &&
//             (p[1] >= std::min(m[1],n[1])-eps) &&
//             (p[1] <= std::max(m[1],n[1])+eps));
// }

// template <typename Scalar>
// bool segment_segment_intersect(
//     const Scalar a[2], 
//     const Scalar b[2], 
//     const Scalar c[2], 
//     const Scalar d[2],
//     const Scalar eps
// ){

//     short t1 = igl::predicates::orient2d(a,b,c);
//     short t2 = igl::predicates::orient2d(b,c,d);
//     short t3 = igl::predicates::orient2d(a,b,d);
//     short t4 = igl::predicates::orient2d(a,c,d);
    
//     // colinear case        
//     if((t1 == 0 && on_segment(a,b,c,eps)) ||
//        (t2 == 0 && on_segment(c,d,b,eps)) ||
//        (t3 == 0 && on_segment(a,b,d,eps)) ||
//        (t4 == 0 && on_segment(c,d,a,eps))) return true;
    
//     // ordinary case
//     return (t1 != t3 && t2 != t4);
// }

// template <typename DerivedV>
// bool segment_segment_intersect(
//     const Eigen::MatrixBase<DerivedV>& a,
//     const Eigen::MatrixBase<DerivedV>& b,
//     const Eigen::MatrixBase<DerivedV>& c,
//     const Eigen::MatrixBase<DerivedV>& d,
//     Eigen::PlainObjectBase<DerivedV>& q
// ){
//     Eigen::RowVector3d x,y;
//     Eigen::RowVector3d dir1, dir2;
//     x<<a[0],a[1],0;
//     dir1<<b[0]-a[0],b[1]-a[1],0;
//     y<<c[0],c[1],0;
//     dir2<<d[0]-c[0],d[1]-c[1],0;
//     double t, u;
//     bool succ = igl::segments_intersect(x,dir1,y,dir2,t,u,1e-16);
//     q.resize(2);
//     q[0] = t*dir1[0]+x[0];
//     q[1] = t*dir1[1]+x[1];
//     return succ;
// }

// template <typename DerivedV>
// bool segment_segment_intersect(
//     const Eigen::MatrixBase<DerivedV>& A,
//     const Eigen::MatrixBase<DerivedV>& B,
//     const Eigen::MatrixBase<DerivedV>& C,
//     const Eigen::MatrixBase<DerivedV>& D, 
//     typename DerivedV::Scalar eps
// ){
//     typedef typename DerivedV::Scalar Scalar;

//     Scalar a[2] = {A(0),A(1)};
//     Scalar b[2] = {B(0),B(1)};
//     Scalar c[2] = {C(0),C(1)};
//     Scalar d[2] = {D(0),D(1)};
    
//     return segment_segment_intersect(a,b,c,d,eps);

// }

// void test_segment_segment_intersect(){
//     Eigen::RowVector2d a(0,0);
//     Eigen::RowVector2d b(1,9);
//     Eigen::RowVector2d c(0,0);
//     Eigen::RowVector2d d(0,1);
//     std::cout<<std::setprecision(20)<<b<<std::endl;
//     std::cout<<std::setprecision(20)<<c<<std::endl;
//     std::cout<<"intersect?: "<<segment_segment_intersect(a,b,c,d,0)<<std::endl;
// }
template bool segment_segment_intersect<Eigen::Matrix<double, 1, -1, 1, 1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> >&, bool);
template short orientation<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&);
//template bool segment_segment_intersect<Eigen::Matrix<double, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&, bool);
//template bool segment_segment_intersect<Eigen::Matrix<double, 1, -1, 1, 1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> >&, bool);
// template short orientation<Eigen::Matrix<double, 3, 2, 0, 3, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, 3, 2, 0, 3, 2> > const&);
// template short orientation<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&);
// template bool segment_segment_intersect<double>(double const*, double const*, double const*, double const*, double);
// template bool segment_segment_intersect<Eigen::Matrix<double, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::Matrix<double, 1, 2, 1, 1, 2>::Scalar);
// template bool segment_segment_intersect<Eigen::Matrix<double, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&, bool);
// template bool segment_segment_intersect<Eigen::Matrix<double, 1, -1, 1, 1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> >&, bool);