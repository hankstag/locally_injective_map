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

// Shamos_Hoey https://cs.stackexchange.com/questions/22443/shamos-hoey-line-segment-intersection-runtime
template <typename DerivedV>
bool is_simple_polygon(const Eigen::MatrixBase<DerivedV>& P){
    typedef typename DerivedV::Scalar Scalar;
    using namespace std;
    // P: [a, b, c, d,..., n]
    // change it to the form of edges
    // E: [a, b; b, c; c, d;...]
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> E;
    E.resize(P.rows(),4);
    for(int i=0;i<P.rows();i++)
        E.row(i)<<P.row(i),P.row((i+1)%P.rows());
    
    // sort ends of segments, pointing up-right
    for(int i=0;i<E.rows();i++){
        if(E(i,0) > E(i,2)){
            std::swap(E(i,2),E(i,0));
            std::swap(E(i,3),E(i,1));
        }else if(E(i,0) == E(i,2) && E(i,1) > E(i,3)){
            swap(E(i,3),E(i,1));
        }
    }
    // build the event queue
    typedef tuple<Scalar,Scalar,Scalar,int> Evt;
    typedef tuple<Scalar,Scalar,Scalar,Scalar,int> Seg;
    auto cmp = [](Evt left, Evt right) { // is left > right? (small one happen first)
        if(get<0>(left) == get<0>(right)){
            if(get<1>(left) == get<1>(right)){
                if(get<2>(left) == get<2>(right)){
                    return get<3>(left) > get<3>(right);
                }else   
                    return get<2>(left) > get<2>(right); // if endpoints meet, deletion happen first
            }else
                return get<1>(left) > get<1>(right); // small y happen first
        }else{  
            return get<0>(left) > get<0>(right);     // small x happen first
        }
    };
    priority_queue<Evt, vector<Evt>, decltype(cmp)> Q(cmp);
    for(int i=0;i<E.rows();i++){
        Q.push(Evt(E(i,0),E(i,1),true,i));  // left endpoint
        Q.push(Evt(E(i,2),E(i,3),false,i)); // right ...
    }
    
    // do a pass on Q, to check the case where there are at least 3 duplicate vertices
    // TODO: this check is problematic
    #define CORNER_CASE
    #ifndef CORNER_CASE
    auto Q_aux = Q;
    int n_dup = 0;
    Scalar o1,o2;
    while(!Q_aux.empty()){
        Evt p = Q_aux.top();
        Q_aux.pop();
        if(o1==std::get<0>(p) && o2==std::get<1>(p)){
            if(++n_dup>2){
                std::cout<<"duplicated "<<n_dup<<std::endl;
                std::cout<<o1<<","<<o2<<std::endl;
                return false;
            }
        }else{
            n_dup = 0;
            o1 = std::get<0>(p); o2 = std::get<1>(p);
        }
    }
    #endif

    struct CMP{
        bool operator() (const Seg& a, const Seg& b) const {
            auto less = [](const Seg& l, const Seg& r){
                // is l less than r?
                // if l or r is vertical
                if(get<0>(l)==get<2>(l) || get<0>(r)==get<2>(r))
                    return min(get<1>(l),get<3>(l)) < min(get<1>(r),get<3>(r));
                Eigen::Matrix<Scalar,3,2> ti;
                ti<<get<0>(r),get<1>(r),
                    get<0>(l),get<1>(l),
                    get<2>(l),get<3>(l);
                switch(orientation(ti)){
                    case 1: return true;
                    case -1: return false;
                    case 0:{ 
                        // move r0 to l0
                        Eigen::Matrix<Scalar,3,2> to;
                        to<<get<0>(l),get<1>(l),
                            get<2>(l),get<3>(l),
                            get<0>(l)+get<2>(r)-get<0>(r),get<1>(l)+get<3>(r)-get<1>(r);
                        short cc = orientation(to);
                        if(cc==1) return true;
                        else if(cc==-1) return false;
                        else return get<4>(l)<get<4>(r);
                    }
                    default: return false;   
                };
            };

            if(get<0>(a)==get<0>(b) && get<1>(a)==get<1>(b) && 
               get<2>(a)==get<2>(b) && get<3>(a)==get<3>(b))
                return get<4>(a)<get<4>(b);
            if(get<0>(b) == get<0>(a) && get<1>(b) == get<1>(a)){
                Eigen::Matrix<Scalar,3,2> ta;
                ta<<get<0>(a),get<1>(a),
                    get<2>(a),get<3>(a),
                    get<2>(b),get<3>(b);
                switch(orientation(ta)){
                    case  1:return true; 
                    case -1:return false;
                    case  0:return get<4>(a) < get<4>(b);
                    default: return false;
                }
            }else if(get<0>(b) < get<0>(a))
                return !less(b,a);
            else
                return less(a,b);
        }
    };

    unordered_set<int> sl; // segment list
    // start line sweeping
    set<Seg,CMP> SL;
    while(!Q.empty()){
        Evt p = Q.top();
        Q.pop();
        int ei = get<3>(p);
        Seg t(E(ei,0),E(ei,1),E(ei,2),E(ei,3),ei);
        if(get<2>(p)){  // left endpoint
            SL.insert(t);
            auto pos = SL.find(t);
            Scalar a[2] = {E(ei,0),E(ei,1)};
            Scalar b[2] = {E(ei,2),E(ei,3)};
            std::vector<typename set<Seg,CMP>::iterator> m;
            if(pos != SL.begin()) m.push_back(prev(pos));
            if(pos != --SL.end()) m.push_back(next(pos));
            for(auto nb: m){
                Scalar c[2] = {get<0>(*nb),get<1>(*nb)};
                Scalar d[2] = {get<2>(*nb),get<3>(*nb)};
                int ej = get<4>(*nb);
                Eigen::Matrix<Scalar, 1, Eigen::Dynamic> _q;
                if(segment_segment_intersect(a,b,c,d,_q,false)){
                    return false;
                }
            }
        }else{ // right endpoint
            auto pos = SL.find(t);
            if(pos != --SL.end() && pos != SL.begin()){
                auto l = prev(pos);
                auto r = next(pos);
                // check whether l and r intersect
                Scalar a[2] = {get<0>(*l),get<1>(*l)};
                Scalar b[2] = {get<2>(*l),get<3>(*l)};
                Scalar c[2] = {get<0>(*r),get<1>(*r)};
                Scalar d[2] = {get<2>(*r),get<3>(*r)};
                int e1 = get<4>(*l), e2 = get<4>(*r);
                Eigen::Matrix<Scalar, 1, Eigen::Dynamic> _q;
                if(segment_segment_intersect(a,b,c,d,_q,false)){
                    return false;
                }
            }
            SL.erase(t);
        }
    }
    return true;
}

template short orientation<Eigen::Matrix<double, 3, 2, 0, 3, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, 3, 2, 0, 3, 2> > const&);
template short orientation<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&);
template bool segment_segment_intersect<Eigen::Matrix<double, 1, -1, 1, 1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> >&, bool);
template bool is_simple_polygon<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&);
template bool segment_segment_intersect<Eigen::Matrix<double, 1, 2, 1, 1, 2> >(Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&, bool);