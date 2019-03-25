#include "is_simple_polygon.h"
#include "segments_intersect.h"
#include <queue>
#include <algorithm>
#include <fstream>
#include <igl/copyleft/cgal/orient2D.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 cPoint;
typedef CGAL::Polygon_2<K> Polygon_2;

typedef std::pair<double,double> Point;
typedef std::tuple<Point,Point,int> Segment;
typedef std::tuple<Point,int,bool> Event; // (point, segment_id, enter/leave)

struct Order{
    bool operator() (const Segment& a, const Segment& b) const {
        double al[2] = {std::get<0>(a).first, std::get<0>(a).second};
        double ar[2] = {std::get<1>(a).first, std::get<1>(a).second};
        double bl[2] = {std::get<0>(b).first, std::get<0>(b).second};
        double br[2] = {std::get<1>(b).first, std::get<1>(b).second};
        int id1 = std::get<2>(a);
        int id2 = std::get<2>(b);
        if(id1 == id2) return false;
        if (al[0] <= bl[0]) {
            short s = igl::copyleft::cgal::orient2D(al,ar,bl);
            if (s != 0)
                return s > 0;
            else {
                if (al[0] == ar[0])     // special case of vertical line.
                    return al[1]<bl[1];
                else{
                    short t = igl::copyleft::cgal::orient2D(al,ar,br);
                    if(t != 0){
                        return t > 0;
                    }else{
                        if(ar[1] == br[1])
                            return ar[0] < br[0];
                        else 
                            return ar[1]<br[1];
                    }
                }
            }
        } else {
            short s = igl::copyleft::cgal::orient2D(bl, br, al);
            if (s != 0)
                return s < 0;
            else
                return igl::copyleft::cgal::orient2D(bl, br, ar) < 0;
        }
    }
};
using SweepList = std::set<Segment,Order>;

void plot_segment(
    const Segment& seg
){
    std::cout<<std::get<2>(seg)<<"("<<std::get<0>(seg).first<<","<<std::get<0>(seg).second<<")"<<" - ("<<std::get<1>(seg).first<<","<<std::get<1>(seg).second<<")\n";
}

void in_tree(){
    for(int i=0;i<segments_in_tree.rows();i++){
        std::cout<<i<<": "<<segments_in_tree(i)<<std::endl;
    }
}

void plot_current_sl(
    const SweepList& seglist
){
    std::cout<<"current segment list: =====> "<<std::endl;
    int i=0;
    for(auto seg: seglist){
        std::cout<<std::get<2>(seg)<<" ("<<std::get<0>(seg).first<<","<<std::get<0>(seg).second<<")"<<" - ("<<std::get<1>(seg).first<<","<<std::get<1>(seg).second<<")\n";
    }
    std::cout<<"current segment list: <===== "<<std::endl;
    std::cout<<std::endl;

}

// check whether Segment a and Segment b intersect each other
bool disjoint_segment_intersect(
    Segment s1,
    Segment s2,
    int n
){
    // std::cout<<"test test"<<std::endl;
    // plot_segment(s1);
    // plot_segment(s2);
    int id1 = std::get<2>(s1);
    int id2 = std::get<2>(s2);
    if(std::abs(id1-id2)==1 || std::abs(id1-id2)==n-1){
        // if s1 and s2 are exactly the same
        if(std::get<0>(s1) == std::get<0>(s2) &&
           std::get<1>(s1) == std::get<1>(s2))
            return true;
        else
            return false;
    }
    double a[2] = {std::get<0>(s1).first,std::get<0>(s1).second};
    double b[2] = {std::get<1>(s1).first,std::get<1>(s1).second};
    double c[2] = {std::get<0>(s2).first,std::get<0>(s2).second};
    double d[2] = {std::get<1>(s2).first,std::get<1>(s2).second};
    bool x = segment_segment_intersect(a,b,c,d,0.0);
    //std::cout<<"intersect? "<<x<<std::endl;
    return x;
}

bool is_simple_polygon(const Eigen::MatrixXd& P){
    // P: [a, b, c, d,...]
    // change it to the form of edges
    // E: [a, b; b, c; c, d;...]

    // build segments list
    std::vector<Segment> segments;
    for(int i=0;i<P.rows();i++){
        int i_1 = (i+1)%P.rows();
        Point p1 = Point(P(i,  0),P(i  ,1));
        Point p2 = Point(P(i_1,0),P(i_1,1));
        if(p1>p2)
            segments.push_back(Segment(p2,p1,segments.size()));
        else 
            segments.push_back(Segment(p1,p2,segments.size()));
    }
    segments_in_tree.setZero(segments.size());

    // build event list (sort segments)
    auto later = [](const Event& a, const Event& b){
        return a > b;
    };
    // for(Segment seg: segments){
    //     std::cout<<std::get<2>(seg)<<"("<<std::get<0>(seg).first<<","<<std::get<0>(seg).second<<")"<<" - ("<<std::get<1>(seg).first<<","<<std::get<1>(seg).second<<")\n";
    // }
    std::cout<<std::setprecision(17)<<P<<std::endl;
    std::priority_queue<Event,std::vector<Event>,decltype(later)> Q(later);
    std::cout<<"#segments "<<segments.size()<<std::endl;
    for(int i=0;i<segments.size();i++){
        Q.push(Event(std::get<0>(segments[i]),i,true ));
        Q.push(Event(std::get<1>(segments[i]),i,false));
    }
    SweepList sl;
    int n = segments.size();
    // start line sweeping
    while(!Q.empty()){
        Event evt = Q.top();
        Q.pop();
        std::cout<<std::get<1>(evt);
        if(std::get<2>(evt))
            std::cout<<" enter"<<std::endl;
        else 
            std::cout<<" leave\n";
        Segment seg=segments[std::get<1>(evt)];
        bool is_enter = std::get<2>(evt);
        int seg_id = std::get<2>(seg);
        if(is_enter){
            plot_current_sl(sl);
            sl.insert(seg);
            plot_current_sl(sl);
            auto pos = sl.find(seg);
            auto prev = (pos == sl.begin()) ? pos: std::prev(pos);
            auto next = std::next(pos);
            if((pos  != sl.begin() && disjoint_segment_intersect(*pos,*prev,n)) ||
               (next != sl.end()   && disjoint_segment_intersect(*pos,*next,n)))
                return false;
            segments_in_tree(seg_id) = true;
            //in_tree();
        }else{
            segments_in_tree(seg_id) = false;
            auto pos = sl.find(seg);
            auto prev = (pos == sl.begin()) ? pos : std::prev(pos);
            auto next = (pos == sl.end()) ? pos : std::next(pos);
            if(pos != sl.begin() && next != sl.end() && disjoint_segment_intersect(*prev,*next,n))
                return false;
            sl.erase(seg);
            // plot_current_sl(sl);
            // in_tree();
        }
    }
    return true;
}

void load_curve(
    Eigen::MatrixXd& polygon,
    const std::string& fileID
){
    std::ifstream myfile;
    myfile.open(fileID);
    polygon.resize(3000,2);
    double x,y;
    int i = 0;
    while (myfile >> x){
        myfile >> y;
        polygon.row(i++) << x,y;
    }
    polygon.conservativeResize(i,2);
}

void test_is_simple_polygon(){
    Eigen::MatrixXd poly(4,2);
    poly<<0.875,-10.875,1,-11,1.125,-11.125,3,-9; 
    load_curve(poly,"complex_polygon2");
    Polygon_2 poly_cgal;
	for(int i=0;i<poly.rows();i++){
		poly_cgal.push_back(cPoint(poly(i,0),poly(i,1)));
    }
    std::cout<<"check simple: "<<is_simple_polygon(poly)<<std::endl;
    std::cout<<"cgal result: "<<poly_cgal.is_simple()<<std::endl;
    exit(0);
}

// Shamos_Hoey https://cs.stackexchange.com/questions/22443/shamos-hoey-line-segment-intersection-runtime
// template <typename DerivedV>
// bool is_simple_polygon(const Eigen::MatrixBase<DerivedV>& P){
//     typedef typename DerivedV::Scalar Scalar;
//     using namespace std;
//     // P: [a, b, c, d,..., n]
//     // change it to the form of edges
//     // E: [a, b; b, c; c, d;...]
//     Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> E;
//     E.resize(P.rows(),4);
//     for(int i=0;i<P.rows();i++)
//         E.row(i)<<P.row(i),P.row((i+1)%P.rows());
    
//     // sort ends of segments, pointing up-right
//     for(int i=0;i<E.rows();i++){
//         if(E(i,0) > E(i,2)){
//             std::swap(E(i,2),E(i,0));
//             std::swap(E(i,3),E(i,1));
//         }else if(E(i,0) == E(i,2) && E(i,1) > E(i,3)){
//             swap(E(i,3),E(i,1));
//         }
//     }
//     // build the event queue
//     typedef tuple<Scalar,Scalar,Scalar,int> Evt;
//     typedef tuple<Scalar,Scalar,Scalar,Scalar,int> Seg;
//     auto cmp = [](Evt left, Evt right) { // is left > right? (small one happen first)
//         if(get<0>(left) == get<0>(right)){
//             if(get<1>(left) == get<1>(right)){
//                 if(get<2>(left) == get<2>(right)){
//                     return get<3>(left) > get<3>(right);
//                 }else   
//                     return get<2>(left) > get<2>(right); // if endpoints meet, deletion happen first
//             }else
//                 return get<1>(left) > get<1>(right); // small y happen first
//         }else{  
//             return get<0>(left) > get<0>(right);     // small x happen first
//         }
//     };
//     priority_queue<Evt, vector<Evt>, decltype(cmp)> Q(cmp);
//     for(int i=0;i<E.rows();i++){
//         Q.push(Evt(E(i,0),E(i,1),true,i));  // left endpoint
//         Q.push(Evt(E(i,2),E(i,3),false,i)); // right ...
//     }
    
//     // do a pass on Q, to check the case where there are at least 3 duplicate vertices
//     // TODO: this check is problematic
//     #define CORNER_CASE
//     #ifndef CORNER_CASE
//     auto Q_aux = Q;
//     int n_dup = 0;
//     Scalar o1,o2;
//     while(!Q_aux.empty()){
//         Evt p = Q_aux.top();
//         Q_aux.pop();
//         if(o1==std::get<0>(p) && o2==std::get<1>(p)){
//             if(++n_dup>2){
//                 std::cout<<"duplicated "<<n_dup<<std::endl;
//                 std::cout<<o1<<","<<o2<<std::endl;
//                 return false;
//             }
//         }else{
//             n_dup = 0;
//             o1 = std::get<0>(p); o2 = std::get<1>(p);
//         }
//     }
//     #endif

//     struct CMP{
//         bool operator() (const Seg& a, const Seg& b) const {
//             auto less = [](const Seg& l, const Seg& r){
//                 // is l less than r?
//                 // if l or r is vertical
//                 if(get<0>(l)==get<2>(l) || get<0>(r)==get<2>(r))
//                     return min(get<1>(l),get<3>(l)) < min(get<1>(r),get<3>(r));
//                 Eigen::Matrix<Scalar,3,2> ti;
//                 ti<<get<0>(r),get<1>(r),
//                     get<0>(l),get<1>(l),
//                     get<2>(l),get<3>(l);
//                 switch(orientation(ti)){
//                     case 1: return true;
//                     case -1: return false;
//                     case 0:{ 
//                         // move r0 to l0
//                         Eigen::Matrix<Scalar,3,2> to;
//                         to<<get<0>(l),get<1>(l),
//                             get<2>(l),get<3>(l),
//                             get<0>(l)+get<2>(r)-get<0>(r),get<1>(l)+get<3>(r)-get<1>(r);
//                         short cc = orientation(to);
//                         if(cc==1) return true;
//                         else if(cc==-1) return false;
//                         else return get<4>(l)<get<4>(r);
//                     }
//                     default: return false;   
//                 };
//             };

//             if(get<0>(a)==get<0>(b) && get<1>(a)==get<1>(b) && 
//                get<2>(a)==get<2>(b) && get<3>(a)==get<3>(b))
//                 return get<4>(a)<get<4>(b);
//             if(get<0>(b) == get<0>(a) && get<1>(b) == get<1>(a)){
//                 Eigen::Matrix<Scalar,3,2> ta;
//                 ta<<get<0>(a),get<1>(a),
//                     get<2>(a),get<3>(a),
//                     get<2>(b),get<3>(b);
//                 switch(orientation(ta)){
//                     case  1:return true; 
//                     case -1:return false;
//                     case  0:return get<4>(a) < get<4>(b);
//                     default: return false;
//                 }
//             }else if(get<0>(b) < get<0>(a))
//                 return !less(b,a);
//             else
//                 return less(a,b);
//         }
//     };

//     unordered_set<int> sl; // segment list
//     // start line sweeping
//     set<Seg,CMP> SL;
//     while(!Q.empty()){
//         Evt p = Q.top();
//         Q.pop();
//         int ei = get<3>(p);
//         Seg t(E(ei,0),E(ei,1),E(ei,2),E(ei,3),ei);
//         if(get<2>(p)){  // left endpoint
//             SL.insert(t);
//             auto pos = SL.find(t);
//             Scalar a[2] = {E(ei,0),E(ei,1)};
//             Scalar b[2] = {E(ei,2),E(ei,3)};
//             std::vector<typename set<Seg,CMP>::iterator> m;
//             if(pos != SL.begin()) m.push_back(prev(pos));
//             if(pos != --SL.end()) m.push_back(next(pos));
//             for(auto nb: m){
//                 Scalar c[2] = {get<0>(*nb),get<1>(*nb)};
//                 Scalar d[2] = {get<2>(*nb),get<3>(*nb)};
//                 int ej = get<4>(*nb);
//                 Eigen::Matrix<Scalar, 1, Eigen::Dynamic> _q;
//                 if(segment_segment_intersect(a,b,c,d,_q,false)){
//                     return false;
//                 }
//             }
//         }else{ // right endpoint
//             auto pos = SL.find(t);
//             if(pos != --SL.end() && pos != SL.begin()){
//                 auto l = prev(pos);
//                 auto r = next(pos);
//                 // check whether l and r intersect
//                 Scalar a[2] = {get<0>(*l),get<1>(*l)};
//                 Scalar b[2] = {get<2>(*l),get<3>(*l)};
//                 Scalar c[2] = {get<0>(*r),get<1>(*r)};
//                 Scalar d[2] = {get<2>(*r),get<3>(*r)};
//                 int e1 = get<4>(*l), e2 = get<4>(*r);
//                 Eigen::Matrix<Scalar, 1, Eigen::Dynamic> _q;
//                 if(segment_segment_intersect(a,b,c,d,_q,false)){
//                     return false;
//                 }
//             }
//             SL.erase(t);
//         }
//     }
//     return true;
// }
