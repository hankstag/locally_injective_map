#include "is_simple_polygon.h"
#include "segments_intersect.h"
#include <queue>
#include <algorithm>

typedef std::pair<double,double> Point;
typedef std::pair<Point,Point> Segment;
typedef std::tuple<Point,int,bool> Event; // (point, segment_id, enter/leave)
struct Order{
    bool operator() (const Segment& a, const Segment& b) const {
        return a < b;
    }
};
using SweepList = std::set<Segment,Order>;

void plot_current_sl(
    const SweepList& seglist
){
    std::cout<<"current segment list: =====> "<<std::endl;
    for(auto seg: seglist){
        std::cout<<"("<<seg.first.first<<","<<seg.first.second<<")"<<" - ("<<seg.second.first<<","<<seg.second.second<<")\n";
    }
    std::cout<<"current segment list: <===== "<<std::endl;

}

// check whether Segment a and Segment b intersect each other
bool is_segment_intersect(
    Segment s1,
    Segment s2
){
    double a[2] = {s1.first.first,s1.first.second};
    double b[2] = {s1.second.first,s1.second.second};
    double c[2] = {s2.first.first,s2.first.second};
    double d[2] = {s2.second.first,s2.second.second};
    return segment_segment_intersect(a,b,c,d,1e-16);
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
            segments.push_back(Segment(p2,p1));
        else 
            segments.push_back(Segment(p1,p2));
    }

    // build event list (sort segments)
    auto later = [](const Event& a, const Event& b){
        return a > b;
    };
    // for(Segment seg: segments){
    //     std::cout<<"("<<seg.first.first<<","<<seg.first.second<<")"<<" - ("<<seg.second.first<<","<<seg.second.second<<")\n";
    // }

    std::priority_queue<Event,std::vector<Event>,decltype(later)> Q(later);
    std::cout<<"#segments "<<segments.size()<<std::endl;
    for(int i=0;i<segments.size();i++){
        Q.push(Event(segments[i].first,i,true));
        Q.push(Event(segments[i].second,i,false));
    }
    SweepList sl;
    auto is_intersected = [&sl](
        const SweepList::iterator& a,
        const SweepList::iterator& b
    ){
        Segment s1 = (*a);
        Segment s2 = (*b);

        return false;
    };
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
        if(is_enter){
            sl.insert(seg);
            auto pos = sl.find(seg);
            auto prev = std::prev(pos);
            auto next = std::next(pos);
            if((prev != sl.begin() && is_segment_intersect(*pos,*prev)) ||
               (next != sl.end()   && is_segment_intersect(*pos,*next)))
                return false;
        }else{
            plot_current_sl(sl);
            auto pos = sl.find(seg);
            auto prev = std::prev(pos);
            auto next = std::next(pos);
            if(prev != sl.begin() && next != sl.end() && is_segment_intersect(*prev,*next))
                return false;
            sl.erase(seg);
        }
    }
    return true;
}

void test_is_simple_polygon(){
    Eigen::MatrixXd P(6,2);
    P<<0,0,0,1,0,3,0,5,0,8,1,9;
    is_simple_polygon(P);
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
