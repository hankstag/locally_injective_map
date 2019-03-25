#include "is_simple_polygon.h"
#include "segments_intersect.h"
#include <queue>
#include <algorithm>
#include <fstream>
#include <igl/copyleft/cgal/orient2D.h>
#include <iostream>

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
        // TODO: is there a more concise way
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
                            return ar[1] < br[1];
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

// check whether Segment a and Segment b intersect each other
bool disjoint_segment_intersect(Segment s1, Segment s2, int n){
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
    return segment_segment_intersect(a,b,c,d,0.0);
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

    // build event list (sort segments)
    auto later = [](const Event& a, const Event& b){
        return a > b;
    };
    std::priority_queue<Event,std::vector<Event>,decltype(later)> Q(later);
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
        Segment seg=segments[std::get<1>(evt)];
        bool is_enter = std::get<2>(evt);
        int seg_id = std::get<2>(seg);
        if(is_enter){
            sl.insert(seg);
            auto pos = sl.find(seg);
            auto prev = (pos == sl.begin()) ? pos: std::prev(pos);
            auto next = std::next(pos);
            if((pos  != sl.begin() && disjoint_segment_intersect(*pos,*prev,n)) ||
               (next != sl.end()   && disjoint_segment_intersect(*pos,*next,n)))
                return false;
        }else{
            auto pos = sl.find(seg);
            auto prev = (pos == sl.begin()) ? pos : std::prev(pos);
            auto next = (pos == sl.end()) ? pos : std::next(pos);
            if(pos != sl.begin() && next != sl.end() && disjoint_segment_intersect(*prev,*next,n))
                return false;
            sl.erase(seg);
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
    Eigen::MatrixXd poly(7,2);
    //poly<<0.875,-10.875,1,-11,1.125,-11.125,3,-9; 
    //poly<<0,0,4,0,4,2,2,2,5,2,4,4,0,4;
    poly<<0,0,4,0,4,4,2,4,2,2,2,4,0,4,0,0;
    // load_curve(poly,"complex_polygon4");
    std::cout<<"check simple: "<<is_simple_polygon(poly)<<std::endl;
    exit(0);
}