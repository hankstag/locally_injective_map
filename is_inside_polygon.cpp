#include "is_inside_polygon.h"

// Ray casting algorithm
// [https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/]
template <typename Scalar>
bool point_in_poly(
    const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& P,
    const Eigen::Matrix<Scalar,1,2>& q
){
    // There must be at least 3 vertices in polygon[]
    if (P.rows() < 3)  return false;
    
    // pick a far right vertex (outside P)
    Scalar r = P.col(0).maxCoeff() + 2.0f;
    Eigen::Matrix<Scalar,1,2> q2,_q;
    q2 << r, q(1);
    int count = 0;
    for(int i=0;i<P.rows();i++){
        Eigen::Matrix<Scalar,1,2> a = P.row(i);
        Eigen::Matrix<Scalar,1,2> b = P.row((i+1)%P.rows());
        if(segment_segment_intersect(a,b,q,q2,_q)){
            count++;
        }
    }
    return count&1;
}

template bool point_in_poly<double>(Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<double, 1, 2, 1, 1, 2> const&);