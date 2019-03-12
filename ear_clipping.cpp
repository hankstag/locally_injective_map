#include "ear_clipping.h"
#include <igl/copyleft/cgal/orient2D.h>

template <typename DerivedP>
bool is_ear(
	const Eigen::MatrixBase<DerivedP>& P,
	const Eigen::VectorXi& RT,
	const Eigen::VectorXi& L,
	const Eigen::VectorXi& R, 
	const int i
){
	typedef typename DerivedP::Scalar Scalar;    
    Eigen::Matrix<Scalar,3,2> T;
    
	// valid neighbor
    int a = L(i), b = R(i);
	if(RT(i) != 0 || RT(a) != 0 || RT(b) != 0) return false;
    double A[2] = {double(P(a, 0)), double(P(a, 1))};
    double B[2] = {double(P(i, 0)), double(P(i, 1))};
    double C[2] = {double(P(b, 0)), double(P(b, 1))};
    if(igl::copyleft::cgal::orient2D(A, B, C)<=0) return false;

    std::vector<int> l={a,i,b};
    int k=b;
    while(R(k)!=b){
        // test whether valid vertex is inside
		if(k == a || k == b){
			k = R(k);
			continue;
		}
        bool inside=true;
        for(int f=0;f<3;f++){
            Eigen::Matrix<Scalar,3,2> S;
            double A[2] = {P(l[f],0),P(l[f],1)};
            double B[2] = {P(l[(f+1)%3],0),P(l[(f+1)%3],1)};
            double C[2] = {P(k,0),P(k,1)};
            if(igl::copyleft::cgal::orient2D(A, B, C)<0){
                inside = false;
                break;
            }
        }
        if(inside)
            return false;
        k=R(k);
    }
    return true;
}

template <typename DerivedP>
void ear_clipping(
	const Eigen::MatrixBase<DerivedP>& P,
	const Eigen::VectorXi& RT,
	Eigen::VectorXi& C,
	Eigen::MatrixXi& eF, // #ear*3 indices of ears
	Eigen::PlainObjectBase<DerivedP>& nP,
	Eigen::VectorXi& nR
){
	Eigen::VectorXi X; // remaining vertices
	Eigen::VectorXi L(P.rows());
	Eigen::VectorXi R(P.rows());
	for(int i=0;i<P.rows();i++){
		L(i) = (i-1+P.rows())%P.rows();
		R(i) = (i+1)%P.rows();
	}
	// initialize ears
	Eigen::VectorXi I(P.rows());
	I.setZero();
	X = I;
	for(int i=0;i<P.rows();i++){
		I(i) = is_ear(P,RT,L,R,i);
	}
	while(I.maxCoeff()==1){  // still have ear
		// find the first ear
		int e = 0;
		for(;e<I.rows() && I(e)!=1;e++);
		assert(I(e)!=0);
		// find valid neighbor
		int a = L(e), b = R(e);
		if(a == b) break;

		// ear cut
		eF.conservativeResize(eF.rows()+1,3);
		eF.bottomRows(1)<<a,e,b;
		L(b) = a;
		L(e) = -1;
		R(a) = b;
		R(e) = -1;

		I(e) = 0; // mark it as non-ear
		// update neighbor's ear status
		I(a) = is_ear(P,RT,L,R,a);
		I(b) = is_ear(P,RT,L,R,b);
		X(e) = 1;

		// corner case where only one edge left
		if(L(a)==b && R(b)==a){
			X(a) = 1;
			X(b) = 1;
		}
	}
	for(int i=0;i<X.rows();i++)
		X(i) = 1-X(i);
	C.resize(X.sum());
	int j=0;
	for(int i=0;i<X.rows();i++)
		if(X(i)==1){
			C(j++) = i;
		}
	igl::slice(P,C,1,nP);
	igl::slice(RT,C,nR);
}

//template void ear_clipping<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, -1, 0, -1, -1>&);
template void ear_clipping<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, -1, 0, -1, -1>&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&);
//template void ear_clipping<Eigen::Matrix<mpfr::mpreal, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<mpfr::mpreal, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, -1, 0, -1, -1>&, Eigen::PlainObjectBase<Eigen::Matrix<mpfr::mpreal, -1, -1, 0, -1, -1> >&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&);