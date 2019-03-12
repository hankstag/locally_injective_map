#ifndef EAR_CLIPPING
#define EAR_CLIPPING

#include <Eigen/Core>
#include <igl/slice.h>
#include "segments_intersect.h"

template <typename DerivedP>
void ear_clipping(
	const Eigen::MatrixBase<DerivedP>& P,
	const Eigen::VectorXi& RT,	
	Eigen::VectorXi& C,
	Eigen::MatrixXi& eF,
	Eigen::PlainObjectBase<DerivedP>& nP,
	Eigen::VectorXi& nR
);

template <typename DerivedP>
bool is_ear(
	const Eigen::MatrixBase<DerivedP>& P,
	const Eigen::VectorXi& RT,
	const Eigen::VectorXi& L,
	const Eigen::VectorXi& R, 
	const int i
);

#endif