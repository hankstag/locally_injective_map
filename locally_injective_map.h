#ifndef LOCALLY_INJECTIVE_MAP
#define LOCALLY_INJECTIVE_MAP

#include "shor.h"
#include "is_simple_polygon.h"
#include <unordered_map>

#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/AABB.h>
#include <igl/harmonic.h>
#include <igl/in_element.h>
#include <igl/triangle_triangle_adjacency.h>

namespace std{
	template <>
	struct hash<std::pair<int,int>>{
		std::size_t operator()(const std::pair<int,int>& k) const{
			using std::size_t;
			using std::hash;
		    return ((hash<int>()(k.first) ^ (hash<int>()(k.second) << 1)));
		}
	};
}

template <typename Scalar>
class PointMap{
public:
	typedef std::pair<int,int> E;
	typedef std::tuple<Scalar,Scalar,Scalar> P3;
	int new_v; // number of new points
	// vertex id of vertex on source mesh (new inserted) on edge [Pair]
	std::unordered_map<E, std::vector<int>> vid_on_edge;
	// position of new vertex in target mesh
	std::unordered_map<E, std::deque<P3>> vpos_t;
	// ... 					  in source mesh
	std::unordered_map<E, std::deque<P3>> vpos_s;
	// edge to face map
	std::unordered_map<E, int > edgeToFace;
	PointMap(){new_v=0;}
	~PointMap(){}
	void add_point(const E& e,const P3& v_t,const P3& v_s){
		E r_e(e.second,e.first);
		vpos_t[e].push_back(v_t);
		vpos_t[r_e].push_front(v_t);
		vpos_s[e].push_back(v_s);
		vpos_s[r_e].push_front(v_s);
	}
	
	void drop_point(const E& e, int k){
		E r_e(e.second,e.first);
		int nk = vpos_s[e].size()-1-k;
		vpos_s[e].erase(vpos_s[e].begin()+k);
		vpos_s[r_e].erase(vpos_s[r_e].begin()+nk);
		vpos_t[e].erase(vpos_t[e].begin()+k);
		vpos_t[r_e].erase(vpos_t[r_e].begin()+nk);
	}

	void set_index(E p,const std::vector<int>& ids){
		vid_on_edge[p] = ids;
		E rev_key(p.second,p.first);
		std::vector<int> rev_ids(ids.rbegin(),ids.rend());
		vid_on_edge[rev_key] = rev_ids;
	}

};

template <typename DerivedV>
void generate_map(
    const Eigen::MatrixXi& Fs,
    const Eigen::MatrixXi& Ft,  
    const Eigen::PlainObjectBase<DerivedV>& harmo_s,
    const Eigen::PlainObjectBase<DerivedV>& harmo_t,
    Eigen::PlainObjectBase<DerivedV>& V,
    Eigen::MatrixBase<DerivedV>& Vd,
    Eigen::MatrixXi& Fn,
    Eigen::PlainObjectBase<DerivedV>& uv);
    
template <typename DerivedV>
void get_image(
	const Eigen::PlainObjectBase<DerivedV>& V,
	const Eigen::MatrixXi& F,
	const Eigen::PlainObjectBase<DerivedV>& Q,
	const igl::AABB<Eigen::Matrix<typename DerivedV::Scalar,Eigen::Dynamic,Eigen::Dynamic>,2>& tree,
	const Eigen::MatrixBase<DerivedV>& Vd,
	Eigen::PlainObjectBase<DerivedV>& img_Q
);

template <typename DerivedV>
void create_poly_from_edges(
    const std::vector<std::pair<typename DerivedV::Scalar,typename DerivedV::Scalar>>& uv_vec,
    PointMap<typename DerivedV::Scalar>& pMap,
    std::vector<std::pair<int,int>> edges,
    std::vector<int>& indices,
    Eigen::PlainObjectBase<DerivedV>& P
);

template <typename DerivedV>
int refine(
    std::vector<int>& neg,
    const igl::AABB<Eigen::Matrix<typename DerivedV::Scalar,Eigen::Dynamic,Eigen::Dynamic>,2>& s_tree,
    const igl::AABB<Eigen::Matrix<typename DerivedV::Scalar,Eigen::Dynamic,Eigen::Dynamic>,2>& t_tree,
    const Eigen::MatrixXi& Fs,
    const Eigen::MatrixXi& Ft,
    const Eigen::PlainObjectBase<DerivedV>& vh_s,
    const Eigen::PlainObjectBase<DerivedV>& vh_t,
    Eigen::PlainObjectBase<DerivedV>& V,
    Eigen::MatrixBase<DerivedV>& Vd,
    Eigen::MatrixXi& Fn,
    Eigen::PlainObjectBase<DerivedV>& uv
);

template <typename DerivedV>
void simplify(PointMap<typename DerivedV::Scalar>& pMap, const Eigen::MatrixXi& Fn, const std::vector<std::pair<typename DerivedV::Scalar,typename DerivedV::Scalar>>& uv_vec);

template <typename DerivedV>
void refineTriangle(
    int triIndex,
    const Eigen::MatrixBase<DerivedV>& V,
    const Eigen::MatrixBase<DerivedV>& Vd,
    const Eigen::PlainObjectBase<DerivedV>& vh_s,
    const Eigen::PlainObjectBase<DerivedV>& vh_t,
    const Eigen::MatrixXi& FF_s,
    const Eigen::MatrixXi& FF_t,
    const Eigen::MatrixXi& Fs,
    const Eigen::MatrixXi& Ft,
    PointMap<typename DerivedV::Scalar>& pMap, 
    std::vector<std::pair<typename DerivedV::Scalar,typename DerivedV::Scalar>>& uv_vec,
    const igl::AABB<Eigen::Matrix<typename DerivedV::Scalar,Eigen::Dynamic,Eigen::Dynamic>,2>& s_tree,
    const igl::AABB<Eigen::Matrix<typename DerivedV::Scalar,Eigen::Dynamic,Eigen::Dynamic>,2>& t_tree
);

template <typename DerivedV>
bool triangulate_face(
    int triIndex, 
    Eigen::MatrixXi& Fn, 
    PointMap<typename DerivedV::Scalar>& pMap, 
    const std::vector<std::pair<typename DerivedV::Scalar,typename DerivedV::Scalar>>& uv_vec
);

void locally_injective_map(
    Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& P,
    const Eigen::VectorXi& R,
    const Eigen::VectorXi& T,     // T(i) is index in V
    Eigen::MatrixXd& UV,
    Eigen::MatrixXi& F_o
);

#endif
