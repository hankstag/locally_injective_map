
#include "locally_injective_map.h"
#include <igl/opengl/glfw/Viewer.h>

using namespace std;

template <typename DerivedV>
void generate_map(
    const Eigen::MatrixXi& Fs,
    const Eigen::MatrixXi& Ft,
    const Eigen::PlainObjectBase<DerivedV>& harmo_s,
    const Eigen::PlainObjectBase<DerivedV>& harmo_t,
    Eigen::PlainObjectBase<DerivedV>& V,
    Eigen::MatrixBase<DerivedV>& Vd, // target triangle mesh
    Eigen::MatrixXi& Fn,
    Eigen::PlainObjectBase<DerivedV>& uv)
{
    typedef typename DerivedV::Scalar Scalar;
    Fn = Fs;
    // calculate reverse map
	igl::AABB<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>,2> aabb_s;
    igl::AABB<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>,2> aabb_t;
  	aabb_s.init(harmo_s,Fs);
	aabb_t.init(harmo_t,Ft);

    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> vd = Vd.leftCols(2);
	get_image(harmo_t,Ft,harmo_s,aabb_t,vd,uv);

    // collect flipped
    std::vector<int> flip_list;
    for(int i=0;i<Fs.rows();i++){
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> tri(3,2);
        tri<<uv.row(Fs(i,0)),uv.row(Fs(i,1)),uv.row(Fs(i,2));
        if(orientation(tri)<=0)
            flip_list.emplace_back(i);
    }
    // refinement
    refine(flip_list, aabb_s, aabb_t, Fs,Ft, harmo_s, harmo_t, V, Vd, Fn, uv);

}

template <typename DerivedV>
void get_image(
	const Eigen::PlainObjectBase<DerivedV>& V,
	const Eigen::MatrixXi& F,
	const Eigen::PlainObjectBase<DerivedV>& Q,
	const igl::AABB<Eigen::Matrix<typename DerivedV::Scalar,Eigen::Dynamic,Eigen::Dynamic>,2>& tree,
	const Eigen::MatrixBase<DerivedV>& Vd,
	Eigen::PlainObjectBase<DerivedV>& img_Q
){
    typedef typename DerivedV::Scalar Scalar;
	img_Q.resize(Q.rows(),Vd.cols());
	Eigen::VectorXi I;
	igl::in_element(V,F,Q,tree,I);
	for(int i=0;i<img_Q.rows();i++){
		int a = F(I(i),0);
		int b = F(I(i),1);
		int c = F(I(i),2);
		Eigen::Matrix<Scalar,3,1> x;
		Eigen::Matrix<Scalar,3,1> rhs;
        rhs<<Q.row(i).transpose(),1.0f;
        Eigen::Matrix<Scalar,3,3> A;
        A<<V(a,0),V(b,0),V(c,0),
           V(a,1),V(b,1),V(c,1),
           1,1,1;
        x = A.colPivHouseholderQr().solve(rhs);
        img_Q.row(i) << Vd.row(a)*x(0)+Vd.row(b)*x(1)+Vd.row(c)*x(2);
	}
}

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
){
    typedef typename DerivedV::Scalar Scalar; 
    typedef std::pair<Scalar,Scalar> P2;   
    int size = neg.size();
    if (size == 0){
        std::cout << "No need to refine.\n";
        return 0;
    }else{
        std::cout << "# of negative triangles: " << size << "\nrefine...\n";
    }
    
    int startSize = V.rows();
    PointMap<Scalar> pMap;
    Eigen::MatrixXi FF_s,FF_t,FFI;
    igl::triangle_triangle_adjacency(Fs,FF_s,FFI);
    igl::triangle_triangle_adjacency(Ft,FF_t,FFI);
    std::vector<int> stack=neg; // store the flipped triangles
    std::vector<bool> isRefined(Fs.rows(),false);
    std::vector<P2> uv_vec(uv.rows());
    for(int i=0;i<uv.rows();i++)
        uv_vec[i] = P2(uv(i,0),uv(i,1));
  
    while(stack.size()>0){
        int f = stack.back();
        stack.pop_back();
        refineTriangle(f, V, Vd, vh_s, vh_t, FF_s, FF_t, Fs, Ft, pMap, uv_vec, s_tree, t_tree);
        isRefined[f] = true;
        for(int i=0;i<3;i++){
            if(FF_s(f,i)<0) continue;
            std::pair<int,int> e1(Fn(FF_s(f,i),0), Fn(FF_s(f,i),1));
            std::pair<int,int> e2(Fn(FF_s(f,i),1), Fn(FF_s(f,i),2));
            std::pair<int,int> e3(Fn(FF_s(f,i),2), Fn(FF_s(f,i),0));
            Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> pg;
            std::vector<int> _list;
            create_poly_from_edges(uv_vec,pMap,{e1,e2,e3},_list,pg);
            if(!(is_simple_polygon(pg) && orientation(pg)>0)){
                if(isRefined[FF_s(f,i)]){
                    //std::cout.precision(30);
                    //std::cout<<pg<<std::endl;
                    //std::cout<<is_simple_polygon(pg)<<","<<dp.is_simple()<<std::endl;
                    //std::cout<<"refined polygon "<<FF_s(f,i)<<" still not simple"<<std::endl;
                    continue;
                    // exit(0);
                }else
                    stack.push_back(FF_s(f,i));
            }
            neg.push_back(FF_s(f,i));
        }
    }
    std::cout << "Done!\nSimplify - trying to reduce number of new points...\n";

#define SIMP 
#ifdef SIMP
    // need to add here simplify function
    simplify(pMap, Fn, uv_vec);
#endif
    int v_num = V.rows();
    uv.conservativeResize(v_num+pMap.new_v,2);
    V.conservativeResize(v_num+pMap.new_v,Eigen::NoChange);
    std::set<std::pair<int,int>> visitedPairs;
    for(auto p: pMap.edgeToFace){
        std::pair<int,int> e = p.first;
        if (visitedPairs.count(e) != 0){
            continue;
        }
        visitedPairs.insert(e);
        visitedPairs.insert(std::make_pair(e.second,e.first));
        auto uv_e = pMap.vpos_t[e];
        auto v_e = pMap.vpos_s[e];

        std::vector<int> v_ids;
        for (int i = 1; i < uv_e.size() - 1; ++i){
            uv.row(v_num)<<std::get<0>(uv_e[i]),std::get<1>(uv_e[i]);
            V.row(v_num)<<std::get<0>(v_e[i]),std::get<1>(v_e[i]),std::get<2>(v_e[i]);
            v_ids.push_back(v_num);
            v_num++;
        }
        // reset all the indices of new points
        pMap.set_index(e,v_ids);
    }
    
    std::vector<bool> isTriangulated(Fs.rows(), false);
    uv_vec.resize(uv.rows());
    for(int i=0;i<uv.rows();i++)
        uv_vec[i] = P2(uv(i,0),uv(i,1));
    for (auto f: neg){
        if (isTriangulated[f])  continue;
        bool result = triangulate_face(f, Fn, pMap, uv_vec);
        assert(result && "triangulation failed");
        isTriangulated[f] = true;
    }
    // drop the (-1,-1,-1) rows
    Eigen::MatrixXi tF=Fn;
	int k=0;
	for(int i=0;i<Fn.rows();i++){
		if(Fn.row(i).sum()!=-3)
			tF.row(k++)<<Fn.row(i);
    }
	tF.conservativeResize(k,3);
	Fn = tF;
    int numOfNewPoints = V.rows() - startSize;
    std::cout << "Done!\n# of new points: " << numOfNewPoints << "\n";
    return 0;
}

template <typename Scalar>
void simplify(PointMap<Scalar>& pMap, const Eigen::MatrixXi& Fn, const std::vector<std::pair<Scalar,Scalar>>& uv_vec)
{
    bool stopFlag = true;
    std::unordered_map<std::pair<int,int>,int> edges_to_check = pMap.edgeToFace;
    std::unordered_map<std::pair<int,int>,int> edgesTmp; // stores the next collection of edges to check

    // whether polygon P is still simple after drop i
    auto check_simple = [](
        const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& P, 
        int i, 
        Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& Q
    ){
        int N = P.rows();
        Eigen::Matrix<Scalar,1,2> a = P.row((i-1+N)%N);
        Eigen::Matrix<Scalar,1,2> b = P.row((i+1)%N);
        bool cross = false;
        for(int e=0;e<P.rows();e++){
            Eigen::Matrix<Scalar,1,2> c(2);
            Eigen::Matrix<Scalar,1,2> d(2);
            Eigen::Matrix<Scalar,1,2> _q;
            c<<P(e,0),P(e,1);
            d<<P((e + 1)%N,0),P((e + 1)%N,1);
            if(segment_segment_intersect(a,b,c,d,_q,false)){
                cross = true;
                break;
            }
        }
        if(!cross){
            Q.resize(N-1,2);
            Q<<P.topRows(i),P.bottomRows(N-1-i);
        }
        return !cross;
    };
    while(!edges_to_check.empty()){
        edgesTmp.clear();
        double loop = 1;
        for(auto edge: edges_to_check){
            std::cerr<<"simplify "<<loop<<"/"<<edges_to_check.size()<<"     ";
            std::cerr<<"\r";
            int a = edge.first.first;
            int b = edge.first.second;
            std::pair<int,int> e1(a,b);
            std::pair<int,int> e2(b,a);
            
            int f1 = pMap.edgeToFace[e1];
            int f2 = pMap.edgeToFace[e2];
            int c = Fn.row(f1).sum()-a-b;
            int d = Fn.row(f2).sum()-b-a;

            Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> p1,p2;
            std::vector<int> _ind1,_ind2;
            create_poly_from_edges(uv_vec,pMap,{e1,std::pair<int,int>(b,c),std::pair<int,int>(c,a)},_ind1,p1);
            create_poly_from_edges(uv_vec,pMap,{e2,std::pair<int,int>(a,d),std::pair<int,int>(d,b)},_ind2,p2);
            int len = pMap.vpos_s[e1].size();
            for (int i = 1; i < len - 1; ++i){
                Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> q1;
                Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> q2;
                if(check_simple(p1,i,q1) && check_simple(p2,len-1-i,q2) && orientation(q1)>0 && orientation(q2)>0){
                    p1 = q1;
                    p2 = q2;
                    pMap.drop_point(e1,i);
                    pMap.new_v--;
                    // no need for -> pMap.drop_point(e1,n-1-i);
                    len--;
                    i--;
                    for(auto p: {std::pair<int,int>(b,c),std::pair<int,int>(c,a),std::pair<int,int>(a,d),std::pair<int,int>(d,b)}){
                        if(pMap.edgeToFace.count(p)!=0){
                            edgesTmp[p]=pMap.edgeToFace[p];
                        }
                    }
                }
            }
            loop++;
        }
        edges_to_check=edgesTmp;    
    }
}
template <typename DerivedV>
void create_poly_from_edges(
    const std::vector<std::pair<typename DerivedV::Scalar,typename DerivedV::Scalar>>& uv_vec,
    PointMap<typename DerivedV::Scalar>& pMap,
    std::vector<std::pair<int,int>> edges,
    std::vector<int>& indices,
    Eigen::PlainObjectBase<DerivedV>& P
){
    typedef typename DerivedV::Scalar Scalar;

    P.resize(0,0);
    for(int i=0;i<edges.size();i++){
        auto e_p = pMap.vpos_t[edges[i]];
        auto e_i = pMap.vid_on_edge[edges[i]];
        
        if(e_p.size()<2)
            e_p = {std::make_tuple(std::get<0>(uv_vec[edges[i].first]),std::get<1>(uv_vec[edges[i].first]),0.0),
                   std::make_tuple(std::get<0>(uv_vec[edges[i].second]),std::get<1>(uv_vec[edges[i].second]),0.0)};
        int init_s = P.rows();
        P.conservativeResize(P.rows()+e_p.size()-1,2);
 
        for(int j=0;j<e_p.size()-1;j++){ // avoid overlap end point
            P.row(init_s+j)<<std::get<0>(e_p[j]),std::get<1>(e_p[j]);
        }
        
        indices.push_back(edges[i].first);
        for(int k=0;k<e_i.size();k++)
            indices.push_back(e_i[k]);
    }
}

// use shor algorithm to triangualte the refined triangle
template <typename Scalar>
bool triangulate_face(
    int triIndex, 
    Eigen::MatrixXi& Fn, 
    PointMap<Scalar>& pMap, 
    const std::vector<std::pair<Scalar,Scalar>>& uv_vec
){

    std::pair<int,int> e1(Fn(triIndex,0), Fn(triIndex,1));
    std::pair<int,int> e2(Fn(triIndex,1), Fn(triIndex,2));
    std::pair<int,int> e3(Fn(triIndex,2), Fn(triIndex,0));
    
    if (e1.first == -1)
        return (true);    // already triangulated
    
    std::vector<int> indices;
    //Polygon_2 poly;
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> P;
    create_poly_from_edges(uv_vec,pMap,{e1,e2,e3},indices,P);

    if (indices.size() == 3)
        return (true); // no need to triangulate
    Eigen::MatrixXi eF;
    if(!is_simple_polygon(P)){
       Eigen::VectorXi local_R;
       Eigen::MatrixXd nP;
       local_R.setZero(P.rows());
       bool succ = Shor_van_wyck(P,local_R,"",nP,eF,false);
       if(!succ){
           std::cout<<std::setprecision(17)<<P<<std::endl;
       }
    }else{
	    Eigen::VectorXi R(P.rows());
        R.setZero();
	    Eigen::VectorXi C_,R_;
	    Eigen::MatrixXd nP_;
        ear_clipping(P,R,C_,eF,nP_,R_);
    }
    int nf=Fn.rows();
    Fn.conservativeResize(Fn.rows()+eF.rows(),Eigen::NoChange);
    for(int i=0; i<eF.rows(); i++){
        Fn.row(nf++)<<indices[eF(i,0)], indices[eF(i,1)], indices[eF(i,2)];
    }
    Fn.row(triIndex)<<-1,-1,-1;    
    return (true);
}

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
){
    typedef typename DerivedV::Scalar Scalar;
    typedef std::tuple<Scalar,Scalar,Scalar> P3;

    std::cout<<"refine "<<triIndex<<std::endl;
    for (int j = 0; j < 3; ++j)    //for each edge in the triangle find the intersections
    {
        int v1 = Fs(triIndex,j);
        int v2 = Fs(triIndex,(j+1)%3);
        std::pair<int,int> e(v1, v2);
        // edge already refined or edge is a boundary        
        if(pMap.vpos_s.count(e) != 0 || pMap.vpos_s.count(std::pair<int,int>(e.second,e.first)) != 0 || FF_s(triIndex,j) == -1){
            continue;
        }
        
        // robustly find intersection points
        Eigen::Matrix<Scalar,1,Eigen::Dynamic> o = vh_s.row(v1);
        Eigen::Matrix<Scalar,1,Eigen::Dynamic> t = vh_s.row(v2);
        
        Eigen::PlainObjectBase<DerivedV> IT;
        Eigen::Matrix<Scalar,1,Eigen::Dynamic> q,nq=o;
        bool found_intersection = true;
        std::unordered_map<std::pair<int,int>, bool> check;
        // find first face
        std::vector<int> fl;
        auto l = t_tree.find(vh_t,Ft,o);
        do{
            found_intersection = false;
            for(int k=0;k<l.size();k++){
                for(int i=0;i<3;i++){
                    Eigen::Matrix<Scalar,1,Eigen::Dynamic> a = vh_t.row(Ft(l[k],i));
                    Eigen::Matrix<Scalar,1,Eigen::Dynamic> b = vh_t.row(Ft(l[k],(i+1)%3));
                    if(segment_segment_intersect(a,b,nq,t,q,true) && check.count(std::pair<int,int>(Ft(l[k],i),Ft(l[k],(i+1)%3)))==0){
                        check[std::pair<int,int>(Ft(l[k],i),Ft(l[k],(i+1)%3))] = true;
                        check[std::pair<int,int>(Ft(l[k],(i+1)%3),Ft(l[k],i))] = true;
                        fl.push_back(l[k]);
                        if(q == t || q == o){ // q==o happens when ab and ot overlap
                            break;
                        }
                        IT.conservativeResize(IT.rows()+1,2);
                        IT.bottomRows(1) << q(0),q(1);
                        nq = q;
                        // special case: if edge (nq,t) passes exactly a or b
                        if(q == a || q == b){
                            l = t_tree.find(vh_t,Ft,q);
                        }else{
                            l = {FF_t(l[k],i)};
                        }
                        found_intersection = true;
                        break;
                    }
                }
            }
        }while(found_intersection);
        pMap.add_point(e,std::make_tuple(get<0>(uv_vec[v1]),get<1>(uv_vec[v1]),0.0),
                         std::make_tuple(V(v1,0),V(v1,1),V(v1,2)));

        Eigen::PlainObjectBase<DerivedV> Ts,Tt; // intersection points
        get_image(vh_s,Fs,IT,s_tree,V,Ts);
        get_image(vh_t,Ft,IT,t_tree,Vd,Tt);
        std::vector<int> v_ids;
        // for each new point find the uv cord , and the point to refine in the original source mesh
        for (int k = 0; k < IT.rows(); ++k){
            P3 p_t = std::make_tuple(Tt(k,0),Tt(k,1),0);
            P3 p_s = std::make_tuple(Ts(k,0),Ts(k,1),Ts(k,2));
            pMap.add_point(e,p_t,p_s);
            pMap.new_v++;
            v_ids.push_back(uv_vec.size());
            uv_vec.push_back(std::make_pair(Tt(k,0),Tt(k,1)));
        }
        pMap.add_point(e,std::make_tuple(get<0>(uv_vec[v2]),get<1>(uv_vec[v2]),0.0),
                         std::make_tuple(V(v2,0),V(v2,1),V(v2,2)));
        pMap.set_index(e,v_ids);
        //now we found the face index of the edge and the opposite edge
        pMap.edgeToFace[e] = triIndex;
        for(int k=0;k<3;k++){
            int fid = FF_s(triIndex,k);
            if(fid == -1)
                continue;
            else{
                for(int s=0;s<3;s++){
                    if(Fs(fid,s) == v2 && Fs(fid,(s+1)%3)== v1){
                        pMap.edgeToFace[std::pair<int,int>(e.second,e.first)] = fid;
                        k=3; // break outer loop
                        break;
                    }
                }
            }
        }
    }
}

// Assume mesh (V,F) is disk topology!
void locally_injective_map(
    Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& P,
    const Eigen::VectorXi& R,
    const Eigen::VectorXi& T,
    Eigen::MatrixXd& uv,
    Eigen::MatrixXi& Fn
){
    
    // [triangulate target domain]
    Eigen::MatrixXd V2;
    Eigen::MatrixXi F2;
    Shor_van_wyck(P,R,"",V2,F2,true);
    
    // [using tutte embedding generate common domain]
    Eigen::VectorXi B1,B2;
    igl::boundary_loop(F ,B1);
    igl::boundary_loop(F2,B2);
    // matching B1 and B2
    // ...
    Eigen::MatrixXd circle;
    Eigen::MatrixXd H1,H2; // two harmonic map as common domain
    igl::map_vertices_to_circle(V,B1,circle);
    igl::harmonic(F ,B1,circle,1,H1);
    igl::harmonic(F2,B2,circle,1,H2);

    // [generate joint common map]
    generate_map(F,F2,H1,H2,V,V2,Fn,uv);
    std::cout<<"F size: "<<Fn.rows()<<"x"<<Fn.cols()<<std::endl;
    std::cout<<"V size: "<<V.rows()<<"x"<<V.cols()<<std::endl;
    std::cout<<"max v index: "<<Fn.maxCoeff()<<std::endl;
    std::cout<<"min v index: "<<Fn.minCoeff()<<std::endl;
}
