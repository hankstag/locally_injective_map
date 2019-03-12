#include "shor.h"
#include "is_inside_polygon.h"
#include <igl/boundary_loop.h>
#include <igl/triangle_triangle_adjacency.h>
#include "util.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/AABB.h>
#include <igl/in_element.h>
#include <igl/copyleft/cgal/orient2D.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

using namespace std;

// set rotation index from initial parametrization
void set_rotation_index(
    const Eigen::MatrixXd& uv,
    const Eigen::MatrixXi& F, 
    Eigen::VectorXi& R
){
    auto get_angle = [&](
        int v,
        int f
    ){
        double angle = -1.0f;
        for(int i=0;i<3;i++){
            if(F(f,i)==v){
                Eigen::Vector2d a = uv.row(F(f,(i+2)%3));
                Eigen::Vector2d b = uv.row(v);
                Eigen::Vector2d c = uv.row(F(f,(i+1)%3));
                auto bc = c - b;
                auto ba = a - b;
                angle = std::acos(bc.dot(ba)/(bc.norm()*ba.norm()));
                break;
            }
        }
        return angle;
    };
    using namespace std;
    Eigen::VectorXi B;
    igl::boundary_loop(F,B);
    R.setZero(B.rows());
    vector<vector<int>> M,MI;
    igl::vertex_triangle_adjacency(uv,F,M,MI);
    for(int i=0;i<B.rows();i++){
        int v = B(i);
        double sum = 0;
        for(int f: M[v]){
            sum += get_angle(v,f);
        }
        int ri = static_cast<int>(std::floor(sum / igl::PI / 2));
        if(ri>0){
            std::cout<<i<<","<<sum*180/igl::PI<<std::endl;
            R(i) = ri;
        }
    }
}

void add_triangle(
	Eigen::MatrixXi& F, 
	int i, int j, 
	std::vector<std::vector<int>>& K,
	std::vector<std::vector<int>>& Q
){
	if(K[i][j]==-1) return;
	int k = K[i][j];
	F.conservativeResize(F.rows()+1,3);
	F.bottomRows(1) << i,k,j;
	add_triangle(F,i,k,K,Q);
	add_triangle(F,k,j,K,Q);
};

// for debug shor

bool weakly_self_overlapping(
	const Eigen::MatrixXd& P,
	const Eigen::VectorXi& R,
	Eigen::MatrixXi& F
){
    std::cout<<"weakly self overlapping test"<<std::endl;
	auto add_to_table = [](
		const Eigen::MatrixXd& P,
		const Eigen::VectorXi& R,
		std::vector<std::vector<Angle>>& FA,
		std::vector<std::vector<Angle>>& LA,
		int i, int k, int j
	){
		const int N = P.rows();
		Eigen::Matrix<double,3,2> tri;
		tri<<P.row(i),P.row(k),P.row(j);
		if(orientation(tri)<=0) return false;
		
		Angle I(tri.row(2),tri.row(0),tri.row(1)); // J I K
		Angle J(tri.row(1),tri.row(2),tri.row(0)); // K J I
		
		Angle k1 = LA[i][k];				
		Angle k2(P.row(i),P.row(k),P.row(j));
		Angle k3 = FA[k][j];
		Angle Fr = I + FA[i][k];
		Angle La = LA[k][j] + J;

		if((Fr.r <= R(i)) && (La.r <= R(j)) && (k1+k2+k3).r == R(k)){
			FA[i][j] = Fr;
			LA[i][j] = La;
			return true;
		}else 
			return false;
	};

	const int N = P.rows();
	std::vector<std::vector<int>> K(N, std::vector<int>(N,-1));
	std::vector<std::vector<int>> Q(N, std::vector<int>(N,0));
	std::vector<std::vector<Angle>> FA(N, std::vector<Angle>(N));
	std::vector<std::vector<Angle>> LA(N, std::vector<Angle>(N));
	bool succ = false;
	int i;
	for(int i=0;i<N;i++){
		Q[i][(i+1)%N]=1;
		FA[i][(i+1)%N] = Angle(P.row((i+1)%N),P.row(i),P.row((i+1)%N));
		LA[i][(i+1)%N] = Angle(P.row(i),P.row((i+1)%N),P.row(i));
	}
    int h = -1;
    std::vector<std::vector<double>> A(N, std::vector<double>(N,-1));
	for(int d=2;d<N && !succ;d++){
		for(i=0;i<N && !succ;i++){
			int j = (i + d)%N;
			int k = (i + 1)%N;
            double min_q = -1.0f;
			while (k!=j){
				if (Q[i][k]==1 && Q[k][j]==1){
					if(add_to_table(P,R,FA,LA,i,k,j)){
						Q[i][j] = 1;
						K[i][j] = k;
						if(i == (j+1)%N){
							succ = true;
							h = i;
							break;
						}
					}
				}
				k = (k + 1)%N;
			}
		}
	}
    std::cout<<"test done"<<std::endl;    
	if(h==-1) return false;
	add_triangle(F,h,(h-1+N)%N,K,Q);
	return true;
}

bool is_convex(
    const Eigen::MatrixXd& P
){
    for(int i=0;i<P.rows();i++){
        int prev = (i-1+P.rows())%P.rows();
        int next = (i+1)%P.rows();
        double a[2] = {P(prev,0),P(prev,1)};
        double b[2] = {P(i,0),P(i,1)};
        double c[2] = {P(next,0),P(next,1)};
        short r = igl::copyleft::cgal::orient2D(a,b,c);
        if(r < 0)
            return false;
    }
    return true;
}

void drop_colinear(
	const Eigen::MatrixXd& P,
	const Eigen::VectorXi& R,
	Eigen::VectorXi& B,
	Eigen::MatrixXd& mP,
	Eigen::VectorXi& mR
){
	int dropped = 0;
	int N = P.rows();
	std::vector<int> Bv;
	for(int i=0;i<N;i++){
		int a = (i-1+N)%N;
		int v = i;
		int b = (i+1)%N;
		Eigen::Matrix<double,3,2> T;
		T<<P.row(a),P.row(v),P.row(b);
		if(R(v)!=0 || orientation(T)!=0){ // non-colinear or rotate index nonzero
			Bv.push_back(v);
		}else
			dropped++;
	}
	igl::list_to_matrix(Bv,B);
	igl::slice(P,B,1,mP);
	igl::slice(R,B,mR);
}

void add_colinear(
	const Eigen::MatrixXd& P,
	const Eigen::MatrixXi& nF,
	const Eigen::VectorXi& B,
	Eigen::MatrixXi& F
){
	int added = 0;
	int sN = nF.maxCoeff()+1; // size of simplified polygon
	F.resizeLike(nF);
	for(int i=0;i<F.rows();i++)
		F.row(i)<<B(nF(i,0)),B(nF(i,1)),B(nF(i,2));
	typedef std::pair<int,int> Edge;
	std::vector<Edge> stk;
	std::set<std::pair<int,int>> sid; // edge vertex pairs needs split
	for(int i=0;i<nF.rows();i++){
		bool split = false;
		for(int k=0;k<3;k++){
			int a = nF(i,k);
			int b = nF(i,(k+1)%3);
			if(	(b-a==1 || a-b == sN-1) &&
				(B(b)-B(a)!=1 && B(a)-B(b)!=P.rows()-1)){ // edge need split
				if(!split){ // do not split a face twice (do it later)	
					stk.push_back(Edge(i,k));
					split = true;
				}
				sid.insert(std::make_pair(B(a),B(b)));
			}
		}
	}
	while(!stk.empty()){
		Edge e = stk.back();
		stk.pop_back();
		int a = F(e.first,e.second);
		int b = F(e.first,(e.second+1)%3);
		int c = F(e.first,(e.second+2)%3);
		int rg = a < b ? b-a : b+P.rows()-a;
		int r = F.rows();
		F.conservativeResize(F.rows()+rg,3);
		added += (rg-1);
		// from r+0 to r+rg-1
		for(int i=0;i<rg;i++){
			F.row(r+i)<<(a+i)%P.rows(),(a+i+1)%P.rows(),c;
		}
		// if bc needs split
		if(sid.find(std::make_pair(b,c))!=sid.end()){
			stk.push_back(Edge(r+rg-1,1)); // its now 2nd edge of face r+rg-1
		}
		// if ac needs split
		if(sid.find(std::make_pair(c,a))!=sid.end()){
			stk.push_back(Edge(r,2)); // its now 3rd edge of face r
		}
		F.row(e.first) << -1,-1,-1; // delete face
	}
	Eigen::MatrixXi tF=F;
	// drop the (-1,-1,-1) rows
	int k=0;
	for(int i=0;i<F.rows();i++)
		if(F.row(i).sum()!=-3)
			tF.row(k++)<<F.row(i);
	tF.conservativeResize(k,3);
	F = tF;
}

void subdivide_polygon(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	std::vector<std::vector<int>>& L
){
	// collect halfedge and corresponding face id
	// initialize L to be collection of all faces(polygons)
	std::map<std::pair<int,int>, int> H;
	Eigen::VectorXi G(F.rows());
	Eigen::MatrixXi FF,FFI;
    igl::triangle_triangle_adjacency(F,FF,FFI);
	L.resize(F.rows());
	for(int i=0;i<F.rows();i++){
		G(i) = i; // every face belongs to itself
		for(int k=0;k<3;k++){
			L[i].push_back(F(i,k));
			H[std::make_pair(F(i,k),F(i,(k+1)%3))] = FF(i,k);
		}
	}

	// traverse all the halfedges
	for(auto h: H){ // [he, pid]
		auto he = h.first;
		auto rhe = std::make_pair(he.second,he.first);
		int p1 = H[he]; // a -> b
		int p2 = H[rhe]; // b -> a
		if(p1 == -1 || p2 == -1) continue;
		
		// up to group root
		while(p1!=G[p1])
			p1 = G[p1];
		while(p2!=G[p2])
			p2 = G[p2];

		// combine p1 and p2
		Eigen::MatrixXd poly(L[p1].size()+L[p2].size()-2,2);
		auto a = std::find(L[p1].begin(),L[p1].end(),he.first);
		auto b = std::find(L[p2].begin(),L[p2].end(),he.second);

		std::vector<int> L1(L[p1].begin(),a);
		std::vector<int> R1(a,L[p1].end());
		std::vector<int> L2(L[p2].begin(),b);
		std::vector<int> R2(b,L[p2].end());
		
		std::vector<int> S;
		S = L1;
		auto c = R2.empty() ? R2.begin() : R2.begin()+1;
		auto d = R1.empty() ? R1.begin() : R1.begin()+1;
		S.insert(S.end(),c,R2.end());
		S.insert(S.end(),L2.begin(),L2.end());
		S.insert(S.end(),d,R1.end());
	    
		for(int i=0;i<poly.rows();i++){
			poly.row(i)<<V.row(S[i]);
        }
		
		// if merged polygon is simple, drop edge he/rhe
		// erase L[p2], add to L[p1]
	    if(is_simple_polygon(poly)){
			H[he]=-1;
			H[rhe]=-1;
			G[p2] = p1; // p2 belongs to p1 now
			L[p1] = S;
			L[p2].clear();
		}
	}
}

void simplify_triangulation(
	const Eigen::MatrixXd& V_i,
    const Eigen::MatrixXd& C,
    const Eigen::VectorXi& cp,
    const Eigen::MatrixXi& E,
	const std::vector<std::vector<int>>& L_i,
    const std::string& flags,
	Eigen::MatrixXd& V,
	Eigen::MatrixXi& F
){
	std::vector<std::vector<int>> L;
	V = V_i;
	F.resize(0,0);
	std::map<std::pair<int,int>, std::vector<int>> Em;
    std::set<std::pair<int,int>> cnnt; // gather connectivity info from E
    std::vector<bool> embed(C.rows(),false); // make sure every C is embeded once

    for(int i=0;i<E.rows();i++){
        cnnt.insert(std::make_pair(E(i,0),E(i,1)));
        cnnt.insert(std::make_pair(E(i,1),E(i,0)));
    }

    // append C to V
    if(C.rows()>0){
        V.conservativeResize(V.rows()+C.rows(),Eigen::NoChange);
        V.bottomRows(C.rows())<<C;
    }

	// calculate the average edge length
	double avl = 0.0f;
	for(int i=0;i<V_i.rows();i++){
		avl += (V_i.row((i+1)%V_i.rows())-V_i.row(i)).norm();
	}
	avl /= V_i.rows();
	// traverse L, split the internal edges
    std::vector<Eigen::MatrixXd> P_set;
    for(int i=0;i<L_i.size();i++){
        if(L_i[i].size()==0) continue;
        Eigen::MatrixXd LLP;
        LLP.resize(L_i[i].size(),2);
        for(int j=0;j<LLP.rows();j++){
            LLP.row(j)<<V.row(L_i[i][j]);
        }
        P_set.push_back(LLP);
    }
    auto count_split = [&](const std::vector<int>& L, int i, std::vector<std::pair<double,double>>& Vc){
        int a = L[i];
        int b = L[(i+1)%L.size()];
        std::vector<int> ops = {a,b};
        // compare left and right edge of [i,i+1]
        double l = (V.row(a) - V.row(b)).norm();
        int op1 = a;
        int op2 = b;
        double l1 = (V.row(L[(i+2)%L.size()]) - V.row(b)).norm();
        double l2 = (V.row(L[(i-1+L.size())%L.size()]) - V.row(a)).norm();
        while(l2/l<1.0){
            double t1,t2;
            if(!Vc.empty()){
                t1 = (V(a,0)+Vc.back().first )/2;
                t2 = (V(a,1)+Vc.back().second)/2;
            }else{
                t1 = (V(a,0)+V(b,0))/2;
                t2 = (V(a,1)+V(b,1))/2;
            }
            Vc.push_back(std::make_pair(t1,t2));
            Eigen::RowVectorXd t(2);
            t<<t1,t2;
            l = (V.row(a) - t).norm();
        }
        std::reverse(Vc.begin(),Vc.end());
        Eigen::RowVectorXd p(2);
        if(Vc.size()!=0){
            p<<Vc.back().first,Vc.back().second;
            l = (V.row(b) - p).norm();
        }else{
            l = (V.row(a) - V.row(b)).norm();
        }
        while(l1/l<1.0){
            double t1,t2;
            if(!Vc.empty()){
                t1 = (V(b,0)+Vc.back().first )/2;
                t2 = (V(b,1)+Vc.back().second)/2;
            }else{
                t1 = (V(a,0)+V(b,0))/2;
                t2 = (V(a,1)+V(b,1))/2;
            }
            Vc.push_back(std::make_pair(t1,t2));
            Eigen::RowVectorXd t(2);
            t<<t1,t2;
            l = (V.row(b) - t).norm();
        }
    };

    std::map<std::pair<int,int>,std::pair<int,int>> EP; // edge to polygon list map
    for(int i=0;i<L_i.size();i++){
        for(int j=0;j<L_i[i].size();j++){
            int a = L_i[i][j];
			int b = L_i[i][(j+1)%L_i[i].size()];
            EP[std::make_pair(a,b)] = std::make_pair(i,j);
        }
    }
    std::vector<std::vector<int>> L_aux(L_i.size());
	for(int i=0;i<L_i.size();i++){
		if(!L_i[i].empty()){
			L.resize(L.size()+1); // conservativeResize
		}
		for(int j=0;j<L_i[i].size();j++){
			int a = L_i[i][j];
			int b = L_i[i][(j+1)%L_i[i].size()];
			L_aux[i].push_back(a);
			if(Em.find(std::make_pair(b,a))!=Em.end()){
				auto tv = Em[std::make_pair(b,a)]; // tmp vec
				std::reverse(tv.begin(),tv.end());
				Em[std::make_pair(a,b)] = tv;
				L_aux[i].insert(L_aux[i].end(),tv.begin(),tv.end());
				continue;
			}
			int vn = V.rows();
			if(b-a!=1 && a-b!=V_i.rows()-1){
                //#define REFINE_EDGE
                #ifndef REFINE_EDGE
                int n = std::max((V.row(b)-V.row(a)).norm() / avl,1.0);
                V.conservativeResize(V.rows()+n-1,2);
                for(int k=0;k<n-1;k++){
					V.row(vn+k)<<V.row(a) + (V.row(b)-V.row(a))*(k+1)/n;
					Em[std::make_pair(a,b)].push_back(vn+k);
					L_aux[i].push_back(vn+k);
				}
                #else
                std::vector<std::pair<double,double>> Vc1,Vc2;
                int i_ = EP[std::make_pair(b,a)].first;
                int j_ = EP[std::make_pair(b,a)].second;
                count_split(L_i[i], j ,Vc1);
                count_split(L_i[i_],j_,Vc2);
                bool use_dual = Vc1.size() < Vc2.size();
                auto op = !use_dual ? Vc1 : Vc2;
                
                V.conservativeResize(V.rows()+op.size(),2);
                if(use_dual) std::reverse(op.begin(),op.end());
                for(int k=0;k<op.size();k++){
                    V.row(vn+k) << op[k].first, op[k].second;
                    Em[std::make_pair(a,b)].push_back(vn+k);
                    L_aux[i].push_back(vn+k);
                }
                #endif
			}
		}
	}
    for(int i=0;i<L_aux.size();i++){
        if(L_aux[i].empty()){
            L_aux.erase(L_aux.begin()+i);
            i--;
        }
    }
    L = L_aux;

	// triangulate subpolygons
	Eigen::MatrixXd LP; // local polygon 
	Eigen::MatrixXi LF; // local faces
	Eigen::MatrixXd LV;

	for(int i=0;i<L.size();i++){
		LP.resize(L[i].size(),2);
		std::unordered_set<int> LPi={-1}; // for those constraint points not in overlapping part
        for(int j=0;j<LP.rows();j++){
			LP.row(j)<<V.row(L[i][j]);
            LPi.insert(L[i][j]);
		}
        
		//display(LP,0,LP.rows()-1);
        Eigen::MatrixXi LE(LP.rows(),2);
		LE<<Eigen::VectorXi::LinSpaced(LP.rows(),0,LP.rows()-1),
			Eigen::VectorXi::LinSpaced(LP.rows(),1,LP.rows());
		LE(LE.rows()-1,1) = 0;

        // embed the extra internal points
        Eigen::MatrixXd LP2 = LP;
        std::vector<int> cin;
        for(int j=0;j<C.rows();j++){
            Eigen::RowVector2d q = C.row(j);
            if(LPi.find(cp(j))!=LPi.end() && !embed[j] && point_in_poly(LP,q)){
                LP2.conservativeResize(LP2.rows()+1,Eigen::NoChange);
                LP2.bottomRows(1) << q;
                cin.push_back(j+V_i.rows());
                embed[j] = true;
            }
        }
        // try to support edge constraints
        for(int j=0;j<cin.size();j++){
            for(int k=j+1;k<cin.size();k++){
                auto pr = std::make_pair(cin[j],cin[k]);
                if(cnnt.find(pr)!=cnnt.end()){
                    LE.conservativeResize(LE.rows()+1,Eigen::NoChange);
                    LE.bottomRows(1) << LP.rows()+j,LP.rows()+k;
                }
            }
        }
        // if(!is_simple_polygon(LP2)){
        //     std::cout<<"not simple"<<std::endl;
        //     display(LP2,0,LP2.rows()-1);
        // }
        //display(LP2,0,LP2.rows()-1);

		igl::triangle::triangulate(LP2,LE,Eigen::MatrixXd(),"YQq33",LV,LF);

		Eigen::MatrixXd nV = LV.bottomRows(LV.rows()-LP2.rows());
		int n_o = V.rows();
		V.conservativeResize(V.rows()+nV.rows(),2);
		V.bottomRows(nV.rows())<<nV;
		for(int f=0;f<LF.rows();f++){
			for(int k=0;k<3;k++){
				if(LF(f,k) >= LP2.rows())
					LF(f,k) += (n_o-LP2.rows());
                else if(LF(f,k) >= LP.rows())
                    LF(f,k) = cin[LF(f,k)-LP.rows()];
                else
					LF(f,k) = L[i][LF(f,k)];
			}
		}
		F.conservativeResize(F.rows()+LF.rows(),3);
		F.bottomRows(LF.rows())<<LF;
	}
}

// the re-implementation of Shor algorithm
bool Shor_van_wyck(
	const Eigen::MatrixXd& P,
	const Eigen::VectorXi& R,
    const std::string flags,
	Eigen::MatrixXd& V,
	Eigen::MatrixXi& F,
    bool do_refine
){

	// [drop colinear points]
	Eigen::VectorXi B; // remaining vertices: non-colinear/rotateindex!=0
	Eigen::MatrixXd mP; // P \ colinear vertices
	Eigen::VectorXi mR;
	drop_colinear(P,R,B,mP,mR);

	// [ear clipping]
	Eigen::VectorXi D;
	Eigen::MatrixXi eF;
	Eigen::MatrixXd nP;
	Eigen::VectorXi nR;
	ear_clipping(mP,mR,D,eF,nP,nR);
    // nP = mP;
    // nR = mR;
    // D = Eigen::VectorXi::LinSpaced(mP.rows(),0,mP.rows()-1);
	// [weakly-self-overlapping test (dynamic programming)]
	Eigen::MatrixXi nF;
	bool succ = (nP.rows()==0) || weakly_self_overlapping(nP,nR,nF);
	if(!succ){
        std::cout<<"shor failed"<<std::endl;
        exit(0);
        return false;
    }
    // [map simplified index to initial polygon]
	for(int i=0;i<nF.rows();i++){
		nF.row(i) << D(nF(i,0)),D(nF(i,1)),D(nF(i,2));
	}
    if(eF.rows()>0){
	    nF.conservativeResize(nF.rows()+eF.rows(),3);
	    nF.block(nF.rows()-eF.rows(),0,eF.rows(),3) = eF;
    }

	// [add back colinear vertices by spliting boundary edges]
	add_colinear(P,nF,B,F);
	if(!do_refine){
        V = P;
        return true;
    }
    V = P;
    //drop (-1,-1,-1) faces
    auto Ft = F;
    int nFt = 0;
    for(int i=0;i<F.rows();i++){
        if(F.row(i).sum()!=-3){
            Ft.row(nFt++) = F.row(i);
        }
    }
    Ft.conservativeResize(nFt,3);
    F = Ft;
	// [simplify mesh (subdivide into small polygons)]
	std::vector<std::vector<int>> L;
	subdivide_polygon(V,F,L);
    // [refine each small polygon]
    // TODO
    return true;
}
