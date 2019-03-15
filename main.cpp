#include <sstream>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/boundary_loop.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/slice.h>
#include "argh.h"

#include "shor.h"
#include "state.h"
#include "is_simple_polygon.h"
void display(const Eigen::MatrixXd& P, int s, int t, std::vector<int>& kt){
	Eigen::MatrixXi edges;
	// s to t
	int ps = s < t ? t - s + 1 : t + P.rows() - s + 1;
	Eigen::MatrixXd bp(ps,3);
	for (int i = 0; i < ps; i++)
	{
		bp(i,0)=P(i,0);
		bp(i,1)=P(i,1);
		bp(i,2)=0;
	}
	edges.resize(ps,2);
	for(int i=0;i<ps;i++){
		edges(i,0)=i;
		edges(i,1)=(i+1)%ps;
	}
	igl::opengl::glfw::Viewer viewer;
	viewer.data().clear();
	viewer.data().set_edges(bp,edges,Eigen::RowVector3d(0,0,0));
    for(int i: kt){
        viewer.data().add_points(P.row(i),Eigen::RowVector3d(1,0,0));
    }
	viewer.core.align_camera_center(bp);
	viewer.launch();
}
int main(int argc, char *argv[]){
    auto cmdl = argh::parser(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);
    if(cmdl[{"-h","-help"}]){
        std::cout<<"Usage: ./param_bin -options"<<std::endl;
        std::cout<<"-in: input disk model"<<std::endl;
        std::cout<<"-poly: target polygon"<<std::endl;
        std::cout<<"-r: rotation index info"<<std::endl;
        std::cout<<"-out: output path"<<std::endl;
        exit(0);
    }
    // test_segment_segment_intersect();
    // test_is_simple_polygon();
    
    int loop, threshold;
    bool extract_bd;
    std::string model_name, ri_file, poly_file;
    cmdl("-in") >> model_name;
    cmdl("-r") >> ri_file;
    cmdl("-poly") >> poly_file;

    Eigen::MatrixXd V,uv;
    Eigen::MatrixXi F,Fuv;
    Eigen::MatrixXd P;
    Eigen::VectorXi R,match;

    State state = State();
    // assume it to be disk
    state.load_mesh(model_name,V,F,uv,Fuv);
    state.load_polygon(poly_file,P);
    
    std::map<int,int> Rmap;
    state.load_rotation_index(ri_file,Rmap);
    // igl::boundary_loop(Fuv,match);
    // // V.conservativeResize(V.rows(),2);
    // igl::slice(uv,match,1,P);
    // // std::cout<<std::setprecision(17)<<P<<std::endl;
    // std::ofstream myfile;
    // myfile.open("genus3_poly_from_input",std::ios_base::app);
    // for(int i=0;i<P.rows();i++){
    //     myfile<<std::setprecision(17)<<P(i,0)<<" "<<P(i,1)<<std::endl;
    // }
    // myfile.close();
    // R.setZero(match.rows());
    R.setZero(P.rows());
    for(int i=0;i<P.rows();i++){
        R(i) = Rmap[i];
    }
    std::vector<int> kt;
    for(int i=0;i<R.rows();i++)
        if(R(i))
            kt.push_back(i);
    display(P,0,P.rows()-1,kt);
    bool succ = Shor_van_wyck(P,R,"",V,F,true);
    if(succ)
        std::cout<<"succ"<<std::endl;
    else
        std::cout<<"failed"<<std::endl;
    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
    viewer.launch();
}
