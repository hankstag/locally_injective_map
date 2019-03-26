#include <sstream>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/boundary_loop.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/slice.h>
#include "argh.h"

#include "shor.h"
#include "state.h"
#include "is_simple_polygon.h"
#include "locally_injective_map.h"

void load_all(
    std::string model, 
    Eigen::MatrixXd& V,
    Eigen::MatrixXd& uv,
    Eigen::MatrixXi& F,
    Eigen::MatrixXd& P,
    Eigen::VectorXi& R,
    Eigen::VectorXi& T
){
    Eigen::MatrixXd CN;
    Eigen::MatrixXi Fuv,FN;
    igl::readOBJ(model,V,uv,CN,F,Fuv,FN);

    // copy face from uv
    std::vector<std::vector<int>> D1,D2;
    igl::boundary_loop(F  ,D1);    
    igl::boundary_loop(Fuv,D2);
    Eigen::VectorXd A1,A2;
    igl::doublearea(V,F,A1);
    igl::doublearea(uv,Fuv,A2);
    std::cout<<"min size 3d "<<A1.minCoeff()<<", min size 2d "<<A2.minCoeff()<<std::endl;
    std::cout<<"model size: #F("<<F.rows()<<"), #V("<<V.rows()<<")"<<std::endl;
    std::cout<<"source #BD("<<D1.size()<<")"<<std::endl;
    std::cout<<"target #BD("<<D2.size()<<")"<<std::endl;
    #define USEUV
    #ifdef USEUV
    Eigen::MatrixXd nV(uv.rows(),3);
    for(int i=0;i<F.rows();i++){
        nV.row(Fuv(i,0)) << V.row(F(i,0));
        nV.row(Fuv(i,1)) << V.row(F(i,1));
        nV.row(Fuv(i,2)) << V.row(F(i,2));
    }
    F = Fuv;
    V = nV;
    #endif
    //fill_in_holes(V,F);
    auto D = D2[0];
    P.resize(D.size(),2);
    T.resize(P.rows());
    set_rotation_index(uv,Fuv,R);
    for(int i=0;i<P.rows();i++){
        P.row(i)<<uv.row(D[i]);
        T(i) = D[i];
    }
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
    
    int loop, threshold;
    bool extract_bd;
    std::string model_name, ri_file, poly_file, model_name2;
    cmdl("-in") >> model_name;
    cmdl("-r") >> ri_file;
    cmdl("-poly") >> poly_file;

    Eigen::MatrixXd V,uv;
    Eigen::MatrixXi F,Fuv;
    Eigen::MatrixXd P;
    Eigen::VectorXi R;

    State state = State();
    // assume it to be disk
    state.load_mesh(model_name,V,F,uv,Fuv);

    state.load_polygon(poly_file,P);
    Eigen::VectorXi bd;
    igl::boundary_loop(F,bd);
    std::map<int,int> Rmap;
    state.load_rotation_index(ri_file,Rmap);
    R.setZero(P.rows());
    for(int i=0;i<P.rows();i++){
        R(i) = Rmap[i];
    }
    Eigen::VectorXi T;
    if(poly_file=="")
        load_all(model_name,V,uv,F,P,R,T);
    Eigen::MatrixXi Fn;
    locally_injective_map(V,F,P,R,T,uv,Fn);
    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(uv, Fn);
    //viewer.data().set_face_based(true);
    viewer.launch();
}
