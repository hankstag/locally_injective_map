#include <sstream>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/boundary_loop.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/slice.h>
#include "argh.h"

#include "load_model.h"
#include "shor.h"
#include "state.h"
#include "is_simple_polygon.h"
#include "locally_injective_map.h"

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
    std::string model_name, ri_file, poly_file, uvfile;
    cmdl("-in") >> model_name;
    cmdl("-uv") >> uvfile;
    cmdl("-r") >> ri_file;
    cmdl("-poly") >> poly_file;
    if(uvfile=="" || ri_file==""){
        uvfile=model_name;
        ri_file=model_name;
    }
    Eigen::MatrixXd V,uv;
    Eigen::MatrixXi F,Fuv;
    Eigen::MatrixXd P;
    Eigen::VectorXi R,T;
    Eigen::VectorXi bd0,bd1;
    std::map<int,int> Rmap;
    load_rotation_index(ri_file,Rmap);
    load_model_with_seam(model_name,V,F,P,bd0);
    // std::pair<int,int> match;
    // load_matching_info(uvfile,match);
    // Eigen::MatrixXd _polygon;
    // load_model_with_seam(uvfile,uv,Fuv,_polygon,bd1);
    // uv.conservativeResize(uv.rows(),2);
    // int id0 = 0,id1 = 0;
    // for(int i=0;i<bd0.rows();i++){
    //     if(bd0(i) == match.first)
    //         id0 = i;
    //     if(bd1(i) == match.second)
    //         id1 = i;
    // }
    // int offset = (id1-id0+bd1.rows())%bd1.rows();
    // std::cout<<"setting rotation index..."<<std::endl;
    // set_rotation_index(uv,Fuv,R,offset);
    R.setZero(P.rows());
    for(int i=0;i<P.rows();i++){
        R(i) = Rmap[bd0(i)]/360;
    }

    igl::opengl::glfw::Viewer vr;
    vr.core.align_camera_center(P);
    for(int i=0;i<P.rows();i++){
        int i_1 = (i+1) % P.rows();
        if(R(i) != 0)
            vr.data().add_points(P.row(i),Eigen::RowVector3d(0,0,0));
        vr.data().add_edges(P.row(i),P.row(i_1),Eigen::RowVector3d(1,0,0));
    }
    vr.launch();
    Eigen::MatrixXi Fn;
    locally_injective_map(V,F,P,R,T,uv,Fn);
    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(uv, Fn);
    //viewer.data().set_face_based(true);
    viewer.launch();
}
