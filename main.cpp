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
#include "is_simple_polygon.h"
#include "locally_injective_map.h"

int main(int argc, char *argv[]){

    auto cmdl = argh::parser(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);
    if(cmdl[{"-h","-help"}]){
        std::cout<<"Usage: ./param_bin -options"<<std::endl;
        std::cout<<"-in: input disk model"<<std::endl;
        exit(0);
    }
    std::string model_name;
    cmdl("-in") >> model_name;

    Eigen::MatrixXd V,uv;
    Eigen::MatrixXi F;
    Eigen::MatrixXd P;
    Eigen::VectorXi R,T;
    
    load_model(model_name,V,uv,F,P,R,T);
    set_rotation_index(uv,F,R,0);
    Eigen::MatrixXi Fn;
    locally_injective_map(V,F,P,R,T,uv,Fn);
    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(uv, Fn);
    //viewer.data().set_face_based(true);
    viewer.launch();
}
