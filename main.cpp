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
    test_is_simple_polygon();
    
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
