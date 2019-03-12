#include <igl/opengl/glfw/Viewer.h>
#include <igl/harmonic.h>
#include <igl/euler_characteristic.h>
#include <igl/boundary_loop.h>
#include "argh.h"
#include "shor.h"

int main(int argc, char *argv[])
{
    auto cmdl = argh::parser(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);
    if(cmdl[{"-h","-help"}]){
        std::cout<<"Usage: ./tutte_test -options"<<std::endl;
        std::cout<<"-in: input model name"<<std::endl;
        exit(0);
    }
    int loop, threshold;
    bool extract_bd;
    std::string model;
    cmdl("-in") >> model;
    Eigen::MatrixXd V;
    Eigen::MatrixXd uv;
    Eigen::MatrixXi F;
    Eigen::MatrixXd CN;
    Eigen::MatrixXi Fuv,FN;
    igl::readOBJ(model,V,uv,CN,F,Fuv,FN);
    int euler_char = igl::euler_characteristic(V, F);
    std::cout<<"Euler character is "<<euler_char<<std::endl;

    Eigen::VectorXi bd,R,cp;
    Eigen::MatrixXd P,C;
    Eigen::MatrixXi E;
    std::vector<std::vector<int>> L;
    igl::boundary_loop(F,bd);
    V.conservativeResize(V.rows(),2);
    igl::slice(V,bd,1,P);
    set_rotation_index(V,F,R);
    bool succ = Shor_van_wyck(P,R,C,cp,E,"",V,F,L,true);
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
