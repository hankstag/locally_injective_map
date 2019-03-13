#include <igl/opengl/glfw/Viewer.h>
#include <igl/harmonic.h>
#include <igl/euler_characteristic.h>
#include <igl/boundary_loop.h>
#include "argh.h"
#include "shor.h"

// void get_matching_from_model(
//     const std::string model,
//     Eigen::VectorXi& match
// ){
//     Eigen::MatrixXd CN,V,uv;
//     Eigen::MatrixXi Fuv,FN,F;
//     igl::readOBJ(model,V,uv,CN,F,Fuv,FN);

//     // copy face from uv
//     std::vector<std::vector<int>> D1,D2;
//     igl::boundary_loop(F  ,D1);
//     igl::boundary_loop(Fuv,D2);
//     std::cout<<"model size: #F("<<F.rows()<<"), #V("<<V.rows()<<")"<<std::endl;
//     std::cout<<"source #BD("<<D1.size()<<")"<<std::endl;
//     std::cout<<"target #BD("<<D2.size()<<")"<<std::endl;
//     #define USEUV
//     #ifdef USEUV
//     Eigen::MatrixXd nV(uv.rows(),3);
//     for(int i=0;i<F.rows();i++){
//         nV.row(Fuv(i,0)) << V.row(F(i,0));
//         nV.row(Fuv(i,1)) << V.row(F(i,1));
//         nV.row(Fuv(i,2)) << V.row(F(i,2));
//     }
//     F = Fuv;
//     V = nV;
//     #endif
//     //fill_in_holes(V,F);
//     auto D = D2[0];
//     Eigen::MatrixXd P;
//     P.resize(D.size(),2);
//     match.resize(P.rows());
//     for(int i=0;i<P.rows();i++){
//         P.row(i)<<uv.row(D[i]);
//         match(i) = D[i];
//     }
// }

int test_flip(
    const Eigen::MatrixXd& V, 
    const Eigen::MatrixXi& F
){
    bool flipped = false;
    int count = 0;
    for(int i=0;i<F.rows();i++){
        if(F.row(i).sum()==0) continue;
        Eigen::Matrix<double,3,2> tri;
        tri<<V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2));
        if(orientation(tri)<=0) {
            count++;
            flipped = true;
        }
    }
    return count;
}

int main(int argc, char *argv[])
{
    auto cmdl = argh::parser(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);
    if(cmdl[{"-h","-help"}]){
        std::cout<<"Usage: ./tutte_test -options"<<std::endl;
        std::cout<<"-uv: uv model"<<std::endl;
        exit(0);
    }
    int loop, threshold;
    bool extract_bd;
    std::string uv_model;
    cmdl("-uv") >> uv_model;
    Eigen::MatrixXd V;
    Eigen::MatrixXd uv;
    Eigen::MatrixXi F;
    Eigen::MatrixXd CN;
    Eigen::MatrixXi Fuv,FN;
    igl::readOBJ(uv_model,V,uv,CN,F,Fuv,FN);
    Eigen::VectorXi bd,R,cp;
    Eigen::MatrixXd P,C;
    Eigen::MatrixXi E;
    std::vector<std::vector<int>> L;
    igl::boundary_loop(F,bd);
    V.conservativeResize(V.rows(),2);
    igl::slice(V,bd,1,P);
    set_rotation_index(V,F,R);
    bool x = test_flip(V, F);
    std::cout<<"flipped "<<x<<std::endl;
    bool succ = Shor_van_wyck(P,R,"",V,F,true);
    //bool succ = Shor_van_wyck(P,R,V,F);
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
