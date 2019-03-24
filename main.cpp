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
    //test_shor();
    // test_angles();
    // exit(0);
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
    std::string model_name, ri_file, poly_file, model_name2;
    cmdl("-in") >> model_name;
    cmdl("-r") >> ri_file;
    cmdl("-poly") >> poly_file;
    cmdl("-in2") >> model_name2;

    Eigen::MatrixXd V,uv;
    Eigen::MatrixXi F,Fuv;
    Eigen::MatrixXd P;
    Eigen::VectorXi R,match;

    State state = State();
    // assume it to be disk
    state.load_mesh(model_name,V,F,uv,Fuv);
    state.load_polygon(poly_file,P);
    Eigen::MatrixXd uv2;
    Eigen::MatrixXi Fuv2;
    state.load_mesh(model_name2,uv2,Fuv2); // target mesh

    Eigen::VectorXi bd;
    igl::boundary_loop(Fuv2,bd);
    
    std::map<int,int> Rmap;
    state.load_rotation_index(ri_file,Rmap);
    igl::boundary_loop(Fuv,match);
    uv2.conservativeResize(uv2.rows(),2);
    for(int i=0;i<match.rows();i++){
        //uv2.row(bd(i))<<uv.row(match((i+8)%match.rows()));
        //uv2.row(bd((i+8)%bd.rows()))<<uv.row(match(i));
    }
    igl::opengl::glfw::Viewer vr;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    vr.plugins.push_back(&menu);
    vr.data().set_mesh(uv2,Fuv2);
    static int pos = 0;
    menu.callback_draw_viewer_menu = [&](){
        // Add new group
        menu.draw_viewer_menu();
        if(ImGui::CollapsingHeader("polygon", ImGuiTreeNodeFlags_DefaultOpen)){
            if (ImGui::InputInt("vid", &pos)){
                ;
            }
            if (ImGui::Button("display", ImVec2(-1,0))){
                vr.data().clear();
                vr.data().set_mesh(uv2,Fuv2);
                vr.data().add_points(uv2.row(pos),Eigen::RowVector3d(1,0,0));
            }
        }
    };
    vr.launch();
    R.setZero(match.rows());
    set_rotation_index(uv2,Fuv2,R);
    // myfile.close();
    // R.setZero(match.rows());
    // R.setZero(P.rows());
    // for(int i=0;i<P.rows();i++){
    //     R(i) = Rmap[i];
    // }
    //R(1270)=0;
    igl::slice(uv2,bd,1,P);
    std::vector<int> kt;
    for(int i=0;i<R.rows();i++)
        if(R(i))
            kt.push_back(i);
    //std::cout<<std::setprecision(17)<<P<<std::endl;
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
