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
    viewer.data().set_mesh(uv,Fn);
    viewer.append_mesh();
    viewer.data().set_mesh(V,Fn);
    unsigned int left_view, right_view;
    int cube_id = viewer.data_list[0].id;
    int sphere_id = viewer.data_list[1].id;
    viewer.callback_init = [&](igl::opengl::glfw::Viewer &)
    {
      viewer.core().viewport = Eigen::Vector4f(0, 0, 640, 800);
      left_view = viewer.core_list[0].id;
      right_view = viewer.append_core(Eigen::Vector4f(640, 0, 640, 800));
      return false;
    };

    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
    {
      if(key == GLFW_KEY_SPACE)
      {
        // By default, when a core is appended, all loaded meshes will be displayed in that core.
        // Displaying can be controlled by calling viewer.data().set_visible().
        viewer.data(cube_id).set_visible(false, left_view);
        viewer.data(sphere_id).set_visible(false, right_view);
        viewer.core_list[0].align_camera_center(uv);
        viewer.core_list[1].align_camera_center(V);
      }
      return false;
    };

    viewer.callback_post_resize = [&](igl::opengl::glfw::Viewer &v, int w, int h) {
      v.core( left_view).viewport = Eigen::Vector4f(0, 0, w / 2, h);
      v.core(right_view).viewport = Eigen::Vector4f(w / 2, 0, w - (w / 2), h);
      return true;
    };
    //viewer.data().set_face_based(true);
    viewer.launch();
}
