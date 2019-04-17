#ifndef LOAD_MODEL
#define LOAD_MODEL
#include <Eigen/Core>
#include <map>

void load_model_with_seam(
    const std::string model,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::MatrixXd& polygon,
    Eigen::VectorXi& bd
);

void load_matching_info(
    std::string fname,
    std::pair<int,int>& match
);

void load_rotation_index(
    const std::string& fname,
    std::map<int,int>& Rmap
);

void load_model(
  std::string model, 
  Eigen::MatrixXd& V,
  Eigen::MatrixXd& uv,
  Eigen::MatrixXi& F,
  Eigen::MatrixXd& P,
  Eigen::VectorXi& R,
  Eigen::VectorXi& T
);

#endif