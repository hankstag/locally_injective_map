#ifndef STATE_H
#define STATE_H

#include <Eigen/Core>
#include <map>

class State{

public:

    State(){};
    ~State(){};

    void load_mesh(const std::string fname, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
    void load_mesh(const std::string fname, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& uv, Eigen::MatrixXi& Fuv);
    void load_polygon(const std::string& fileID, Eigen::MatrixXd& polygon);
    void load_rotation_index(const std::string& fname, std::map<int,int>& Rmap);

};

#endif