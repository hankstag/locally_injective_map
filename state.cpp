#include "state.h"
#include <string>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>

void State::load_polygon(
    const std::string& fileID,
    Eigen::MatrixXd& polygon
){
    std::ifstream myfile;
    myfile.open(fileID);
    if (!myfile.is_open()){
        std::cout << "Could not open file: " << fileID << std::endl;
        return;
    }
    double x,y;
    std::vector<std::complex<double>> polygon_vec;
    while (myfile >> x){
        myfile >> y;
        polygon_vec.push_back(std::complex<double>(x,y));
    }
    polygon.resize(polygon_vec.size(),2);
    int i=0;
    for(auto c: polygon_vec){
        polygon.row(i++)<<c.real(),c.imag();
    }
}

void State::load_rotation_index(
    const std::string& fname,
    std::map<int,int>& Rmap
){
    // Open file, and check for error
    FILE * obj_file = fopen(fname.c_str(),"r");
    if(NULL==obj_file){
        fprintf(stderr,"IOError: %s could not be opened...\n",
                fname.c_str());
        return;
    }
    #define LINE_MAX_S 2048
    char line[LINE_MAX_S];
    int line_no = 1;
    while (fgets(line, LINE_MAX_S, obj_file) != NULL){
        char type[LINE_MAX_S];
        // Read first word containing type
        if(sscanf(line, "%s",type) == 1){
            char * l = &line[strlen(type)];
            if(strlen(type) >= 1){
                if (type[0] == 'r'){ // index of point constraints
                    int x,i;
                    sscanf(l,"%d%d\n",&x,&i);
                    Rmap[x]=i;
                }
            }
        }
    }
}

void State::load_mesh(
    const std::string fname,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F
){
    auto split = [](const std::string& s, char delimiter){
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(s);
        while (std::getline(tokenStream, token, delimiter))
        {
            tokens.push_back(token);
        }
        return tokens;
    };
    std::string x = split(fname,'.').back();
    if(x == "obj"){
        Eigen::MatrixXd uv;
        Eigen::MatrixXd CN;
        Eigen::MatrixXi Fuv,FN;
        igl::readOBJ(fname,V,uv,CN,F,Fuv,FN);
    }else if(x == "off")
        igl::readOFF(fname,V,F);
}

void State::load_mesh(
    const std::string fname,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::MatrixXd& uv,
    Eigen::MatrixXi& Fuv
){
    auto split = [](const std::string& s, char delimiter){
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(s);
        while (std::getline(tokenStream, token, delimiter))
        {
            tokens.push_back(token);
        }
        return tokens;
    };
    std::string x = split(fname,'.').back();
    if(x == "obj"){
        Eigen::MatrixXd CN;
        Eigen::MatrixXi FN;
        igl::readOBJ(fname,V,uv,CN,F,Fuv,FN);
    }else if(x == "off")
        igl::readOFF(fname,V,F);
}