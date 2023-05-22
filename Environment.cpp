#include <iostream>
#include <fstream>
#include <string>
#include "class.h"

Environment::Environment(std::string pdb){
    //Acoords.reserve(30000);
    std::fstream pdb_file;
    pdb_file.open(pdb, std::ios::in);
    std::cout << "help!";
    if(pdb_file.is_open()){
        std::string temp;
        while(getline(pdb_file, temp)){
            //std::cout << temp << std::endl;
            if(temp.substr(0,4) == "ATOM"){
                append_name(temp.substr(12, 16));
                append_residue(temp.substr(17, 20));
                append_coords(temp.substr(30,38),temp.substr(38, 46),temp.substr(46, 54));
                std::cout << "atom!" << std::endl;
            }
            else if(temp.substr(0, 6) == "HETATM"){
                //class name.HETATM.append(string)
                std::cout << "hetatm!" << std::endl;
            }
        }
    }
    
    
}

void Environment::append_name(std::string a){
    Aatom_name.push_back(a);
}
void Environment::append_Hname(std::string a){
    Hatom_name.push_back(a);
}
void Environment::append_residue(std::string a){
    Aresidue.push_back(a);
}
void Environment::append_Hresidue(std::string a){
    Hresidue.push_back(a);
}
void Environment::append_coords(std::string a, std::string b, std::string c){
    float x = std::stof(a);
    float y = std::stof(b);
    float z = std::stof(c);
    std::cout << x << " " << y << " " << z << std::endl;
    Acoords.push_back(x);
    Acoords.push_back(y);
    Acoords.push_back(z);
}

void Environment::append_Hcoords(std::string a, std::string b, std::string c){
    float x = std::stof(a);
    float y = std::stof(b);
    float z = std::stof(c);
    Hcoords.push_back(x);
    Hcoords.push_back(y);
    Hcoords.push_back(z);
}
