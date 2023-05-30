#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "parm.h"

parm::parm(std::string parm7file){
    std::fstream nfile;
    nfile.open(parm7file, std::ios::in);
    if(nfile.is_open()){
        std::string temp;
        std::string current_flag;
        while(getline(nfile, temp)){
            if(temp.substr(0,5) == "%FLAG"){
                std::istringstream splitstring(temp);
                splitstring >> current_flag;
                splitstring >> current_flag;
                std::cout << current_flag << std::endl;

            }
            if(current_flag == "ATOM_NAME" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                std::string value;
                for(int i = 0; i < temp.size(); i+= 4){
                    std::cout << temp.substr(i, 4) << std::endl;
                    ATOM_NAME.push_back(temp.substr(i, 4));
                }
            }
            else if(current_flag == "CHARGE" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    CHARGE.push_back(value);
                }
            }
           // else if(current_flag == "")
        }


    }
}