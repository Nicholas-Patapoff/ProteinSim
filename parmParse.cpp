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
            else if(current_flag == "ATOMIC_NUMBER" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG" ){
                int value;
                std::istringstream line(temp);
                while(line >> value){
                    ATOMIC_NUMBER.push_back(value);
                }
            }
            else if(current_flag == "MASS" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    MASS.push_back(value);
                }
            }
            else if(current_flag == "ATOM_TYPE_INDEX" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                int value;
                std::istringstream line(temp);
                while(line >> value){
                    ATOM_TYPE_INDEX.push_back(value);
                }
            }
            else if(current_flag == "NONBONDED_PARM_INDEX" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                int value;
                std::istringstream line(temp);
                while(line >> value){
                    NONBONDED_PARM_INDEX.push_back(value);
                }
            }
            else if(current_flag == "RESIDUE_LABEL" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                std::string value;
                std::istringstream line(temp);
                while(line >> value){
                    RESIDUE_LABEL.push_back(value);
                }
            }
            else if(current_flag == "RESIDUE_POINTER" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                int value;
                std::istringstream line(temp);
                while(line >> value){
                    RESIDUE_POINTER.push_back(value);
                }
            }
            else if(current_flag == "BOND_FORCE_CONSTANT" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    BOND_FORCE_CONSTANT.push_back(value);
                }
            }
            else if(current_flag == "BOND_EQUIL_VALUE" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    BOND_EQUIL_VALUE.push_back(value);
                }
            }
            else if(current_flag == "ANGLE_FORCE_CONSTANT" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    ANGLE_FORCE_CONSTANT.push_back(value);
                }
            }
            else if(current_flag == "ANGLE_EQUIL_VALUE" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    ANGLE_EQUIL_VALUE.push_back(value);
                }
            }
            else if(current_flag == "DIHEDRAL_FORCE_CONSTANT" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    DIHEDRAL_FORCE_CONSTANT.push_back(value);
                }
            }
            else if(current_flag == "DIHEDRAL_PERIODICITY" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    DIHEDRAL_PERIODICITY.push_back(value);
                }
            }
            else if(current_flag == "DIHEDRAL_PHASE" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    DIHEDRAL_PHASE.push_back(value);
                }
            }
            else if(current_flag == "SCEE_SCALE_FACTOR" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    SCEE_SCALE_FACTOR.push_back(value);
                }
            }
            else if(current_flag == "SCNB_SCALE_FACTOR" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    SCNB_SCALE_FACTOR.push_back(value);
                }
            }
            else if(current_flag == "SOLTY" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    SOLTY.push_back(value);
                }
            }
            else if(current_flag == "LENNARD_JONES_ACOEF" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    LENNARD_JONES_ACOEF.push_back(value);
                }
            }
            else if(current_flag == "LENNARD_JONES_BCOEF" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    LENNARD_JONES_BCOEF.push_back(value);
                }
            }

           
        }


    }
}