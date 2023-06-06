#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <variant>
#include "parm.h"

parm::parm(std::string parm7file){
    std::fstream nfile;
    nfile.open(parm7file, std::ios::in);
    if(nfile.is_open()){
        std::string temp;
        std::string current_flag;
        std::string current_format;
        while(getline(nfile, temp)){
            if(temp.substr(0,5) == "%FLAG"){ // this is good but problematic for files which have comments between the flags/formats and thier identifiers
                std::string format = temp;
                std::istringstream splitstring(temp);
                splitstring >> current_flag; // saves as "%flag"
                splitstring >> current_flag; // saves as "flag_type"
                getline(nfile, temp);
                std::cout << temp << std::endl;
                temp = temp.substr(7, 10);
                std::istringstream splitstring2(temp);
                splitstring2 >> current_format; //this saves as "%format"
                std::cout << temp << std::endl;

                if(current_format == "(20a4)" || current_format == "(1a80)"){
                    std::cout << current_format << std::endl;
                    values.emplace(current_flag, std::vector<T>());    
                }
                else if(current_format == "(5E16.8)"){
                    std::cout << current_format << std::endl;

                    values.emplace(current_flag, std::vector<T>());

                }
                else if(current_format == "(10I8)" || current_format == "(1I8)"){
                    std::cout << current_format << std::endl;

                    values.emplace(current_flag, std::vector<T>());
                }
            } 
            
            else if(current_format == "(20a4)"){
                std::unordered_map<std::string, std::vector<T> >::iterator test = values.find(current_flag);
                std::vector<T>& vec = test->second;

                for(int i = 0; i < temp.size() - 1; i+= 4){
                    vec.push_back(temp.substr(i, 4));
                }
                for (const auto& element : vec) {
                    if (std::holds_alternative<int>(element)) {
                        std::cout << "Integer value: " << std::get<int>(element) << std::endl;
                    } else if (std::holds_alternative<float>(element)) {
                        std::cout << "Float value: " << std::get<float>(element) << std::endl;
                    } else if (std::holds_alternative<std::string>(element)) {
                        std::cout << "String value: " << std::get<std::string>(element) << std::endl;
                    }
            }
            }

            else if(current_format == "(5E16.8)"){
                for(int i = 0; i < temp.size() - 1; i+= 16){
                    //values[current_flag].push_back(std::stof(temp.substr(i, 16)));
                }
            }

            else if(current_format == "(10I8)"){
                for(int i = 0; i < temp.size() - 1; i+= 8){
                    //this->values[current_flag].push_back(std::stoi(temp.substr(i, 8)));
                }

            }

            else if(current_format == "(1a80)"){
                //values[current_flag].push_back(temp);
            }

            else if(current_format == "(1I8)"){
                //values[current_flag].push_back(std::stoi(temp));
            }






/*


                
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
            else if(current_flag == "BONDS_WITHOUT_HYDROGEN" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                int value;
                std::istringstream line(temp);
                while(line >> value){
                    BONDS_WITHOUT_HYDROGEN.push_back(value);
                }
            }
            else if(current_flag == "ANGLES_INC_HYDROGEN" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                int value;
                std::istringstream line(temp);
                while(line >> value){
                    ANGLES_INC_HYDROGEN.push_back(value);
                }
            }
            else if(current_flag == "DIHEDRALS_INC_HYDROGEN" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                int value;
                std::istringstream line(temp);
                while(line >> value){
                    DIHEDRALS_INC_HYDROGEN.push_back(value);
                }
            }
            else if(current_flag == "DIHEDRALS_WITHOUT_HYDROGEN" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                int value;
                std::istringstream line(temp);
                while(line >> value){
                    DIHEDRALS_WITHOUT_HYDROGEN.push_back(value);
                }
            }
            else if(current_flag == "EXCLUDED_ATOMS_LIST" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                int value;
                std::istringstream line(temp);
                while(line >> value){
                    EXCLUDED_ATOMS_LIST.push_back(value);
                }
            }
            else if(current_flag == "HBOND_ACOEF" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    HBOND_ACOEF.push_back(value);
                }
            }
            else if(current_flag == "HBOND_BCOEF" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    HBOND_BCOEF.push_back(value);
                }
            }
            else if(current_flag == "HBCUT" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    HBCUT.push_back(value);
                }
            }
            else if(current_flag == "AMBER_ATOM_TYPE" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                std::string value;
                std::istringstream line(temp);
                while(line >> value){
                    AMBER_ATOM_TYPE.push_back(value);
                }
            }
            else if(current_flag == "AMBER_ATOM_TYPE" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                std::string value;
                std::istringstream line(temp);
                while(line >> value){
                    AMBER_ATOM_TYPE.push_back(value);
                }
            }
            else if(current_flag == "TREE_CHAIN_CLASSIFICATION" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                std::string value;
                std::istringstream line(temp);
                while(line >> value){
                    TREE_CHAIN_CLASSIFICATION.push_back(value);
                }
            }
            else if(current_flag == "JOIN_ARRAY" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                int value;
                std::istringstream line(temp);
                while(line >> value){
                    JOIN_ARRAY.push_back(value);
                }
            }
            else if(current_flag == "IROTAT" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                int value;
                std::istringstream line(temp);
                while(line >> value){
                    IROTAT.push_back(value);
                }
            }
            else if(current_flag == "SOLVENT_POINTERS" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                int value;
                std::istringstream line(temp);
                while(line >> value){
                    SOLVENT_POINTERS.push_back(value);
                }
            }
           else if(current_flag == "ATOMS_PER_MOLECULE" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                int value;
                std::istringstream line(temp);
                while(line >> value){
                    ATOMS_PER_MOLECULE.push_back(value);
                }
            }
            else if(current_flag == "BOX_DIMENSIONS" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    BOX_DIMENSIONS.push_back(value);
                }
            }
            else if(current_flag == "CAP_INFO" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                int value;
                std::istringstream line(temp);
                while(line >> value){
                    CAP_INFO.push_back(value);
                }
            }
            else if(current_flag == "CAP_INTO2" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    CAP_INTO2.push_back(value);
                }
            }
            else if(current_flag == "CAP_INFO" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                RADIUS_SET.push_back(temp);
            }
            else if(current_flag == "RADII" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    RADII.push_back(value);
                }
            }
            else if(current_flag == "IPOL" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                int value;
                std::istringstream line(temp);
                while(line >> value){
                    IPOL.push_back(value);
                }
            }
            else if(current_flag == "POLARIZABILITY" && temp.substr(0,7) != "%FORMAT" && temp.substr(0,5) != "%FLAG"){
                float value;
                std::istringstream line(temp);
                while(line >> value){
                    POLARIZABILITY.push_back(value);
                }
            }
        }
*/

        } 
    }
}