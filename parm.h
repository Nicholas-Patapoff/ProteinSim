#ifndef parm_H
#define parm_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class parm{
private:

public: 
parm(std::string pdb);
std::vector<std::string> ATOM_NAME;
std::vector<float> CHARGE;
std::vector<int> ATOMIC_NUMBER;
std::vector<float> MASS;
std::vector<int> ATOM_TYPE_INDEX;
std::vector<int> NONBONDED_PARM_INDEX;
std::vector<std::string> RESIDUE_LABEL;
std::vector<int> RESIDUE_POINTER;
std::vector<float> BOND_FORCE_CONSTANT;
std::vector<float> BOND_EQUIL_VALUE;
std::vector<float> ANGLE_FORCE_CONSTANT;
std::vector<float> ANGLE_EQUIL_VALUE;
std::vector<float> DIHEDRAL_FORCE_CONSTANT;
std::vector<float> DIHEDRAL_PERIODICITY;
std::vector<float> DIHEDRAL_PHASE;
std::vector<float> SCEE_SCALE_FACTOR;
std::vector<float> SCNB_SCALE_FACTOR;
std::vector<float> SOLTY; //not currently implemented in AMBER
std::vector<float> LENNARD_JONES_ACOEF;
std::vector<float> LENNARD_JONES_BCOEF;


};
#endif