#ifndef parm_H
#define parm_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <variant>
#include <unordered_map>

class parm{
private:

public: 
using T = std::variant<int, float, std::string>;
std::unordered_map<std::string, std::vector<T> > values;
//std::unordered_map<std::string, std::vector<T> >::iterator test;
parm(std::string pdb);





/*
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
std::vector<int> BONDS_WITHOUT_HYDROGEN;
std::vector<int> ANGLES_INC_HYDROGEN;
std::vector<int> DIHEDRALS_INC_HYDROGEN;
std::vector<int> DIHEDRALS_WITHOUT_HYDROGEN;
std::vector<int> EXCLUDED_ATOMS_LIST;
std::vector<float> HBOND_ACOEF;
std::vector<float> HBOND_BCOEF;
std::vector<float> HBCUT;
std::vector<std::string> AMBER_ATOM_TYPE;
std::vector<std::string> TREE_CHAIN_CLASSIFICATION;
std::vector<int> JOIN_ARRAY;
std::vector<int> IROTAT;
std::vector<int> SOLVENT_POINTERS;
std::vector<int> ATOMS_PER_MOLECULE;
std::vector<float> BOX_DIMENSIONS;
std::vector<int> CAP_INFO;
std::vector<float> CAP_INTO2;
std::vector<std::string> RADIUS_SET;
std::vector<float> RADII;
std::vector<int> IPOL;
std::vector<float> POLARIZABILITY;
*/
};


#endif