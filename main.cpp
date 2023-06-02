#include <iostream>
#include <fstream> 
#include <string>
#include "Env.h"
#include "parm.h"
int main()
{

Environment PISUM_SATIVUM("output.pdb");
parm test("output.parm7"); 
std::cout<< "Completed!" << std::endl;
std::cout << test.ATOM_NAME.size() << " ATOM_NAME"<< std::endl;
std::cout << test.CHARGE.size() << " CHARGE" << std::endl;
std::cout << test.ATOMIC_NUMBER.size() << " ATOMIC_NUMBER" << std::endl;
std::cout << test.MASS.size() << " MASS" << std::endl;
std::cout << test.ATOM_TYPE_INDEX.size() << " ATOM_TYPE_INDEX" << std::endl;
std::cout << test.NONBONDED_PARM_INDEX.size() << " NONBUNDED_PARM_INDEX" << std::endl;
std::cout << test.RESIDUE_LABEL.size() << " RESIDUE_LABEL" << std::endl;
std::cout << test.RESIDUE_POINTER.size() << " RESIDUE_POINTER" << std::endl;
std::cout << test.BOND_FORCE_CONSTANT.size() << " BOND_FORCE_CONSTANT" << std::endl;
std::cout << test.BOND_EQUIL_VALUE.size() << " BOND_EQUIL_VALUE" << std::endl;
std::cout << test.ANGLE_FORCE_CONSTANT.size() << " ANGLE_FORCE_CONSTANT" << std::endl;
std::cout << test.ANGLE_EQUIL_VALUE.size() << " ANGLE_EQUIL_VALUE" << std::endl; 
std::cout << test.DIHEDRAL_FORCE_CONSTANT.size() << " DIH_FORCE_CONST" << std::endl;
std::cout << test.DIHEDRAL_PERIODICITY.size() << " DIH_PERIODICITY" << std::endl;
std::cout << test.DIHEDRAL_PHASE.size() << " DIH_PHASE" << std::endl;
std::cout << test.SCEE_SCALE_FACTOR.size() << " SCEE_SCALE" << std::endl;
std::cout << test.SCNB_SCALE_FACTOR.size() << " SCNB_SCALE" << std::endl;
std::cout << test.SOLTY.size() << " SOLTY" << std::endl;
std::cout << test.LENNARD_JONES_ACOEF.size() << " LJAC" << std::endl;
std::cout << test.LENNARD_JONES_BCOEF.size() << " LJAB" << std::endl;
std::cout << test.BONDS_WITHOUT_HYDROGEN.size() << " BONDS_WITHOUT_H" << std::endl; 
std::cout << test.ANGLES_INC_HYDROGEN.size() << " ANGLES_INC_HYDROGEN" << std::endl; 
std::cout << test.DIHEDRALS_WITHOUT_HYDROGEN.size() << " DIHEDRALS_WITHOUT_HYDROGEN" << std::endl; 
std::cout << test.EXCLUDED_ATOMS_LIST.size() << " EXCLUDED_ATOMS_LIST" << std::endl; 
std::cout << test.HBOND_ACOEF.size() << " HBOND_ACOEF" << std::endl; 
std::cout << test.HBOND_BCOEF.size() << " HBOND_BCOEF" << std::endl; 
std::cout << test.AMBER_ATOM_TYPE.size() << " AMBER_ATOM_TYPE" << std::endl; 
std::cout << test.TREE_CHAIN_CLASSIFICATION.size() << " TREE_CHAIN_CLASSIFICATION" << std::endl; 
std::cout << test.JOIN_ARRAY.size() << " JOIN_ARRAY" << std::endl; 
std::cout << test.IROTAT.size() << " IROTAT" << std::endl; 
std::cout << test.SOLVENT_POINTERS.size() << " SOLVENT_POINTERS" << std::endl; 
std::cout << test.ATOMS_PER_MOLECULE.size() << " ATOMS_PER_MOLECULE" << std::endl; 
std::cout << test.BOX_DIMENSIONS.size() << " BOX_DIMENSIONS" << std::endl; 
std::cout << test.CAP_INFO.size() << " CAP_INFO" << std::endl; 
std::cout << test.CAP_INTO2.size() << " CAP_INTO2" << std::endl; 
std::cout << test.RADIUS_SET.size() << " RADIUS_SET" << std::endl; 
std::cout << test.RADII.size() << " RADII" << std::endl; 
std::cout << test.IPOL.size() << " IPOL" << std::endl; 
std::cout << test.POLARIZABILITY.size() << " POLARIZABILITY" << std::endl; 




return 0;
}
