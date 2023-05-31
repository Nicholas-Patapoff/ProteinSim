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
//for(int i = 0; i < test.ATOM_NAME.size(); i++){
//std::cout << test.ATOM_NAME[i] << "/";
//}



return 0;
}
