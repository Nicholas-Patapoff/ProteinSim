#include <iostream>
#include <fstream> 
#include <string>
#include <variant>
#include <unordered_map>
#include "Env.h"
#include "parm.h"
#include "Sim.h"
int main()
{

Environment PISUM_SATIVUM("GLY2.pdb");
parm test("GLY.parm7"); 

simulation small_probe(PISUM_SATIVUM, test, 1, "GLY");

small_probe.force_additions();
small_probe.update_coord(0.01, 50000, 4);



std::cout<< "Completed!" << std::endl;

return 0;
}





