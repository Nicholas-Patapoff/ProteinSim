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

simulation small_probe(PISUM_SATIVUM, test, 1, "GLYDIH");

small_probe.force_additions();
small_probe.update_coord(0.01, 2000, 20);



std::cout<< "Completed!" << std::endl;

return 0;
}





