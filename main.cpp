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

simulation small_probe(PISUM_SATIVUM, test, 1, "temp");

small_probe.force_additions();
small_probe.update_coord(0.0001, 20000, 2000);



std::cout<< "Completed!" << std::endl;

return 0;
}





