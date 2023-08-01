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

Environment PISUM_SATIVUM("chignolin_ext.pdb");
parm test("chig_ext.parm7"); 

simulation small_probe(PISUM_SATIVUM, test, 1, "chig");

small_probe.force_additions();
small_probe.update_coord(0.01, 1000, 2);



std::cout<< "Completed!" << std::endl;

return 0;
}





