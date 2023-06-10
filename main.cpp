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

Environment PISUM_SATIVUM("output.pdb");
parm test("output.parm7"); 

simulation small_probe(PISUM_SATIVUM, test, 1);

small_probe.force_additions();
small_probe.update_coord(0.000001, 10000000);



std::cout<< "Completed!" << std::endl;

return 0;
}
