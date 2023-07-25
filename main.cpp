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

Environment PISUM_SATIVUM("chignolintemp.pdb");
parm test("chigparmtemp.parm7"); 

simulation small_probe(PISUM_SATIVUM, test, 1);

small_probe.force_additions();
small_probe.update_coord(0.0001, 100001);



std::cout<< "Completed!" << std::endl;

return 0;
}





