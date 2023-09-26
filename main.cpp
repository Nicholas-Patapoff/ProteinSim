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

Environment coords("test_prot_ext.pdb");
parm data("test_prot.parm7"); 

simulation small_probe(coords, data, 1, "test_prot");

small_probe.update_coord(0.01, 20000, 5);



std::cout<< "Completed!" << std::endl;

return 0;
}





