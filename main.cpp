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
std::cout<< "Completed!" << std::endl;

simulation small_probe(PISUM_SATIVUM, test, 1);

small_probe.force_additions();

for(int i = 0; i < 20; i ++){
    small_probe.exports(i);
    small_probe.update_coord(1);
}


return 0;
}
