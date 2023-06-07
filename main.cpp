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
std::vector<float> item = small_probe.spring_force(0, 1, 1, 1);
for(int i = 0; i < item.size(); i++){
    std::cout<< item[i] << std::endl;
}

small_probe.random_vel();
for(int i = 0; i < 20; i ++){
    small_probe.exports(i);
    small_probe.update_coord(1);
}


return 0;
}
