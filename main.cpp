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
std::cout << test.ATOM_NAME.size() << std::endl;
std::cout << test.CHARGE.size() << std::endl;

//for(int i = 0; i < test.ATOM_NAME.size(); i++){
//std::cout << test.ATOM_NAME[i] << "/";
//}



return 0;
}
