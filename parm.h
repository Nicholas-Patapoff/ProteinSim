#ifndef parm_H
#define parm_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class parm{
private:

public: 
parm(std::string pdb);
std::vector<std::string> ATOM_NAME;
std::vector<float> CHARGE;
};
#endif