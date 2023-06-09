#ifndef Sim_H
#define Sim_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include "Env.h"
#include "parm.h"

class simulation{
private:
std::unique_ptr<Environment> coord;
std::unique_ptr<parm> top;
void displacement_vect(std::vector<float>& d, int atom1, int atom2);
//float spring_force;
void magnitude(std::vector<float>& object, float& mag);
void unit_vector(float& mag, std::vector<float> d, std::vector<float>& unitv);
using T = std::variant<int, float, std::string>;
public:
std::vector<float> velocities;
std::vector<float> forces;



simulation(Environment& coord, parm& top, float step);

void update_coord(float step_size, int frames);
void exports(int count);
void spring_force(int atom1, int atom2, float kval, float eq);
void force_additions();
void VerletAlg(float& step); 

};


#endif