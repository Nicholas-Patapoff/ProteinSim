#ifndef Sim_H
#define Sim_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include<memory>
#include "Env.h"
#include "parm.h"

class simulation{
private:
std::unique_ptr<Environment> coord;
std::unique_ptr<parm> top;
void displacement_vect(std::vector<float>& d, std::vector<float>& atom1, std::vector<float>& atom2);
//float spring_force;
void magnitude(std::vector<float>& object, float& mag);
void unit_vector(float& mag, std::vector<float> d, std::vector<float>& unitv);

public:
std::vector<float> velocities;


simulation(Environment& coord, parm& top, float step);
void random_vel();
void update_coord(float step_size);
void exports(int count);
std::vector<float> spring_force();

};


#endif