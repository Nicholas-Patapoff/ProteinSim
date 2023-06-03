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

public:
std::vector<float> velocities;


simulation(Environment& coord, parm& top, float step);
void random_vel();
void update_coord(float step_size);
void exports(int count);

};


#endif