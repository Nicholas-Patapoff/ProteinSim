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
void magnitude(std::vector<float>& object, float& mag);
void unit_vector(float& mag, std::vector<float> d, std::vector<float>& unitv);
void theta_from_dot(int& atom1, int& atom2, int& atom3, float& theta);
void dot(std::vector<float>& disp1, std::vector<float>& disp2, float& val);
void cross(std::vector<float>& vect1, std::vector<float>& vect2, std::vector<float>& cprod);
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
void angle_force(int atom1, int atom2, int atom3, float k, float eq);

};


#endif