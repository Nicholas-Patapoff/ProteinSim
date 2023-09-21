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
using T = std::variant<int, float, std::string>;
std::unique_ptr<Environment> coord;
std::unique_ptr<parm> top;
std::fstream temp_file;

void theta_from_dot(int& atom1, int& atom2, int& atom3, float& theta);
void DHtheta_from_dot(std::vector<float>& nplane1, std::vector<float>& nplane2, float np1mag, float np2mag, float& theta);
void DHrotTheta_from_dot(std::vector<float>& disp1, std::vector<float>& disp2, float& mag_disp1, float& mag_disp2, float& theta);
void DH_LJF(int atom1, int atom4, float SCNBF, float LJA, float LJB);

public:
std::vector<float> velocities;
std::vector<float> forces;



simulation(Environment& coord, parm& top, float step, std::string export_name);
void update_coord(float step_size, int frames, int export_step);
void exports();
void spring_force(int atom1, int atom2, float kval, float eq);
void force_additions(std::vector<T>& BWoutH, std::vector<T>& BIH, std::vector<T>& BForceC, std::vector<T>& BEQV, std::vector<T>& AWoutH,
 std::vector<T>& AIH,std::vector<T>& AForceC, std::vector<T>& AEQV, std::vector<T>& DincH, std::vector<T>& DWoutH, std::vector<T>& DForceC, 
 std::vector<T>& DPeriod , std::vector<T>& DPhase, std::vector<T>& SCEE_SF, std::vector<T>& SCNB_SF, std::vector<T>& LJAC, std::vector<T>& LJBC,std::vector<T>& ATI,
  std::vector<T>& NBPIndex,std::vector<std::vector<int> >& excluded);
void VerletAlg(float& step, std::vector<T>& BWoutH, std::vector<T>& BIH, std::vector<T>& BForceC, std::vector<T>& BEQV, std::vector<T>& AWoutH,
 std::vector<T>& AIH,std::vector<T>& AForceC, std::vector<T>& AEQV, std::vector<T>& DincH, std::vector<T>& DWoutH, std::vector<T>& DForceC, 
 std::vector<T>& DPeriod , std::vector<T>& DPhase, std::vector<T>& SCEE_SF, std::vector<T>& SCNB_SF, std::vector<T>& LJAC, std::vector<T>& LJBC,std::vector<T>& ATI,
  std::vector<T>& NBPIndex, std::vector<std::vector<int> >& excluded); 
void angle_force(int atom1, int atom2, int atom3, float k, float eq);
void dihedral_force(int atom1, int atom2, int atom3, int atom4, float k, float period, float sceef, float scnbf, float phase);
void LennardJ_force(int atom1, int atom2, float LJA, float LJB);
};


#endif