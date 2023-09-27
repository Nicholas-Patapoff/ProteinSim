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
using T = std::variant<int, double, std::string>;
std::unique_ptr<Environment> coord;
std::unique_ptr<parm> top;
std::fstream temp_file;

void theta_from_dot(int& atom1, int& atom2, int& atom3, double& theta);
void DHtheta_from_dot(std::vector<double>& nplane1, std::vector<double>& nplane2, double np1mag, double np2mag, double& theta);
void DHrotTheta_from_dot(std::vector<double>& disp1, std::vector<double>& disp2, double& mag_disp1, double& mag_disp2, double& theta);
void DH_LJF(int atom1, int atom4, double SCNBF, double LJA, double LJB);

public:
std::vector<double> velocities;
std::vector<double> forces;



simulation(Environment& coord, parm& top, double step, std::string export_name);
void update_coord(double step_size, int frames, int export_step);
void exports();
void spring_force(int atom1, int atom2, double kval, double eq);
void force_additions(std::vector<T>& BWoutH, std::vector<T>& BIH, std::vector<T>& BForceC, std::vector<T>& BEQV, std::vector<T>& AWoutH,
 std::vector<T>& AIH,std::vector<T>& AForceC, std::vector<T>& AEQV, std::vector<T>& DincH, std::vector<T>& DWoutH, std::vector<T>& DForceC, 
 std::vector<T>& DPeriod , std::vector<T>& DPhase, std::vector<T>& SCEE_SF, std::vector<T>& SCNB_SF, std::vector<T>& LJAC, std::vector<T>& LJBC,std::vector<T>& ATI,
  std::vector<T>& NBPIndex,std::vector<std::vector<int> >& excluded);
void VerletAlg(double& step, std::vector<T>& BWoutH, std::vector<T>& BIH, std::vector<T>& BForceC, std::vector<T>& BEQV, std::vector<T>& AWoutH,
 std::vector<T>& AIH,std::vector<T>& AForceC, std::vector<T>& AEQV, std::vector<T>& DincH, std::vector<T>& DWoutH, std::vector<T>& DForceC, 
 std::vector<T>& DPeriod , std::vector<T>& DPhase, std::vector<T>& SCEE_SF, std::vector<T>& SCNB_SF, std::vector<T>& LJAC, std::vector<T>& LJBC,std::vector<T>& ATI,
  std::vector<T>& NBPIndex, std::vector<std::vector<int> >& excluded); 
void angle_force(int atom1, int atom2, int atom3, double k, double eq);
void dihedral_force(int atom1, int atom2, int atom3, int atom4, double k, double period, double sceef, double scnbf, double phase);
void LennardJ_force(int atom1, int atom2, double LJA, double LJB);
void improper_dih_force(int atom1, int atom2, int atom3, int atom4, double k, double period, double sceef, double scnbf, double phase);
};


#endif