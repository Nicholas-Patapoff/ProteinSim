#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include "Env.h"
#include "parm.h"
#include "Sim.h"
#include "math.h"


simulation::simulation(Environment& temp1, parm& temp2, float step){

    coord = std::make_unique<Environment>(temp1);
    top = std::make_unique<parm>(temp2);



    velocities = std::vector<float>(coord->Acoords.size(), 0);
}


void simulation::update_coord(float step_size){

for(int i = 0; i < velocities.size(); i++){
    coord->Acoords[i] +=  velocities[i] * step_size;
}
//update git


}

std::vector<float> simulation::spring_force(int atom1, int atom2, float k, float eq){

    std::vector<float> disp;
    float dmag;
    std::vector<float> dunit_vect;
    std::vector<float> forces;
    displacement_vect(disp, atom1, atom2);
    magnitude(disp, dmag);
    unit_vector(dmag, disp, dunit_vect);
    for(int i = 0; i < disp.size(); i ++){
        forces.push_back(-k * (disp[i] - eq * dunit_vect[i]));
    }
    return forces;
}

void simulation::unit_vector(float& mag, std::vector<float> d, std::vector<float>& unitv){
    for(int i = 0; i < d.size(); i++){
        unitv.push_back(d[i]/mag);
    }
}

void simulation::magnitude(std::vector<float>& object, float& mag){
    float temp = 0;
    for(int i = 0; i < object.size(); i++){
        temp+= object[i]*object[i];
    }
    temp = std::sqrt(temp);
    mag = temp;
}

void simulation::displacement_vect(std::vector<float>& d, int atom1, int atom2){
    for(int i = 0; i < 3; i++){
        d.push_back(coord->Acoords[atom1 * 3 + i] - coord->Acoords[atom2 * 3 + i]);
    }
}

void simulation::exports(int count){
    std::fstream temp_file;
    temp_file.open("coord_data/" + std::to_string(count), std::ios::out);
    for(int i = 0; i < coord->Acoords.size() - 1; i+=3){
        temp_file << std::to_string(coord->Acoords[i]) << " " << std::to_string(coord->Acoords[i+1]) << " " << std::to_string(coord->Acoords[i+2]) << "\n";
    }
}