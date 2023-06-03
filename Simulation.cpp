#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include "Env.h"
#include "parm.h"
#include "Sim.h"


simulation::simulation(Environment& temp1, parm& temp2, float step){

coord = std::make_unique<Environment>(temp1);
top = std::make_unique<parm>(temp2);

velocities = std::vector<float>(coord->Acoords.size(), 0);


}

void simulation::random_vel(){

for(int i = 0; i < velocities.size(); i++){
    velocities[i] = rand() % 7 - 3; 
}

}

void simulation::update_coord(float step_size){

for(int i = 0; i < velocities.size(); i++){
    coord->Acoords[i] +=  velocities[i] * step_size;
}
//update git


}

void simulation::exports(int count){
    std::fstream temp_file;
    temp_file.open("coord_data/" + std::to_string(count), std::ios::out);
    for(int i = 0; i < coord->Acoords.size() - 1; i+=3){
        temp_file << std::to_string(coord->Acoords[i]) << " " << std::to_string(coord->Acoords[i+1]) << " " << std::to_string(coord->Acoords[i+2]) << "\n";
    }

}