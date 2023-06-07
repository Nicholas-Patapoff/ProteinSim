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

std::vector<float> simulation::spring_force(){
std::vector<float> coorda = {0,0,0};
std::vector<float> coordb = {-1, -4, 3};
std::vector<float>& coord1 = coorda;
std::vector<float>& coord2 = coordb;
std::vector<float> d;
float dmag;
std::vector<float> dunit_vect;
std::vector<float> forces;
float ideal = 3;
float k = 1;
displacement_vect(d, coord1, coord2);
magnitude(d, dmag);
unit_vector(dmag, d, dunit_vect);
for(int i = 0; i < d.size(); i ++){
    forces.push_back(-k * (d[i] - ideal * dunit_vect[i]));
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
std::cout << temp << std::endl;
temp = std::sqrt(temp);
std::cout << temp << std::endl;
mag = temp;
}

void simulation::displacement_vect(std::vector<float>& d, std::vector<float>& atom1, std::vector<float>& atom2){
    for(int i = 0; i < atom1.size(); i++){
        d.push_back(atom2[i] - atom1[i]);
    }
    


}

void simulation::exports(int count){
    std::fstream temp_file;
    temp_file.open("coord_data/" + std::to_string(count), std::ios::out);
    for(int i = 0; i < coord->Acoords.size() - 1; i+=3){
        temp_file << std::to_string(coord->Acoords[i]) << " " << std::to_string(coord->Acoords[i+1]) << " " << std::to_string(coord->Acoords[i+2]) << "\n";
    }

}