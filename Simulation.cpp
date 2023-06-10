#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include "Env.h"
#include "parm.h"
#include "Sim.h"
#include <math.h>
#include <thread>

simulation::simulation(Environment& temp1, parm& temp2, float step){

    coord = std::make_unique<Environment>(temp1);
    top = std::make_unique<parm>(temp2);



    velocities = std::vector<float>(coord->Acoords.size(), 0);
    forces = std::vector<float>(coord->Acoords.size(), 0);
    
    }


void simulation::update_coord(float step_size, int frames){
    for(int i = 0; i < frames; i++){
        VerletAlg(step_size);
        if(i % 5000 == 0) {
            std::cout << i << std::endl;
            std::cout << "Forces: "<< forces[6] << " | " << forces[7] << " | " << forces[8] <<std::endl;
            std::cout << "coords: "<< coord->Acoords[6] << " | " << coord->Acoords[7] << " | " << coord->Acoords[8] <<std::endl;

        }
        //std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }

}

void simulation::force_additions(){ //all forces on protien for bond/FF interactions
    //find the value vector for all bonds that do not include hydrogen
    std::unordered_map<std::string, std::vector<T> >::iterator bwh = top->values.find("BONDS_WITHOUT_HYDROGEN");
        std::vector<T>& BWoutH = bwh->second;
    //find the value vector for all bonds that inolve atleast 1 hydrogen
    std::unordered_map<std::string, std::vector<T> >::iterator bih = top->values.find("BONDS_INC_HYDROGEN");
        std::vector<T>& BIH = bih->second;

    std::unordered_map<std::string, std::vector<T> >::iterator bfc = top->values.find("BOND_FORCE_CONSTANT");
        std::vector<T>& BForceC = bfc->second;

    std::unordered_map<std::string, std::vector<T> >::iterator bev = top->values.find("BOND_EQUIL_VALUE");
        std::vector<T>& BEQV = bev->second;
    for(int i = 0; i < BWoutH.size(); i+=3){
        spring_force( std::get<int>(BWoutH[i]) / 3,std::get<int>(BWoutH[i + 1]) / 3, std::get<float>(BForceC[std::get<int>(BWoutH[i + 2]) - 1]), std::get<float>(BEQV[std::get<int>(BWoutH[i + 2]) - 1]) );
    }
    for(int i = 0; i < BIH.size(); i+=3){ // this is a bond which include hydrogen. I want to stiffen interaction so as to reduce computational complexity. I can set up a method for selecting a model and integrate it into there
        //spring_force( std::get<int>(BIH[i]) / 3,std::get<int>(BIH[i + 1]) / 3, std::get<float>(BForceC[std::get<int>(BIH[i + 2]) - 1]), std::get<float>(BEQV[std::get<int>(BIH[i + 2]) - 1]) );
    }
}


void simulation::spring_force(int atom1, int atom2, float k, float eq){
    std::vector<float> disp;
    float dmag;
    std::vector<float> dunit_vect;

    displacement_vect(disp, atom1, atom2);
    magnitude(disp, dmag);
    unit_vector(dmag, disp, dunit_vect);

    for(int i = 0; i < disp.size(); i ++){
        forces[atom2 * 3 + i] += -0.5 * (-k * (disp[i] - eq * dunit_vect[i]));
        forces[atom1 * 3 + i] += -forces[atom2 * 3 + i];
    }

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

void simulation::VerletAlg(float& step_size){
    std::unordered_map<std::string, std::vector<T> >::iterator ms = top->values.find("MASS");
        std::vector<T>& Mass = ms->second;
for(int atom = 0; atom < velocities.size(); atom++){
    velocities[atom] = velocities[atom] + (forces[atom] * step_size/(2 * std::get<float>(Mass[atom/(int)3])));
}
for(int atom = 0; atom < coord->Acoords.size(); atom++){
    coord->Acoords[atom] = coord->Acoords[atom] + velocities[atom]*step_size;
}
forces.assign(forces.size(), 0);
force_additions();
for(int atom = 0; atom < velocities.size(); atom++){
   velocities[atom] = velocities[atom] + forces[atom]*step_size / (2 * std::get<float>(Mass[atom/(int)3])  );//6.02214086e+26
}


}

void simulation::exports(int count){
    std::fstream temp_file;
    temp_file.open("coord_data/" + std::to_string(count), std::ios::out);
    for(int i = 0; i < coord->Acoords.size() - 1; i+=3){
        temp_file << std::to_string(coord->Acoords[i]) << " " << std::to_string(coord->Acoords[i+1]) << " " << std::to_string(coord->Acoords[i+2]) << "\n";
    }
}