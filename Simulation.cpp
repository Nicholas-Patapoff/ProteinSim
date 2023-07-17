#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include "Env.h"
#include "parm.h"
#include "Sim.h"
#include <math.h>
#include <thread>
#include <iomanip>

simulation::simulation(Environment& temp1, parm& temp2, float step){

    coord = std::make_unique<Environment>(temp1);
    top = std::make_unique<parm>(temp2);



    velocities = std::vector<float>(coord->Acoords.size(), 0);
    forces = std::vector<float>(coord->Acoords.size(), 0);
    
    }


void simulation::update_coord(float step_size, int frames){
    int counting = 0;
    for(int i = 0; i < frames; i++){
        VerletAlg(step_size);
        
        if(i % 50 == 0) {
            std::cout << i << std::endl;
            std::cout << "vel: " << velocities[6] << " " << velocities[7] << " " << velocities[8] << std::endl;
            std::cout << "vel: " << velocities[18] << " " << velocities[19] << " " << velocities[20] << std::endl;
            std::cout << "vel: " << velocities[21] << " " << velocities[22] << " " << velocities[23] << std::endl;
            std::cout << "forces6: " << forces[6] << " " << forces[7] << " " << forces[8] << std::endl;
            std::cout << "forces18: " << forces[18] << " " << forces[19] << " " << forces[20] << std::endl;
            std::cout << "forces21: " << forces[21] << " " << forces[22] << " " << forces[23] << std::endl;
            float total_f = 0;
            for(int i = 0; i < forces.size(); i++){
                total_f += forces[i];
            }
            if(total_f != 0){
                //std::cout << "not all add to 0" << std::endl;
                std::cout << total_f << std::endl;
            }
                    std::this_thread::sleep_for(std::chrono::milliseconds(0));

            //exports(counting);
            counting += 1;

        }
    }
}

void simulation::force_additions(){ //all forces on protien for bond/FF interactions

    //bond forces
    std::unordered_map<std::string, std::vector<T> >::iterator bwh = top->values.find("BONDS_WITHOUT_HYDROGEN");
        std::vector<T>& BWoutH = bwh->second;
    
    std::unordered_map<std::string, std::vector<T> >::iterator bih = top->values.find("BONDS_INC_HYDROGEN");
        std::vector<T>& BIH = bih->second;

    std::unordered_map<std::string, std::vector<T> >::iterator bfc = top->values.find("BOND_FORCE_CONSTANT");
        std::vector<T>& BForceC = bfc->second;

    std::unordered_map<std::string, std::vector<T> >::iterator bev = top->values.find("BOND_EQUIL_VALUE");
        std::vector<T>& BEQV = bev->second;
    for(int i = 0; i < BWoutH.size() - 1; i+=3){ //BWoutH.size()
        spring_force( std::get<int>(BWoutH[i]) / 3,std::get<int>(BWoutH[i + 1]) / 3, std::get<float>(BForceC[std::get<int>(BWoutH[i + 2]) - 1]), std::get<float>(BEQV[std::get<int>(BWoutH[i + 2]) - 1]) );

    }
    for(int i = 0; i < BIH.size(); i+=3){ // this is a bond which include hydrogen. I want to stiffen interaction so as to reduce computational complexity. I can set up a method for selecting a model and integrate it into there
        spring_force( std::get<int>(BIH[i]) / 3,std::get<int>(BIH[i + 1]) / 3, std::get<float>(BForceC[std::get<int>(BIH[i + 2]) - 1]), std::get<float>(BEQV[std::get<int>(BIH[i + 2]) - 1]) );
    }


    //angle forces
    std::unordered_map<std::string, std::vector<T> >::iterator awh = top->values.find("ANGLES_WITHOUT_HYDROGEN");
        std::vector<T>& AWoutH = awh->second;

    std::unordered_map<std::string, std::vector<T> >::iterator aih = top->values.find("ANGLES_INC_HYDROGEN");
        std::vector<T>& AIH = aih->second;

    std::unordered_map<std::string, std::vector<T> >::iterator afc = top->values.find("ANGLE_FORCE_CONSTANT");
        std::vector<T>& AForceC = afc->second;

    std::unordered_map<std::string, std::vector<T> >::iterator aev = top->values.find("ANGLE_EQUIL_VALUE");
        std::vector<T>& AEQV = aev->second;
    for(int i = 0; i < AWoutH.size(); i+=4){
        angle_force( std::get<int>(AWoutH[i]) / 3,std::get<int>(AWoutH[i + 1]) / 3, std::get<int>(AWoutH[i + 2]) / 3, std::get<float>(AForceC[std::get<int>(AWoutH[i + 3]) - 1]), std::get<float>(AEQV[std::get<int>(AWoutH[i + 3]) - 1]));

    }

    for(int i = 0; i < AIH.size() - 1; i+=4){ 
        angle_force( std::get<int>(AIH[i]) / 3,std::get<int>(AIH[i + 1]) / 3, std::get<int>(AIH[i + 2]) / 3, std::get<float>(AForceC[std::get<int>(AIH[i + 3]) - 1]), std::get<float>(AEQV[std::get<int>(AIH[i + 3]) - 1]));

    }

    
}

void simulation::angle_force(int atom1, int atom2, int atom3, float k, float eq){ //this needs to be optimized so as not to make disp and mags two times
    //first find the change of potential with respect the change of angle : 2k(thetaEQ - current_theta)
    //then find the change of potential with respect to the change of bond length
    // we must then make a nomalized vector which is in plane to atoms 1,2,3 and points in a direction orth to vector 1-2, and then 2-3
    // from there we can apply forces via F12 = - du/dr dot orth_vect12 which: -2k(thetaEQ - current_theta)/mag(vector1-2) dot orth_vect12.
    //F23 is the same 
    //then take the total of both as a negative ( - f12 - f23) and then we obtain force on atom 2. since -FA + -FC = FB;
    
    std::vector<float> disp1, disp2, dispac;
    displacement_vect(disp1, atom1, atom2);
    displacement_vect(disp2, atom3, atom2);
    displacement_vect(dispac, atom1, atom3);
    
    float magba, magbc;
    magnitude(disp1, magba);
    magnitude(disp2, magbc);

    std::vector<float> unitv_ba, unitv_bc;
    unit_vector(magba, disp1, unitv_ba);
    unit_vector(magbc, disp2, unitv_bc);

    float theta; 
    theta_from_dot(atom1, atom2, atom3, theta);
    std::cout<< theta << std::endl;

    std::vector<float> force_ba, force_bc;



    
}


void simulation::theta_from_dot(int& atom1, int& atom2, int& atom3, float& theta){
    std::vector<float> disp1, disp2;
    displacement_vect(disp1, atom1, atom2);
    displacement_vect(disp2, atom3, atom2);
    float magba;
    magnitude(disp1, magba);
    float magbc;
    magnitude(disp2, magbc);
    float dotac = 0;
    dot(disp1, disp2, dotac);
    std::vector<float> ba_unitvect, bc_unitvect;


    theta = acos(dotac/(magba*magbc));


}
void simulation::cross(std::vector<float>& vect1, std::vector<float>& vect2, std::vector<float>& cprod){

    cprod.push_back(vect1[1] * vect2[2] - vect1[2] * vect2[1]);
    cprod.push_back(-(vect1[0] * vect2[2] - vect1[2] * vect2[0]));
    cprod.push_back(vect1[0] * vect2[1] - vect1[1] * vect2[0]);

}

void simulation::dot(std::vector<float>& disp1, std::vector<float>& disp2, float& val){
    for(int i = 0; i < 3; i++){
        val += disp1[i] * disp2[i];
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
        float force = -0.5 * (-k * (disp[i] - eq * dunit_vect[i]));
        forces[atom2 * 3 + i] += force;
        forces[atom1 * 3 + i] -= force;

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
   velocities[atom] = velocities[atom] + forces[atom]*step_size / (2 * std::get<float>(Mass[atom/(int)3]) );//6.02214086e+26
}


}

void simulation::exports(int count){
    std::fstream temp_file;
    temp_file.open("temp.crd", std::ios::out);
    if(!temp_file.is_open()){
        std::cout << "failed" << std::endl;
    }
    temp_file << std::left << std::setw(20) << "title" << "\n";
    std::cout << count << std::endl;
    for(int i = 1; i < coord->Acoords.size(); i++){

        temp_file << std::fixed << std::right << std::setw(8) << std::setprecision(3) << (coord->Acoords[i - 1]);
        if(i % 9 == 0){
            temp_file << std::endl;
        }
            }
}



