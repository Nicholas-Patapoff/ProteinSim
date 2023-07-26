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

simulation::simulation(Environment& temp1, parm& temp2, float step, std::string export_name){

    coord = std::make_unique<Environment>(temp1);
    top = std::make_unique<parm>(temp2);



    velocities = std::vector<float>(coord->Acoords.size(), 0);
    forces = std::vector<float>(coord->Acoords.size(), 0);
    


    temp_file.open("/home/nich/Documents/GitHub/coord_vis/" + export_name +".crd", std::ios::out);
    temp_file << std::endl;
    }


void simulation::update_coord(float step_size, int frames, int export_step){
    for(int i = 0; i <= frames; i++){
        VerletAlg(step_size);
        
        if(i % export_step == 0) {
            std::cout<< (float) i / frames * 100 << "% complete!" << std::endl;
            float total_f = 0;
            for(int i = 0; i < forces.size(); i++){
                total_f += forces[i];
            }
            exports();
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

    //Dihedral forces
    std::unordered_map<std::string, std::vector<T> >::iterator dih = top->values.find("DIHEDRALS_INC_HYDROGEN");
        std::vector<T>& DincH = dih->second;

    std::unordered_map<std::string, std::vector<T> >::iterator dwh = top->values.find("DIHEDRALS_WITHOUT_HYDROGEN");
        std::vector<T>& DWoutH = dwh->second;

    std::unordered_map<std::string, std::vector<T> >::iterator dfc = top->values.find("DIHEDRAL_FORCE_CONSTANT");
        std::vector<T>& DForceC = dfc->second;

    std::unordered_map<std::string, std::vector<T> >::iterator dp = top->values.find("DIHEDRAL_PERIODICITY");
        std::vector<T>& DPeriod = dp->second;

    std::unordered_map<std::string, std::vector<T> >::iterator dph = top->values.find("DIHEDRAL_PHASE");
        std::vector<T>& DPhase = dph->second;

    std::unordered_map<std::string, std::vector<T> >::iterator scee = top->values.find("SCEE_SCALE_FACTOR");
        std::vector<T>& SCEE_SF = scee->second;      

    std::unordered_map<std::string, std::vector<T> >::iterator scnb = top->values.find("SCNB_SCALE_FACTOR");
        std::vector<T>& SCNB_SF = scnb->second;  

    for(int i = 0; i < DWoutH.size(); i+=5){
        dihedral_force( std::get<int>(DWoutH[i]) / 3,std::get<int>(DWoutH[i + 1]) / 3, std::get<int>(DWoutH[i + 2]) / 3, std::get<int>(DWoutH[i + 3]) / 3, std::get<float>(DForceC[std::get<int>(DWoutH[i + 4]) - 1]), 
        std::get<float>(DPeriod[std::get<int>(DWoutH[i + 4]) - 1]), std::get<float>(SCEE_SF[std::get<int>(DWoutH[i + 4]) - 1]), std::get<float>(SCNB_SF[std::get<int>(DWoutH[i + 4]) - 1]), std::get<float>(DPhase[std::get<int>(DWoutH[i + 4]) - 1]));

    }
    for(int i = 0; i < DincH.size(); i+=5){
        dihedral_force( std::get<int>(DincH[i]) / 3,std::get<int>(DincH[i + 1]) / 3, std::get<int>(DincH[i + 2]) / 3, std::get<int>(DincH[i + 3]) / 3, std::get<float>(DForceC[std::get<int>(DincH[i + 4]) - 1]), 
        std::get<float>(DPeriod[std::get<int>(DincH[i + 4]) - 1]), std::get<float>(SCEE_SF[std::get<int>(DincH[i + 4]) - 1]), std::get<float>(SCNB_SF[std::get<int>(DincH[i + 4]) - 1]), std::get<float>(DPhase[std::get<int>(DincH[i + 4]) - 1]));

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

    std::vector<float> orth_abc, inp_ba, inp_bc;
    cross(disp1, disp2, orth_abc);
    cross(disp1, orth_abc, inp_ba);
    cross(disp2, orth_abc, inp_bc);

    float magba, magbc, mag_inp_ba, mag_inp_bc;
    magnitude(disp1, magba);
    magnitude(disp2, magbc);
    magnitude(inp_ba, mag_inp_ba);
    magnitude(inp_bc, mag_inp_bc);

    std::vector<float> unit_inp_ba, unit_inp_bc;
    unit_vector(mag_inp_ba, inp_ba, unit_inp_ba);
    unit_vector(mag_inp_bc, inp_bc, unit_inp_bc);

    float theta; 
    theta_from_dot(atom1, atom2, atom3, theta);

    float force_ba, force_bc;
    for(int i = 0; i < unit_inp_ba.size(); i++){
        force_ba =  k * (eq - theta)/magba * unit_inp_ba[i];
        force_bc =  -k * (eq - theta)/magbc * unit_inp_bc[i];
        forces[atom1 * 3 + i] += force_ba;
        forces[atom3 * 3 + i] += force_bc;
        forces[atom2 * 3 + i] += -force_ba - force_bc;
        
    }    

}

void simulation::dihedral_force(int atom1, int atom2, int atom3, int atom4, float k, float period, float sceef, float scnbf, float phase){
//torsion potential is defined as U= 0.5[A1(1 + cos(theta)) + A2(1 - cos2theta)) + A3(1 + cost(3theta)) + A4]
//the partial derivative of tors w/ respect to pos of atom a (pa) is: dU/dpa = du/dtheta * dtheta/dpa
//dU/dpa = +/- 0.5(-A1sin(theta) + 2A2sing(2theta) - 3A3sin(3theta))
//dtheta/dpa = 1/(mag(end_vector) * sin(theta))
//steps to find forces on end atoms: 1.create unit vectors whose direction is orthogonal abc and bcd 2. find the angle between the planes 
//3. calculate the forces for atoms a and d usign (0.5/mag(ab)sin(theta)) * (1A1sin(theta) - 2A2sin(2theta) + 3A3(sin(theta)) * norm vector orth to plane. 

//1.
std::vector<float> dispba, dispbc, dispcb, dispcd;
displacement_vect(dispba, atom2, atom1);
displacement_vect(dispbc, atom2, atom3);
displacement_vect(dispcb, atom3, atom2);
displacement_vect(dispcd, atom3, atom4);

std::vector<float> orthabc, orthbcd;
cross(dispba, dispbc, orthabc);
cross(dispcd, dispcb, orthbcd);

float orthabcmag, orthbcdmag, normabcmag, normbcdmag, magab, magcd;
magnitude(orthabc, orthabcmag);
magnitude(orthbcd, orthbcdmag);
magnitude(dispba, magab);
magnitude(dispcd, magcd);

std::vector<float> normabc, normbcd;
unit_vector(orthabcmag, orthabc, normabc);
unit_vector(orthbcdmag, orthbcd, normbcd);
magnitude(normabc,normabcmag);
magnitude(normbcd,normbcdmag);

//2
float theta;
DHtheta_from_dot(normabc, normbcd, normabcmag, normbcdmag, theta);

//3
std::vector<float> fba_direction, fcd_direction, nfba_direction, nfcd_direction;
cross(dispba, normabc, fba_direction);
cross(dispcd, normbcd, fcd_direction);

float forcea, forced, mag_badirection, mag_cddirection;
magnitude(fba_direction, mag_badirection);
magnitude(fcd_direction, mag_cddirection);
unit_vector(mag_badirection, fba_direction, nfba_direction);
unit_vector(mag_cddirection, fcd_direction, nfcd_direction);

forcea = (-k * period * sin(period * theta + phase))/ magab;//update for new mag as shown in paper. will be distance to origin. 
forced = (-k * period * sin(period * theta + phase))/ magcd;


for(int i = 0; i < 3; i++){
        forces[atom1 * 3 + i] +=  forcea * normabc[i];
        forces[atom4 * 3 + i] +=  forced * normbcd[i];
        std::cout << atom1 << " " <<forces[atom1 * 3 + i] << " ";
        std::cout << forces[atom4 * 3 + i] << " ";
    }    
    std::cout << std::endl;
}




void simulation::DHtheta_from_dot(std::vector<float>& nplane1, std::vector<float>& nplane2, float np1mag, float np2mag, float& theta){
    float dot12 = 0;
    dot(nplane1, nplane2, dot12);
    float cosine_value = dot12 / (np1mag * np2mag);
    cosine_value = std::max(-1.0f, std::min(1.0f, cosine_value));
    theta = acos(cosine_value);
    std::cout << " theta then dot: " << dot12 << " " << np1mag << " " << np2mag << " " << theta << std::endl;
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
   velocities[atom] = velocities[atom] + forces[atom]*step_size / (2 * std::get<float>(Mass[atom/(int)3]) );
   
}
//std::cout << coord->Acoords[0] << " " << coord->Acoords[1] << " " << coord->Acoords[2] << std::endl;

}

void simulation::exports(){
    if(!temp_file.is_open()){
        std::cout << "failed open external file" << std::endl;
    }
    for(int i = 0; i < coord->Acoords.size(); i++){
        if(i % 10 == 0 && i != 0){
            temp_file << std::endl;
        }
        temp_file << std::fixed << std::right << std::setw(8) << std::setprecision(3) << (coord->Acoords[i]);
        if(i == coord->Acoords.size() - 1){
            temp_file << std::endl;
            for(int z = 1; z < 4; z++){
                temp_file << std::fixed << std::right << std::setw(8) << std::setprecision(3) << (1.000);
            }
            
        }
        
    }
    temp_file << std::endl;

}



