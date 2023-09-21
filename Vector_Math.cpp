#include "Vector_Math.h"
#include "Env.h"
#include <math.h>

void vect_add(std::vector<float>& v1, std::vector<float>& v2, std::vector<float>& product){

    for(int i = 0; i < v1.size(); i++){
        product[i] = v1[i] + v2[i]; 
    }

}

void resize(std::vector<float>& vect, float scale){
    for(int i = 0; i < vect.size(); i++){
        vect[i] *= scale;
    }
}

void magnitude(std::vector<float>& object, float& mag){
    float temp = 0;
    for(int i = 0; i < object.size(); i++){
        temp+= object[i] *object[i];
    }
    temp = std::sqrt(temp);
    mag = temp;
}

void displacement_vect(std::vector<float>& d, std::vector<float>& Acoords , int atom1, int atom2){
    for(int i = 0; i < 3; i++){
        d.push_back(Acoords[atom1 * 3 + i] - Acoords[atom2 * 3 + i]);
    }
}

void unit_vector(float& mag, std::vector<float> d, std::vector<float>& unitv){
    
    for(int i = 0; i < d.size(); i++){
        unitv.push_back(d[i]/mag);
    }
}

void dot(std::vector<float>& disp1, std::vector<float>& disp2, float& val){
    for(int i = 0; i < 3; i++){
        val += disp1[i] * disp2[i];
    }
}

void cross(std::vector<float>& vect1, std::vector<float>& vect2, std::vector<float>& cprod){

    cprod.push_back(vect1[1] * vect2[2] - vect1[2] * vect2[1]);
    cprod.push_back(-(vect1[0] * vect2[2] - vect1[2] * vect2[0]));
    cprod.push_back(vect1[0] * vect2[1] - vect1[1] * vect2[0]);

}
