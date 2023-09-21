#ifndef VECTOR_MATH_H
#define VECTOR_MATH_H

#include <vector>
#include <memory>

void vect_add(std::vector<float>& v1, std::vector<float>& v2, std::vector<float>& product);
void resize(std::vector<float>& vect, float scale);
void magnitude(std::vector<float>& object, float& mag);
void displacement_vect(std::vector<float>& d, std::vector<float>& Acoords, int atom1, int atom2);
void unit_vector(float& mag, std::vector<float> d, std::vector<float>& unitv);
void dot(std::vector<float>& disp1, std::vector<float>& disp2, float& val);
void cross(std::vector<float>& vect1, std::vector<float>& vect2, std::vector<float>& cprod);




#endif