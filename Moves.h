#ifndef MOVES_H
#define MOVES_H

#include <vector>
#include <cmath>
#include <cstdint>  // For uint64_t

void setSpacer(double* arr, double dCore, int valency);

void rotateX(double* arr, double phi, int valency);

void rotateY(double* arr, double phi, int valency);

void rotateZ(double* arr, double phi, int valency);

double random(double factor, uint64_t *seed);  // Updated to use uint64_t for 64-bit seed

void moveSpacer(double* arr, double r, double theta, double phi, int valency);

double pegRete(int nPeg);

double stretch(double z0, double r, double z, int nPeg);

double stretchBound(double z0, double r, int nPeg);

double crossTest(const double* pos1, const double* pos2, int s1, int s2, int e1, int e2);

#endif // MOVES_H