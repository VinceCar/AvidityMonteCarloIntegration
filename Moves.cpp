#include "Moves.h"
#include <iostream>
#include <cstddef> 
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstdint>  // For uint64_t

// Custom Linear Congruential Generator using 64-bit unsigned integers
uint64_t lcg_random(uint64_t *seed) {
    *seed = (6364136223846793005ULL * (*seed) + 1ULL);  // Example LCG for 64-bit
    return *seed;
}

// Generate a random double in the range [0, factor) using a 64-bit seed
double random_double(uint64_t *seed, double factor) {
    return factor * ((double)lcg_random(seed) / 18446744073709551615.0);
}

// Wrapper for random double
double random(double factor, uint64_t *seed) {
    return random_double(seed, factor);
}

// Other functions (unchanged except for using `seq` where needed)
void setSpacer(double* arr, double dCore, int valency)
{
    // Constants that don’t change during the loop
    double pi = 3.14159265359;
    double factor = 2.0 / valency * pi;  // Precompute factor to reduce register use

    // Loop over the valency
    for (int i = 0; i < valency; i++) {
        double angle = i * factor;  // Only compute this once per iteration
        double cos_val = cos(angle);  // Precompute cosine
        double sin_val = sin(angle);  // Precompute sine
        
        arr[i * 3 + 0] = dCore * cos_val;
        arr[i * 3 + 1] = dCore * sin_val;
        arr[i * 3 + 2] = 0.0;
    }
    return;
}

void rotateX(double* arr, double phi, int valency)
{
    for (int i = 0; i < valency; i++) {
        double y = cos(phi) * arr[i*3+1] + sin(phi) * arr[i*3+2];
        double z = -sin(phi) * arr[i*3+1] + cos(phi) * arr[i*3+2];
        arr[i*3+1] = y;
        arr[i*3+2] = z;
    }
    return;
}

void rotateY(double* arr, double phi, int valency)
{
    for (int i = 0; i < valency; i++) {
        double x = cos(phi) * arr[i*3+0] - sin(phi) * arr[i*3+2];
        double z = sin(phi) * arr[i*3+0] + cos(phi) * arr[i*3+2];
        arr[i*3+0] = x;
        arr[i*3+2] = z;
    }
    return;
}

void rotateZ(double* arr, double phi, int valency)
{
    for (int i = 0; i < valency; i++) {
        double x = cos(phi) * arr[i*3+0] + sin(phi) * arr[i*3+1];
        double y = -sin(phi) * arr[i*3+0] + cos(phi) * arr[i*3+1];
        arr[i*3+0] = x;
        arr[i*3+1] = y;
    }
    return;
}

void moveSpacer(double* arr, double r, double theta, double phi, int valency)
{
    for (int i = 0; i < valency; i++) {
        arr[i*3+0] += r * sin(theta) * cos(phi);
        arr[i*3+1] += r * sin(theta) * sin(phi);
        arr[i*3+2] += r * cos(theta);
    }
    return;
}

double pegRete(int nPeg)
{
    return 0.4 * pow(nPeg - 1, 0.6);
}

double stretch(double z0, double r, double z, int nPeg)
{
    double rete = pegRete(nPeg);
    double retes = pegRete(nPeg);
    double prob = 1.0 / (3.031013847264815 * pow(rete, 3)) * exp(-1.5 * pow(r / rete, 2)) * 
                  (exp(-1.5 * pow((z - z0) / retes, 2)) - exp(-1.5 * pow((z + z0) / retes, 2)));
    return prob < 0 ? 0 : prob;
}

double stretchBound(double z0, double r, int nPeg)
{
    double rete = pegRete(nPeg);
    double factor = 0.448695;
    return 4.0 / 3.14159265359 * factor / pow(rete, 3) * exp(-1.5 * (pow(r, 2) + pow(z0, 2)) / pow(rete, 2)) * z0 / pow(rete, 2);
}

double crossTest(const double* pos1, const double* pos2, int s1, int s2, int e1, int e2)
{
    double dist12 = pow(pos1[s1*3+0] - pos2[e2*3+0], 2) + pow(pos1[s1*3+1] - pos2[e2*3+1], 2) + pow(pos1[s1*3+2] - pos2[e2*3+2], 2);
    double dist21 = pow(pos1[s2*3+0] - pos2[e1*3+0], 2) + pow(pos1[s2*3+1] - pos2[e1*3+1], 2) + pow(pos1[s2*3+2] - pos2[e1*3+2], 2);
    double dist11 = pow(pos1[s1*3+0] - pos2[e1*3+0], 2) + pow(pos1[s1*3+1] - pos2[e1*3+1], 2) + pow(pos1[s1*3+2] - pos2[e1*3+2], 2);
    double dist22 = pow(pos1[s2*3+0] - pos2[e2*3+0], 2) + pow(pos1[s2*3+1] - pos2[e2*3+1], 2) + pow(pos1[s2*3+2] - pos2[e2*3+2], 2);
    
    return (dist12 + dist21 < 0.96 * (dist11 + dist22)) ? 1.0 : 0.0;
}
