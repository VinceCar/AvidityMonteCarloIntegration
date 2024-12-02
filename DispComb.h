#ifndef DISPCOMB_H
#define DISPCOMB_H

#include <vector>
#include <iostream>
#include <cstddef> 

// Function to generate all combinations
unsigned int* generateCombinations(unsigned int* arr, unsigned int n, unsigned int k, uint64_t numCombinations);

unsigned int* generateDispositions(unsigned int* arr, unsigned int n, unsigned int k, uint64_t numResults);

uint64_t binomialCoefficient(unsigned int n, unsigned int k);

uint64_t permutation(unsigned int n, unsigned int k);

#endif // DISPCOMB_H
