#ifndef PARSEINPUT_H
#define PARSEINPUT_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <thread>

// Struct for ParseInput
struct ParseInput {
    unsigned int valency_lig, valency_rec, max_k;
    uint64_t nstep, ran_seed, ran_seed_init;
    double drec, dcore;
    double* kD;   // Pointer to hold dynamic array of floats
    unsigned int* npeg;    // Pointer to hold dynamic array of integers
    // Loop indices
    uint64_t i, j, k, ii, jj, m, nb, lf, lint_ord, rint_ord;

    // Random number generation variables
    double phi0, phirot, theta0, r0, z;
    
    // Constants
    double pi = 3.141592653589793;

    // Spacer and receptor positions
    double* spacer;        // Pointer to dynamically allocated array
    double* receptor;      // Pointer to dynamically allocated array

    // Movement limits
    double movespacermax;
    double* stretchspacermax;   // Pointer to dynamically allocated array
    
    // Free ligand unit properties
    double* theta;         // Pointer to dynamically allocated array
    double* r;             // Pointer to dynamically allocated array
    double* phi;           // Pointer to dynamically allocated array
    double* freeLU;        // Pointer to dynamically allocated array
    double* end;           // Pointer to dynamically allocated array
    double* rete;
    // Bound ligand unit properties
    double* boundLU;       // Pointer to dynamically allocated array
    
    // Control flags
    double testspacer, test, testcross;
    
    // Variables for combination and interaction
    unsigned int nbound, level;
    uint64_t* rows_comb_lint;   // Pointer to dynamically allocated array
    uint64_t* rows_disp_rint;   // Pointer to dynamically allocated array
    uint64_t* cumsum_lint; // Pointer to dynamically allocated array
    uint64_t* cumsum_rint; // Pointer to dynamically allocated array
    unsigned int* lint_comb_container;   // Pointer to dynamically allocated array
    unsigned int* rint_disp_container;   // Pointer to dynamically allocated array

    unsigned int lint_comb_id;
    unsigned int rint_disp_id;    
    // Energy and partition function variables
    double frees, bound;
    double* Q; 
    double* partitionfunction;  // Pointer to dynamically allocated array
};

// Function to parse the configuration file
ParseInput parseInputFile(const std::string& filename);

// Helper function to split a string by a delimiter
char** split(const std::string& str, char delimiter, int& count);

// Function to read filenames from a text file
std::vector<std::string> readFilenames(const std::string& listFile);

// Function to read and parse files concurrently using threads
ParseInput* parseFilesInParallel(const std::string& listFile);

void writeResultsToFile(const ParseInput* results, std::size_t count, const std::string& filename);

#endif // PARSEINPUT_H
