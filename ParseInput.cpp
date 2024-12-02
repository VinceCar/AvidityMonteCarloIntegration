#include "ParseInput.h"
#include "ThreadPool.h"  // Assuming the above code is saved in ThreadPool.h
#include "DispComb.h"
#include "Moves.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <cctype>   // For std::isspace
#include <cstdlib>  // For std::strtof, std::strtol
#include <cstring>  // For strtok
#include <thread>
#include <mutex>

// Helper function to check if a character is not a space
bool notSpace(unsigned char ch) {
    return !std::isspace(ch);
}

// Function to split a string by a delimiter and return an array of tokens
char** split(const std::string& str, char delimiter, int& count) {
    count = 0;
    std::istringstream tokenStream(str);
    std::string token;
    
    // First, count the number of tokens
    while (std::getline(tokenStream, token, delimiter)) {
        ++count;
    }

    // Allocate an array of C-style strings
    char** tokens = new char*[count];

    // Parse again and fill the array
    tokenStream.clear();
    tokenStream.str(str);
    int index = 0;
    while (std::getline(tokenStream, token, delimiter)) {
        tokens[index] = new char[token.size() + 1];
        std::strcpy(tokens[index], token.c_str());
        ++index;
    }
    return tokens;
}

// Helper function to trim whitespace from the beginning and end of a string
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(' ');
    if (std::string::npos == first) {
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, last - first + 1);
}

// Function to parse the configuration file
ParseInput parseInputFile(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    ParseInput Collect;
    Collect.kD = nullptr;  // Ensure pointers are initialized to null
    Collect.npeg = nullptr;

    while (std::getline(infile, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        // Split the line by '=' into key and value
        size_t delimiterPos = line.find('=');
        if (delimiterPos == std::string::npos) {
            continue;  // Skip invalid lines
        }

        std::string key = trim(line.substr(0, delimiterPos));
        std::string value = trim(line.substr(delimiterPos + 1));

        try {
            if (key == "valency_lig") {
                Collect.valency_lig = static_cast<unsigned int>(std::stoul(value));
            } else if (key == "valency_rec") {
                Collect.valency_rec = static_cast<unsigned int>(std::stoul(value));
            } else if (key == "nstep") {
                Collect.nstep = std::stoll(value);
            } else if (key == "drec") {
                Collect.drec = std::stof(value);
            } else if (key == "dcore") {
                Collect.dcore = std::stof(value);
            } else if (key == "ran_seed") {
                Collect.ran_seed = std::stoll(value);
            } else if (key == "kD") {
                int count;
                char** values = split(value, ',', count);  // Split the kD values by commas
                Collect.kD = new double[count];
                for (int i = 0; i < count; i++) {
                    Collect.kD[i] = std::strtof(values[i], nullptr);  // Convert each value to float
                    delete[] values[i];  // Free memory for each value
                }
                delete[] values;  // Free the array of pointers
            } else if (key == "npeg") {
                int count;
                char** values = split(value, ',', count);  // Split the npeg values by commas
                Collect.npeg = new unsigned int[count];
                for (int i = 0; i < count; i++) {
                    Collect.npeg[i] = static_cast<unsigned int>(std::stoul(values[i]));  // Convert each value to int
                    delete[] values[i];  // Free memory for each value
                }
                delete[] values;  // Free the array of pointers
            } else {
            }
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument for key: " << key << ", value: " << value << std::endl;
            throw;
        } catch (const std::out_of_range& e) {
            std::cerr << "Value out of range for key: " << key << ", value: " << value << std::endl;
            throw;
        }
    }

    infile.close();

    if (!Collect.kD || !Collect.npeg) {
        throw std::runtime_error("kD or npeg arrays are not allocated properly.");
    }

    unsigned int ii, jj;
    Collect.ran_seed_init=Collect.ran_seed;
    Collect.max_k = std::min(Collect.valency_lig, Collect.valency_rec);

    // Dynamically allocate 1D arrays
    Collect.rete = (double*) malloc(Collect.valency_lig * sizeof(double));
    Collect.stretchspacermax = (double*) malloc(Collect.valency_lig * sizeof(double));
    Collect.phi = (double*) malloc(Collect.valency_lig * sizeof(double));
    Collect.theta = (double*) malloc(Collect.valency_lig * sizeof(double));
    Collect.r = (double*) malloc(Collect.valency_lig * sizeof(double));
    Collect.partitionfunction = (double*) malloc(Collect.valency_lig * sizeof(double));
    Collect.Q = (double*) malloc(Collect.valency_lig * sizeof(double));
    Collect.freeLU = (double*) malloc(Collect.valency_lig * sizeof(double));

    //Dynamically allocate linearised 2D arrays for spacer and boundLU
    Collect.spacer = (double*) malloc(Collect.valency_lig * 3 * sizeof(double));        // 3 elements for each "row"
    Collect.receptor = (double*) malloc(Collect.valency_rec * 3 * sizeof(double));        // 3 elements for each "row"
    Collect.end = (double*) malloc(Collect.valency_lig * 3 * sizeof(double));        // 3 elements for each "row"

    Collect.boundLU = (double*) malloc(Collect.valency_lig * Collect.valency_rec * sizeof(double));

    // Declare other variables
    double dummymaxspacer;
   
    //Initialize support vectors containing indeces of ligands and receptors
    unsigned int* vec_lig = (unsigned int*) malloc(Collect.valency_lig * sizeof(unsigned int));
    unsigned int* vec_rec = (unsigned int*) malloc(Collect.valency_rec * sizeof(unsigned int));


    for (Collect.i = 0; Collect.i < Collect.valency_lig; Collect.i++) {
        vec_lig[Collect.i] = Collect.i;  // Manually implementing std::iota
    }
    for (Collect.i = 0; Collect.i < Collect.valency_rec; Collect.i++) {
        vec_rec[Collect.i] = Collect.i;  // Manually implementing std::iota
    }

    // Memory to store row sizes for combinations and dispositions
    Collect.rows_comb_lint = (uint64_t*) malloc(Collect.max_k * sizeof(uint64_t));
    Collect.rows_disp_rint = (uint64_t*) malloc(Collect.max_k * sizeof(uint64_t));
    // Initialize variables to count the total number of elements
    Collect.cumsum_lint = (uint64_t*) malloc((Collect.max_k+1) * sizeof(uint64_t));
    Collect.cumsum_rint = (uint64_t*) malloc((Collect.max_k+1) * sizeof(uint64_t));
    Collect.cumsum_lint[0]=0;
    Collect.cumsum_rint[0]=0;

    for (ii = 1; ii <= Collect.max_k; ii++)
    {
        Collect.rows_comb_lint[ii-1] = binomialCoefficient(Collect.valency_lig, ii);
        Collect.rows_disp_rint[ii-1] = permutation(Collect.valency_rec, ii) ;
	Collect.cumsum_lint[ii] = Collect.cumsum_lint[ii-1] + Collect.rows_comb_lint[ii-1] * ii;
        Collect.cumsum_rint[ii] = Collect.cumsum_rint[ii-1] + Collect.rows_disp_rint[ii-1] * ii;
    }

    Collect.lint_comb_container = (unsigned int*) malloc(Collect.cumsum_lint[Collect.max_k] * sizeof(unsigned int));
    Collect.rint_disp_container = (unsigned int*) malloc(Collect.cumsum_rint[Collect.max_k] * sizeof(unsigned int));
    
    uint64_t counter_lint = 0;
    uint64_t counter_rint = 0;

    for (ii = 1; ii <= Collect.max_k; ii++)
    {
        unsigned int* genComb = generateCombinations(vec_lig, Collect.valency_lig, ii, Collect.rows_comb_lint[ii-1]);
        unsigned int* genDisp = generateDispositions(vec_rec, Collect.valency_rec, ii, Collect.rows_disp_rint[ii-1]);
        counter_lint = 0;
        counter_rint = 0;
        for (jj = Collect.cumsum_lint[ii-1]; jj <= Collect.cumsum_lint[ii]-1; jj++)
        {
            Collect.lint_comb_container[jj] = genComb[counter_lint];
            counter_lint += 1;
        }
        for (jj = Collect.cumsum_rint[ii-1]; jj <= Collect.cumsum_rint[ii]-1; jj++)
        {
            Collect.rint_disp_container[jj] = genDisp[counter_rint];
            counter_rint += 1;
        }
        free(genComb);
        free(genDisp);
    }

    Collect.movespacermax = 0.0;
//start PEG-loop, you are doing many simulations based on the pegloop
    for (ii = 0; ii < Collect.valency_lig; ii++)
    {
        Collect.rete[ii]=pegRete(Collect.npeg[ii]);
        Collect.stretchspacermax[ii]=0.3 * Collect.npeg[ii];
        dummymaxspacer = Collect.dcore + Collect.stretchspacermax[ii];
        if (Collect.movespacermax<dummymaxspacer) Collect.movespacermax=dummymaxspacer;
    }
//PARTITION FUNTION
    for(Collect.i=0; Collect.i<Collect.valency_lig; Collect.i++) Collect.partitionfunction[Collect.i]=0.0;
    for(Collect.i=0; Collect.i<Collect.valency_lig; Collect.i++) Collect.Q[Collect.i]=0.0;

    return Collect;
}

// Function to read filenames from a text file
std::vector<std::string> readFilenames(const std::string& listFile) {
    std::ifstream infile(listFile);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file: " + listFile);
    }

    std::vector<std::string> filenames;
    std::string filename;
    while (std::getline(infile, filename)) {
        if (!filename.empty()) {
            filenames.push_back(filename);
        }
    }

    infile.close();
    return filenames;
}

// Mutex for thread safety when adding to the array
std::mutex vec_mutex;
ParseInput* parseFilesInParallel(const std::string& listFile) {
// Function to parse files in parallel using threads
    std::vector<std::string> filenames = readFilenames(listFile);
    size_t numFiles = filenames.size();
    
    // Create a thread pool with the number of threads based on available hardware concurrency
    ThreadPool pool(std::thread::hardware_concurrency());
    
    // Dynamically allocate an array of ParseInput objects
    ParseInput* parsedInputs = new ParseInput[numFiles];

    // Mutex for thread safety when writing to parsedInputs


    // Vector to store the futures of the results
    std::vector<std::future<void>> futures;

    // Lambda function to parse a file and store the result in the array
    auto parseTask = [&](const std::string& filename, int index) {
        // Parse the file and store the result
        ParseInput result = parseInputFile(filename);

        // Use mutex to safely store the result in the parsedInputs array
        std::lock_guard<std::mutex> lock(vec_mutex);
        parsedInputs[index] = result;
    };

    // Submit tasks to the thread pool to process each file
    for (size_t i = 0; i < filenames.size(); ++i) {
        futures.emplace_back(
            pool.enqueue([&, i] {
                parseTask(filenames[i], i);
            })
        );
    }

    // Wait for all tasks to complete
    for (auto& future : futures) {
        future.get();  // Ensure each task completes
    }

    return parsedInputs;
}

// Function to write the results to a file
void writeResultsToFile(const ParseInput* results, std::size_t count, const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // Loop through each ResultsStruct object
    for (std::size_t idx = 0; idx < count; ++idx) {
        const ParseInput& result = results[idx];

        // Write dcore, drec, valency_lig, valency_rec
        outFile << result.dcore << " " << result.drec << " " << result.valency_lig << " " << result.valency_rec << " ";

        // Write all npeg entries (size determined by valency_lig)
        for (unsigned int i = 0; i < result.valency_lig; i++) {
            outFile << result.npeg[i] << " ";
        }

        // Write all kD entries (size determined by valency_lig)
        for (unsigned int i = 0; i < result.valency_lig; i++) {
            outFile << result.kD[i] << " ";
        }

        // Write all Q entries (size determined by valency_lig)
        for (unsigned int i = 0; i < result.valency_lig; i++) {
            outFile << result.Q[i] << " ";
        }

        // End the line for this result
        outFile << std::endl;
    }

    outFile.close();
    std::cout << "Results written to " << filename << " successfully." << std::endl;
}
