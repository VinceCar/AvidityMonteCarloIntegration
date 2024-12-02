#include "ParseInput.h"
#include "DispComb.h"
#include "Moves.h"
#include <iostream>
#include <fstream>
#include <string>
#include <openacc.h>
#include <cmath>
#include <omp.h>  // Include OpenMP header
#include <chrono>
#include <iomanip>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << std::endl;
        return 1;
    }

    std::string inputFilePath = argv[1];
    std::string outputFilePath = argv[2];

    try {
        ParseInput* Collector = parseFilesInParallel(inputFilePath);
        std::vector<std::string> filenames = readFilenames(inputFilePath);
        size_t numConfigs = filenames.size();
        
    // Get the current time as a time_point
    auto now = std::chrono::system_clock::now();

    // Convert to a time_t, which represents calendar time
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);

    // Convert the time to a human-readable format
    std::tm* local_time = std::localtime(&now_time);
    unsigned int i, aa, aaa, bb, cc, dd, nbound, ff, gg, lf, ee;
    uint64_t lint_ord, rint_ord;


    // Print the time in a readable format
    std::cout << "Current Time: " << std::put_time(local_time, "%Y-%m-%d %H:%M:%S") << std::endl;
    // Perform the main parallel computation on the GPU
    #pragma omp parallel for num_threads(72) private( i, aa, aaa, bb, cc, dd, nbound, ff, gg, ee, lint_ord, rint_ord, lf)            
    for (i = 0; i <  numConfigs ; i++) {

        if (i < numConfigs) 
        {
            for(aa=0;aa<Collector[i].nstep;aa++)
            {

                //set ligand core and receptor positions
                //the spacer positions are the points, where the spacer is connected to the ligand core
                setSpacer(Collector[i].spacer,Collector[i].dcore,Collector[i].valency_lig);  // the initial position of the core is just flat right at the origin; the core position will 
                setSpacer(Collector[i].receptor,Collector[i].drec,Collector[i].valency_rec); // receptor positions are generated with the same function; the receptor positions will not be

                //rotate ligand core
                Collector[i].phi0=random(2.0*Collector[i].pi/Collector[i].valency_lig,&Collector[i].ran_seed);
                rotateZ(Collector[i].spacer,Collector[i].phi0,Collector[i].valency_lig);
                Collector[i].phi0=random(Collector[i].pi,&Collector[i].ran_seed);
                Collector[i].phirot=Collector[i].phi0;
                rotateY(Collector[i].spacer,Collector[i].phi0,Collector[i].valency_lig);
                Collector[i].phi0=random(2.0*Collector[i].pi,&Collector[i].ran_seed);
                rotateX(Collector[i].spacer,Collector[i].phi0,Collector[i].valency_lig);

                //move ligand core relative to the center of the receptor
                Collector[i].theta0=random(Collector[i].pi/2.0,&Collector[i].ran_seed);
                Collector[i].phi0=random(2.0*Collector[i].pi,&Collector[i].ran_seed);
                Collector[i].r0=random(Collector[i].movespacermax,&Collector[i].ran_seed);
                moveSpacer(Collector[i].spacer,Collector[i].r0,Collector[i].theta0,Collector[i].phi0,Collector[i].valency_lig);

                //check that none of the core corners are below 0 (i.e. that the core does not intersect with the receptor)
                Collector[i].testspacer=0;
                for(aaa=0;aaa<Collector[i].valency_lig;aaa++)if(Collector[i].spacer[aaa*3+2]<0.0)Collector[i].testspacer=1;

                // everything occurs if the configuration is not rejected
                if(Collector[i].testspacer==0)
                {
                    //unbound ligand units
                    for(bb=0;bb<Collector[i].valency_lig;bb++)
                    {
                        Collector[i].theta[bb]=random(Collector[i].pi,&Collector[i].ran_seed);
                        Collector[i].r[bb]=random(Collector[i].stretchspacermax[bb],&Collector[i].ran_seed);
                        Collector[i].phi[bb]=random(2.0*Collector[i].pi,&Collector[i].ran_seed);
                        Collector[i].z=Collector[i].spacer[bb*3+2]+Collector[i].r[bb]*cos(Collector[i].theta[bb]);
                        if(Collector[i].z>0)
                        {
                            Collector[i].freeLU[bb]=sin(Collector[i].theta[bb])*pow(Collector[i].r[bb],2)*stretch(Collector[i].spacer[bb*3+2],Collector[i].r[bb]*sin(Collector[i].theta[bb]),Collector[i].z,Collector[i].npeg[bb]);
                        }
                        else {Collector[i].freeLU[bb]=0.0;}
                        Collector[i].end[bb*3+0]=Collector[i].r[bb]*sin(Collector[i].theta[bb])*cos(Collector[i].phi[bb]);
                        Collector[i].end[bb*3+1]=Collector[i].r[bb]*sin(Collector[i].theta[bb])*sin(Collector[i].phi[bb]);
                        Collector[i].end[bb*3+2]=Collector[i].r[bb]*cos(Collector[i].theta[bb]);
                    }
                    //bound ligand
                    for(cc=0;cc<Collector[i].valency_lig;cc++)
                    {
                        for(dd=0;dd<Collector[i].valency_rec;dd++)
                        {
                            Collector[i].r[cc]=sqrt(pow(Collector[i].receptor[dd*3+0]-Collector[i].spacer[cc*3+0],2)+pow(Collector[i].receptor[dd*3+1]-Collector[i].spacer[cc*3+1],2));
                            Collector[i].boundLU[cc*Collector[i].valency_rec+dd]=stretchBound(Collector[i].spacer[cc*3+2],Collector[i].r[cc],Collector[i].npeg[cc]);
                        }
                    }

                    // Continue here with the cycle over all lint and rint
                    for(nbound = 1; nbound <= Collector[i].max_k; nbound++)
                    {
                        for(lint_ord = 0; lint_ord < Collector[i].rows_comb_lint[nbound -1]; lint_ord++)
                        {
                            for(rint_ord = 0; rint_ord < Collector[i].rows_disp_rint[nbound -1]; rint_ord++)
                            {
                               Collector[i].frees=1.0;
                               for(lf = 0; lf<Collector[i].valency_lig; lf++)
                               {
                                    Collector[i].test=0;
                                    for(ee=0;ee<nbound;ee++)
                                    {
                                        Collector[i].lint_comb_id=Collector[i].cumsum_lint[nbound -1]+lint_ord*nbound+ee;
                                        if(lf==Collector[i].lint_comb_container[Collector[i].lint_comb_id])Collector[i].test=1;
                                    }
                                    if(Collector[i].test==0) {
                                        Collector[i].frees*=Collector[i].freeLU[lf]*2.0*Collector[i].pi*Collector[i].pi*Collector[i].stretchspacermax[lf]; 
                                    }
                                }
                                Collector[i].bound=1.0;
                                for (ff = 0; ff < nbound; ff++) {
                                    Collector[i].lint_comb_id = Collector[i].cumsum_lint[nbound -1] + lint_ord * nbound + ff;
                                    Collector[i].rint_disp_id = Collector[i].cumsum_rint[nbound -1] + rint_ord * nbound + ff;
                                    Collector[i].bound *= Collector[i].boundLU[Collector[i].lint_comb_container[Collector[i].lint_comb_id] * Collector[i].valency_lig + Collector[i].rint_disp_container[Collector[i].rint_disp_id]] / Collector[i].kD[Collector[i].lint_comb_container[Collector[i].lint_comb_id]];
                                }
                                Collector[i].testcross=0.0;
                                if (Collector[i].testcross == 0) {
                                    Collector[i].partitionfunction[nbound -1] += sin(Collector[i].phirot) * sin(Collector[i].theta0) * pow(Collector[i].r0, 2) * Collector[i].bound * Collector[i].frees;
                                }
                            }
                        }
                    }
                }
            }
            for(gg=0;gg<Collector[i].max_k;gg++)Collector[i].Q[gg]=Collector[i].pi/2.0*Collector[i].pi*Collector[i].pi*Collector[i].movespacermax*pow(1.66,gg)*Collector[i].partitionfunction[gg]/Collector[i].nstep;
        }
    }
    now = std::chrono::system_clock::now();

    // Convert to a time_t, which represents calendar time
    now_time = std::chrono::system_clock::to_time_t(now);

    // Convert the time to a human-readable format
    local_time = std::localtime(&now_time);

    // Print the time in a readable format
    std::cout << "Current Time: " << std::put_time(local_time, "%Y-%m-%d %H:%M:%S") << std::endl;
 
        std::cout << "Test" << std::endl;
        std::cout << numConfigs << std::endl;
        std::cout << "Test passed: All files parsed and contents verified successfully." << std::endl;
        writeResultsToFile(Collector, numConfigs, outputFilePath);


        // Clean up dynamically allocated memory on the CPU
        delete[] Collector;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
