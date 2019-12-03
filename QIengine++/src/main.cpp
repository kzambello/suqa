#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <cstring>
#include <stdio.h>
#include <bits/stdc++.h>
#include <cmath>
#include <cassert>
#include "include/Rand.hpp"
#include <chrono>
#include "include/io.hpp"
#include "include/parser.hpp"
#include "include/suqa_gates.hpp"
#include "include/qms.hpp"

using namespace std;

void print_banner(){
    printf("\n"
"                                          \n" 
"    ███████╗██╗   ██╗ ██████╗  █████╗     \n" 
"    ██╔════╝██║   ██║██╔═══██╗██╔══██╗    \n" 
"    ███████╗██║   ██║██║   ██║███████║    \n" 
"    ╚════██║██║   ██║██║▄▄ ██║██╔══██║    \n" 
"    ███████║╚██████╔╝╚██████╔╝██║  ██║    \n" 
"    ╚══════╝ ╚═════╝  ╚══▀▀═╝ ╚═╝  ╚═╝    \n" 
"                                          \n" 
"\nSimulator for Universal Quantum Algorithms\n");
}



/* Hamiltonian
 *
 * H = 1/4 (1 + X1 X0 + X2 X0 + X2 X1)
 *
 */


// simulation parameters
double beta;
double h;

void init_state(){
    qms::gState.resize(qms::Dim);
    std::fill_n(qms::gState.begin(), qms::gState.size(), 0.0);
    qms::gState[0] = TWOSQINV; 
    qms::gState[3] = -TWOSQINV; 
}

arg_list args;

int main(int argc, char** argv){
    if(argc < 6){
        cout<<"arguments: <beta> <h> <metro steps> <reset each> <num_E_qbits> <output file path> [--max-reverse <max reverse attempts>=20] [--seed <seed>=random] [--PE-time <factor for time in PE (coeff. of 2pi)>=1.0] [--PE-steps <steps of PE evolution>=10] [--X-mat-stem <stem for X measure matrix>] [--record-reverse]"<<endl;
        exit(1);
    }

    parse_arguments(args, argc, argv);

    beta = args.beta;
    h = args.h;
    qms::metro_steps = (uint)args.metro_steps;
    qms::reset_each = (uint)args.reset_each;
    qms::ene_qbits = (uint)args.ene_qbits;
    string outfilename(args.outfile);
    qms::max_reverse_attempts = (uint)args.max_reverse_attempts;
    qms::t_PE_factor = args.pe_time_factor;
    qms::t_phase_estimation = qms::t_PE_factor*8.*atan(1.0); // 2*pi*t_PE_factor
    qms::n_phase_estimation = args.pe_steps;
    qms::Xmatstem = args.Xmatstem;
    qms::record_reverse= args.record_reverse;
    qms::iseed = args.seed;
    if(qms::iseed>0)
        qms::rangen.set_seed(qms::iseed);
    
    qms::iseed = qms::rangen.get_seed();

    qms::state_qbits = 3;
    qms::nqubits = qms::state_qbits + 2*qms::ene_qbits + 1;
    qms::Dim = (uint)pow(2, qms::nqubits);
    qms::ene_levels = (uint)pow(2, qms::ene_qbits);
    qms::state_levels = (uint)pow(2, qms::state_qbits);
    
    // Banner
    print_banner();
    cout<<"arguments:\n"<<args<<endl;

    auto t_start = std::chrono::high_resolution_clock::now();

    // Initialization of utilities
    qms::fill_rphase(qms::ene_qbits+1);
    qms::fill_bitmap();
    qms::fill_W_utils(beta, qms::t_PE_factor);
    if(qms::Xmatstem!="")
        qms::init_measure_structs();

    // Initialization:
    // known eigenstate of the system: psi=0, E_old = 0
    
    init_state();

    uint perc_mstep = qms::metro_steps/20;

    bool take_measure;
    uint s0 = 0U;
    for(uint s = 0U; s < qms::metro_steps; ++s){
        take_measure = (s>s0 and (s-s0)%qms::reset_each ==0U);
        int ret = qms::metro_step(take_measure);

        if(ret<0){ // failed rethermalization, reinitialize state
            init_state();
            //ensure new rethermalization
            s0 = s+1; 
        }
        if(s%perc_mstep==0){
#ifdef NDEBUG
            cout<<("\riteration: "+to_string(s)+"/"+to_string(qms::metro_steps));
            cout.flush();
#else
            cout<<("iteration: "+to_string(s)+"/"+to_string(qms::metro_steps))<<endl;
#endif
        }
    }
    cout<<endl;

    cout<<"all fine :)\n"<<endl;

    FILE * fil = fopen(outfilename.c_str(), "w");

    fprintf(fil, "# it E%s\n",(qms::Xmatstem!="")?" A":"");

    for(uint ei = 0; ei < qms::E_measures.size(); ++ei){
        if(qms::Xmatstem!=""){
            fprintf(fil, "%d %.16lg %.16lg\n", ei, qms::E_measures[ei], qms::X_measures[ei]);
        }else{
            fprintf(fil, "%d %.16lg\n", ei, qms::E_measures[ei]);
        }
    }
    fclose(fil);

    if(qms::record_reverse){
        FILE * fil_rev = fopen((outfilename+"_revcounts").c_str(), "w");


        for(uint i = 0; i < qms::reverse_counters.size(); ++i){
            fprintf(fil_rev, "%d %d\n", i, (int)qms::reverse_counters[i]);
        }
        fclose(fil_rev);
    }

    cout<<"\n\tSuqa!\n"<<endl;

    auto t_end = std::chrono::high_resolution_clock::now();
    double secs_passed = (1./1000.)*std::chrono::duration<double, std::milli>(t_end-t_start).count();
	cout<<"All [DONE] in "<<secs_passed<<" seconds"<<endl;

    return 0;
}
