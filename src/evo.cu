#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <cstring>
#include <stdio.h>
//#include <bits/stdc++.h>
//#include <unistd.h>
#include <cmath>
#include <cassert>
#include <chrono>
#include "io.hpp"
#include "suqa.cuh"
#include "system.cuh"
#include "Rand.hpp"


using namespace std;

void self_plaquette(const bmReg& qr0, const bmReg& qr1, const bmReg& qr2, const bmReg& qr3);

int main(int argc, char** argv){
    if(argc<4){
       printf("usage: %s <total_steps> <trotter_stepsize> <outfile>\n",argv[0]); 
       exit(1);
    }


   // g_beta = stod(argv[1]); // (extern) def in src/system.cu
    t_param = 1.0;
    U_param = 1.0;
    mu_param = 1.0;
    int total_steps = atoi(argv[1]);
    double trotter_stepsize = stod(argv[2]);
    string outfilename = argv[3];

	suqa::setup(syst_qbits);


    //printf("arguments:\n g_beta = %.16lg\n total_steps = %d\n trotter_stepsize = %.16lg\n outfile = %s\n", g_beta, total_steps, trotter_stepsize, outfilename.c_str());
    printf("arguments:\n total_steps = %d\n trotter_stepsize = %.16lg\n outfile = %s\n\n", total_steps, trotter_stepsize, outfilename.c_str());


    // init
    pcg rangen;
//    rangen.set_seed(12345); //time(NULL));

    FILE * outfile;

    init_state_testevo();
    DEBUG_CALL(printf("initial state:\n"));
    DEBUG_READ_STATE();

    for(uint ii=0; ii<=(uint)total_steps; ++ii){
        double t = ii*trotter_stepsize;
        DEBUG_CALL(printf("time %.16lg\n", t));

        // legit way to perform measures
//        for(uint hit=0; hit<(uint)num_hits; ++hit){
//            printf("\thit %u\n", hit);
//            init_state(state, Dim, g_beta);
//            evolution(state, t, ii);
//            self_plaquette(state, bm_qlink1, bm_qlink0, bm_qlink2, bm_qlink0);
//            std::vector<uint> c(3);
//            suqa::measure_qbit(state,bm_qlink1[0],c[0],rangen.doub());
//            suqa::measure_qbit(state,bm_qlink1[1],c[1],rangen.doub());
//            suqa::measure_qbit(state,bm_qlink1[2],c[2],rangen.doub());
//            uint plaq_idx = 4*c[2]+2*c[1]+c[0];
//            double plaq_tmp = (plaq_idx==0)? 2.0 : ((plaq_idx==2)? -2.0 : 0.0);
//            plaq_val += plaq_tmp;
//            plaq_val_std += plaq_tmp*plaq_tmp;
//        }
//        plaq_val /=(double)num_hits;
//        plaq_val_std = sqrt((plaq_val_std/(double)num_hits - plaq_val*plaq_val)/(double)(num_hits-1));
//        fprintf(outfile, "%.16lg %d %.16lg %.16lg\n", t, num_hits, plaq_val, plaq_val_std);

        init_state_testevo();
//		suqa::apply_h(bm_spin[rangen.randint(0,3)]);
//        printf("random number= %d\n", rangen.randint(0,3));
//        DEBUG_READ_STATE();
        evolution(t, ii);
        DEBUG_CALL(printf("after evolution by t=%lg:\n",t));
        DEBUG_READ_STATE();

        //suqa::apply_h(state,  bm_spin[rangen.randint(0,3)]);
	
        double p0000=0, p0001=0, p0010=0, p0011=0, p0100=0, p0101=0, p0110=0, p0111=0, p1000=0, p1001=0, p1010=0, p1011=0, p1100=0, p1101=0, p1110=0, p1111=0;
        suqa::prob_filter(bm_spin, {0U,0U,0U,0U}, p0000);
        suqa::prob_filter(bm_spin, {0U,0U,0U,1U}, p0001);
        suqa::prob_filter(bm_spin, {0U,0U,1U,0U}, p0010);
        suqa::prob_filter(bm_spin, {0U,0U,1U,1U}, p0011);
        suqa::prob_filter(bm_spin, {0U,1U,0U,0U}, p0100);
        suqa::prob_filter(bm_spin, {0U,1U,0U,1U}, p0101);
        suqa::prob_filter(bm_spin, {0U,1U,1U,0U}, p0110);
        suqa::prob_filter(bm_spin, {0U,1U,1U,1U}, p0111);
        suqa::prob_filter(bm_spin, {1U,0U,0U,0U}, p1000);
        suqa::prob_filter(bm_spin, {1U,0U,0U,1U}, p1001);
        suqa::prob_filter(bm_spin, {1U,0U,1U,0U}, p1010);
        suqa::prob_filter(bm_spin, {1U,0U,1U,1U}, p1011);
        suqa::prob_filter(bm_spin, {1U,1U,0U,0U}, p1100);
        suqa::prob_filter(bm_spin, {1U,1U,0U,1U}, p1101);
        suqa::prob_filter(bm_spin, {1U,1U,1U,0U}, p1110);
        suqa::prob_filter(bm_spin, {1U,1U,1U,1U}, p1111);
        DEBUG_CALL(printf("p0000 = %.12lg; p1000 = %.12lg\n", p0000, p1000));
        outfile = fopen(outfilename.c_str(), "a");
        if (ii == 0) {
           fprintf(outfile, "#t, p0000, p0001, p0010, p0011, p0100, p0101, p0110, p0111, p1000, p1001, p1010, p1011, p1100, p1101, p1110, p1111\n");
        }
        fprintf(outfile, "%.6lg %.6lg %.6lg %.6lg %.6lg %.6lg %.6lg %.6lg %.6lg %.6lg %.6lg %.6lg %.6lg %.6lg %.6lg %.6lg %.6lg\n", t, p0000,p0001,p0010,p0011,p0100,p0101,p0110,p0111,p1000,p1001,p1010,p1011,p1100,p1101,p1110,p1111);

        fclose(outfile);
    }
    
    suqa::clear();

    return 0;
}
