#include "system.cuh"
#include "Rand.hpp"
#include <vector>

double g_beta;
double t_param;
double U_param;
double mu_param;

void init_state(){
    suqa::init_state();
//	suqa::apply_x(bm_spin[1]);
//	suqa::apply_h(bm_spin[1]);
//	suqa::apply_cx(bm_spin[1], bm_spin[0]);
}

void init_state_testevo(){
    size_t dim = 1 << syst_qbits;

    std::vector<double> re_coeff(dim, 0.0);
    std::vector<double> im_coeff(dim, 0.0);

    double re_coeff_1 = 0.435162146493599;
    double re_coeff_2 = 0.557345410189304;

    re_coeff[5] = re_coeff_1;
    re_coeff[10] = re_coeff_1;
    re_coeff[6] = re_coeff_2;
    re_coeff[9] = re_coeff_2;

    suqa::init_state(re_coeff,im_coeff);
}


void exp_it_id_x_x(const bmReg& q, uint pos_id, double phase_t){

	suqa::apply_pauli_TP_rotation({q[(pos_id+1)%4],q[(pos_id+2)%4]}, {PAULI_X,PAULI_X}, phase_t);

}

// applies exp(-iHt)
void evolution(const double& t, const int& n){
    //(void)n;
    double N = (double) n;

    double fac_t = -t_param/2.0;
    double fac_U = U_param/4.0;
    double fac_mu = mu_param/2.0 - U_param/4.0;
    double fac_phase = (U_param/2.0 - 2.0*mu_param);

    fac_t = fac_t/N;
    fac_U = fac_U/N;
    fac_mu = fac_mu/N;
    fac_phase = fac_phase/N;

for (int it=0; it<N; it++) {
    // alpha term
    suqa::apply_pauli_TP_rotation({bm_spin[0],bm_spin[1]}, {PAULI_X,PAULI_X}, -t*fac_t);
    suqa::apply_pauli_TP_rotation({bm_spin[0],bm_spin[1]}, {PAULI_Y,PAULI_Y}, -t*fac_t);
    suqa::apply_pauli_TP_rotation({bm_spin[2],bm_spin[3]}, {PAULI_X,PAULI_X}, -t*fac_t);
    suqa::apply_pauli_TP_rotation({bm_spin[2],bm_spin[3]}, {PAULI_Y,PAULI_Y}, -t*fac_t);

    // beta term
    suqa::apply_pauli_TP_rotation({bm_spin[0],bm_spin[2]}, {PAULI_Z,PAULI_Z}, -t*fac_U);
    suqa::apply_pauli_TP_rotation({bm_spin[1],bm_spin[3]}, {PAULI_Z,PAULI_Z}, -t*fac_U);

    // gamma term
    suqa::apply_pauli_TP_rotation({bm_spin[0]}, {PAULI_Z}, -t*fac_mu);
    suqa::apply_pauli_TP_rotation({bm_spin[1]}, {PAULI_Z}, -t*fac_mu);
    suqa::apply_pauli_TP_rotation({bm_spin[2]}, {PAULI_Z}, -t*fac_mu);
    suqa::apply_pauli_TP_rotation({bm_spin[3]}, {PAULI_Z}, -t*fac_mu);

    // global phase
//    for (int s=0; s<4; s++) {
//    suqa::apply_x(bm_spin[s]);
//    suqa::apply_u1(bm_spin[s], -t*fac_phase);
//    suqa::apply_x(bm_spin[s]);
//    suqa::apply_u1(bm_spin[s], -t*fac_phase);
//    }

}

//	for (uint iii=0; iii<4; ++iii){
//		exp_it_id_x_x(bm_spin, iii, -t);
// 	}
}

// qsa specifics
void qsa_init_state(){
    suqa::init_state();
    suqa::apply_h(bm_spin[0]);
    suqa::apply_h(bm_spin[1]);
    suqa::apply_h(bm_spin[2]);
    suqa::apply_cx(bm_spin[0], bm_spin_tilde[0]);
    suqa::apply_cx(bm_spin[1], bm_spin_tilde[1]);
    suqa::apply_cx(bm_spin[2], bm_spin_tilde[2]);
}

void evolution_szegedy(const double& t, const int& n){
    (void)n;
      DEBUG_CALL(std::cout<<"before evolution_szegedy()"<<std::endl);
      DEBUG_READ_STATE();
      DEBUG_CALL(std::cout<<"apply evolution_szegedy()"<<std::endl);
    for (uint i = 0; i < 3; i++) {
      suqa::apply_pauli_TP_rotation({bm_spin_tilde[(0+i)%3],bm_spin_tilde[(1+i)%3]}, {PAULI_X,PAULI_X}, -t);
      DEBUG_CALL(std::cout<<"apply pauli 1 it "<<i<<std::endl);
      DEBUG_READ_STATE();
      suqa::apply_pauli_TP_rotation({bm_spin[(0+i)%3],bm_spin[(1+i)%3]}, {PAULI_X,PAULI_X}, t);
      DEBUG_CALL(std::cout<<"apply pauli 2 it "<<i<<std::endl);
      DEBUG_READ_STATE();
    }
}

void evolution_measure(const double& t, const int& n){
    (void)n;
  for (uint i = 0; i < 3; i++) {
    suqa::apply_pauli_TP_rotation({bm_spin[(0+i)%3],bm_spin[(1+i)%3]}, {PAULI_X,PAULI_X}, -t);
  }

}
void evolution_tracing(const double& t, const int& n){
    (void)n;
  for (uint i = 0; i < 3; i++) {
    suqa::apply_pauli_TP_rotation({bm_spin_tilde[(0+i)%3],bm_spin_tilde[(1+i)%3]}, {PAULI_X,PAULI_X}, -t);
  }

}

/* Measure facilities */
const uint op_bits = 4; // 2^op_bits is the number of eigenvalues for the observable
const bmReg bm_op = bm_spin; // where the measure has to be taken
const std::vector<double> op_vals = {0.0, 1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 3.0, 1.0, 2.0, 2.0, 3.0, 2.0, 3.0, 3.0, 4.0}; // eigvals
 
// change basis to the observable basis somewhere in the system registers
void apply_measure_rotation(){
//	suqa::apply_h(bm_spin[0]);
//	suqa::apply_h(bm_spin[1]);
//	suqa::apply_h(bm_spin[2]);
//	suqa::apply_cx(bm_spin[0], bm_spin[1]);
//	suqa::apply_cx(bm_spin[0], bm_spin[2]);
//	suqa::apply_u1(bm_spin[2], M_PI*0.5);
}

// inverse of the above function
void apply_measure_antirotation(){
 	apply_measure_rotation();
}

// map the classical measure recorded in creg_vals
// to the corresponding value of the observable;
// there is no need to change it
double get_meas_opvals(const uint& creg_vals){
    return op_vals[creg_vals];
}

// actually perform the measure
// there is no need to change it
double measure_X(pcg& rgen){

    std::vector<uint> classics(op_bits);
    
    apply_measure_rotation();

    std::vector<double> rdoubs(op_bits);
    for(auto& el : rdoubs){
        el = rgen.doub();
    }
    suqa::measure_qbits(bm_op, classics, rdoubs);

    apply_measure_antirotation();

    uint meas = 0U;
    for(uint i=0; i<op_bits; ++i){
        meas |= (classics[i] << i);
    }

    return get_meas_opvals(meas);

    return 0.0;
}

/* Moves facilities */

std::vector<double> C_weigthsums = {1./4, 2./4, 3./4, 1.0};

void apply_C(const uint &Ci){
    if(Ci>3)
        throw std::runtime_error("ERROR: wrong move selection");
    //suqa::apply_x(bm_spin[Ci]);
    suqa::apply_h(bm_spin[Ci]);
}

void apply_C_inverse(const uint &Ci){
    apply_C(Ci);
}

void qsa_apply_C(const uint &Ci){
  if(Ci>2) throw std::runtime_error("ERROR: wrong move selection");
  suqa::apply_h(bm_spin_tilde[Ci]);
// suqa::apply_h(state,bm_spin_tilde[(Ci+1)%3]);


  // suqa::apply_h(state,bm_spin_tilde);
}

void qsa_apply_C_inverse(const uint &Ci){
  if(Ci>2) throw std::runtime_error("ERROR: wrong move selection");
  //suqa::apply_h(state,bm_spin_tilde);
  //suqa::apply_h(state,bm_spin_tilde[(Ci+1)%3]);
  suqa::apply_h(bm_spin_tilde[Ci]);
}

std::vector<double> get_C_weigthsums(){ return C_weigthsums; }

