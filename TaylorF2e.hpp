/*
 * TaylorF2e.hpp
 *
 *  Created on: Feb 26, 2019
 *      Author: blakemoore
 */

#ifndef TAYLORF2E_HPP_
#define TAYLORF2E_HPP_

#include "Phase_Coef.hpp"
#include "Amps.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <complex>

namespace std {

class TaylorF2e {
public:
	TaylorF2e();
//	TaylorF2e(double e, double p, double M, double eta, double psi, double phi, double thet, double iot, double bet, double f0, double fend, double df);
	TaylorF2e(double M, double eta, double e, double p, double ampre, double ampim, double f0, double fend, double df);
	TaylorF2e(double M, double eta, double e, double p, double ampmag, double f0, double fend, double df);
	TaylorF2e(double M_in, double eta_in, double e_in, double ampmag, double f0_in, double fend_in, double df_in);
	TaylorF2e(double M_in, double eta_in, double e_in, double thet_in, double phi_in, double psi_in, double iot_in, double lamc_in, double lc_in, double tc_in, double D_in, double f0_in, double fend_in, double df_in);
	virtual ~TaylorF2e();

	void init_interps(int N);
	double get_fn_e(double e);
	double get_fw_e(double e);
	double get_y_e(double e);
	double get_e_fn(double fn);
	double fin_cond(double e);
	double get_p_e(double e);
	double fourier_f_minus_e(double e, int j);
	double fourier_f_plus_e(double e, int j);
	double fourier_f_s_e(double e, int s);
	double get_e_fin();
	double stat_e_s(double f, int s);
	double cond_minus(double& e, double& f, int& j);
	double cond_plus(double e, double f, int j);
	double stat_e_j_minus(double& f, int& j);
	double amplookup_j(double& e, double& y, int& j);
	double amplookup_s(double e, int s);
	double e_j_10hz(int& j);
	complex<double> h_j_plus(double f, int j);
	complex<double> h_j_minus(double& f, int& j);
	complex<double> h_s(double f, int j);
	void make_scheme();
	double raw_freq_to_disc_ind(double f);
	void make_F2e_min();
	vector<vector<complex<double>>> get_F2e_min();
	void set_e_fin();
	void make_F2e_summed();
	vector<complex<double>> get_F2e_min_plus();
	vector<complex<double>> get_F2e_min_cross();
	complex<double> h_j_minus_no_sky(double& f, int& j);
	void make_F2e_min_plus_cross();
	int calls;
	int count;
	int calls_guess;


private:
	double msun;
	double e0;
	double p0;
	double M;
	double eta;
	double psi;
	double phi;
	double thet;
	double iot;
	double bet;
	double F_p;
	double F_c;
	double f0;
	double fend;
	double df;
	double t_c;
	double l_c;
	double lam_c;
	complex<double> Q;
	complex<double> over_amp;
	double y0;
	vector<double> C_vec;
	vector<double> t_vec;
	vector<double> lam_vec;
	vector<double> l_vec;
	vector<double> y_vec;
	gsl_spline *spline_fn_e;
	gsl_spline *spline_y_e;
	gsl_spline *spline_e_fn;
	gsl_spline *spline_fw_e;
	gsl_interp_accel *acc_fn_e;
	gsl_interp_accel *acc_y_e;
	gsl_interp_accel *acc_e_fn;
	gsl_interp_accel *acc_fw_e;
	double elast;
	double F_n_last;
	double F_n_in;
	double DL;
	double e_fin;
	int j_min_max;
	int j_plus_max;
	int s_max;
	vector<vector<int>> j_min_range;
	vector<vector<int>> j_plus_range;
	vector<vector<int>> s_range;
	vector<vector<complex<double>>> F2_min;
	vector<double> phase_container;
	vector<double> e_10hz;
	double e_stat_last;
	double toff;
	double loff;
	double lamoff;
	vector<complex<double>> F2_summed;
	vector<complex<double>> F2_min_plus;
	vector<complex<double>> F2_min_cross;
};

} /* namespace std */

#endif /* TAYLORF2E_HPP_ */
