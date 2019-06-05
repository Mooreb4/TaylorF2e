/*
 * TaylorF2e.cpp
 *
 *  Created on: Feb 26, 2019
 *      Author: blakemoore
 */

#include "TaylorF2e.hpp"
namespace std {

TaylorF2e::TaylorF2e() {
	// TODO Auto-generated constructor stub
	msun = 4.925502303934785*pow(10, -6);
	e0 = 0.1;
	p0 = 50;
	y0 = 1/sqrt(p0);
	M = 20*msun;
	eta = 0.25;
	psi = 3./7.*M_PI;
	phi = 3./7.*M_PI;
	thet = 3./7.*M_PI;
	iot = 3./7.*M_PI;
	bet = 3./7.*M_PI;
	C_vec = Cvec(e0, y0, eta);
	t_vec = tvec(y0, eta, C_vec);
	lam_vec = lamvec(y0, eta, C_vec);
	l_vec = lvec(y0, eta, C_vec);
	y_vec = yvec(y0, eta, C_vec);
	acc_fn_e = gsl_interp_accel_alloc ();
	acc_fw_e = gsl_interp_accel_alloc ();
	acc_y_e = gsl_interp_accel_alloc ();
	acc_e_fn = gsl_interp_accel_alloc ();
	elast = 0;
	F_p = 1./2.*(1+cos(thet)*cos(thet))*cos(2*phi)*cos(2*psi)-cos(thet)*sin(2*phi)*sin(2*psi);
	F_c = 1./2.*(1+cos(thet)*cos(thet))*cos(2*phi)*cos(2*psi - M_PI/4.)-cos(thet)*sin(2*phi)*sin(2*psi - M_PI/4.);
	Q = -(F_p*(1+cos(iot))/2 + 1i*cos(iot)*F_c)*(cos(2*bet) + 1i*sin(2*bet));
	DL = 1;
	phase_container.resize(3);
}

TaylorF2e::TaylorF2e(double e_in, double p_in, double M_in, double eta_in, double psi_in, double phi_in, double thet_in, double iot_in, double bet_in, double f0_in, double fend_in, double df_in) {
	// TODO Auto-generated constructor stub
	msun = 4.925502303934785*pow(10, -6);
	e0 = e_in;
	p0 = p_in;
	y0 = 1/sqrt(p0);
	M = M_in*msun;
	eta = eta_in;
	psi = psi_in*M_PI;
	phi = phi_in*M_PI;
	thet = thet_in*M_PI;
	iot = iot_in*M_PI;
	bet = bet_in*M_PI;
	f0 = f0_in; 				//the f0 to begin the waveform
	fend = fend_in;				//the final f at which to compute waveform
	df = df_in;					//the df of the waveform
	C_vec = Cvec(e0, y0, eta);			//computes the constants of integration coming from y(e0) = y0
	t_vec = tvec(y0, eta, C_vec);		//computes the constants related to t(e)
	lam_vec = lamvec(y0, eta, C_vec);	//computes the constants related to lamda(e)
	l_vec = lvec(y0, eta, C_vec); 		//computes the constants related to l(e)
	y_vec = yvec(y0, eta, C_vec);		//computes the constants related to y(e)
	acc_fn_e = gsl_interp_accel_alloc (); 	//alloc the interpolation accelerators for the interpolations of F_n(e), F_w(e), y(e), and e(F_n)
	acc_fw_e = gsl_interp_accel_alloc ();
	acc_y_e = gsl_interp_accel_alloc ();
	acc_e_fn = gsl_interp_accel_alloc ();
	elast = 0; 								//initialize the value which tracks the e corresponding to choice of waveform truncations
	F_p = 1./2.*(1+cos(thet)*cos(thet))*cos(2*phi)*cos(2*psi)-cos(thet)*sin(2*phi)*sin(2*psi); //antenna functions
	F_c = 1./2.*(1+cos(thet)*cos(thet))*cos(2*phi)*cos(2*psi - M_PI/4.)-cos(thet)*sin(2*phi)*sin(2*psi - M_PI/4.);
	Q = -(F_p*(1+cos(iot))/2 + 1i*cos(iot)*F_c)*(cos(2*bet) + 1i*sin(2*bet)); //overall amplitude factors
	DL = 1;
	over_amp = Q*sqrt(10*M_PI*eta)*M*M/DL;
	calls = 0; 								//trackers for the stationary phase inversion
	count = 0;
	e_stat_last = e0;						//this holds the last value of the stationary phase inversion, which is used for an initial guess in the last
	phase_container.resize(3);				//holds the current value of the phase
}

TaylorF2e::TaylorF2e(double M_in, double eta_in, double e_in, double p_in, double ampre, double ampim, double f0_in, double fend_in, double df_in){
	msun = 4.925502303934785*pow(10, -6);
	e0 = e_in;
	p0 = p_in;
	y0 = 1/sqrt(p0);
	M = M_in*msun;
	eta = eta_in;
	f0 = f0_in; 				//the f0 to begin the waveform
	fend = fend_in;				//the final f at which to compute waveform
	df = df_in;					//the df of the waveform
	C_vec = Cvec(e0, y0, eta);			//computes the constants of integration coming from y(e0) = y0
	t_vec = tvec(y0, eta, C_vec);		//computes the constants related to t(e)
	lam_vec = lamvec(y0, eta, C_vec);	//computes the constants related to lamda(e)
	l_vec = lvec(y0, eta, C_vec); 		//computes the constants related to l(e)
	y_vec = yvec(y0, eta, C_vec);		//computes the constants related to y(e)
	acc_fn_e = gsl_interp_accel_alloc (); 	//alloc the interpolation accelerators for the interpolations of F_n(e), F_w(e), y(e), and e(F_n)
	acc_fw_e = gsl_interp_accel_alloc ();
	acc_y_e = gsl_interp_accel_alloc ();
	acc_e_fn = gsl_interp_accel_alloc ();
	elast = 0; 								//initialize the value which tracks the e corresponding to choice of waveform truncations
	over_amp = (ampre + 1i*ampim)*sqrt(10*M_PI*eta)*M*M;
	e_stat_last = e0;						//this holds the last value of the stationary phase inversion, which is used for an initial guess in the last
	phase_container.resize(3);
}
TaylorF2e::TaylorF2e(double M_in, double eta_in, double e_in, double p_in, double ampmag, double f0_in, double fend_in, double df_in){
	msun = 4.925502303934785*pow(10, -6);
	e0 = e_in;
	p0 = p_in;
	y0 = 1/sqrt(p0);
	eta = eta_in;
	M = M_in*msun/pow(eta, 3./5.);
	f0 = f0_in; 				//the f0 to begin the waveform
	fend = fend_in;				//the final f at which to compute waveform
	df = df_in;					//the df of the waveform
	C_vec = Cvec(e0, y0, eta);			//computes the constants of integration coming from y(e0) = y0
	t_vec = tvec(y0, eta, C_vec);		//computes the constants related to t(e)
	lam_vec = lamvec(y0, eta, C_vec);	//computes the constants related to lamda(e)
	l_vec = lvec(y0, eta, C_vec); 		//computes the constants related to l(e)
	y_vec = yvec(y0, eta, C_vec);		//computes the constants related to y(e)
	acc_fn_e = gsl_interp_accel_alloc (); 	//alloc the interpolation accelerators for the interpolations of F_n(e), F_w(e), y(e), and e(F_n)
	acc_fw_e = gsl_interp_accel_alloc ();
	acc_y_e = gsl_interp_accel_alloc ();
	acc_e_fn = gsl_interp_accel_alloc ();
	elast = 0; 								//initialize the value which tracks the e corresponding to choice of waveform truncations
	over_amp = ampmag*sqrt(10*M_PI*eta)*M*M;
	e_stat_last = e0;						//this holds the last value of the stationary phase inversion, which is used for an initial guess in the last
	phase_container.resize(3);
}
TaylorF2e::TaylorF2e(double M_in, double eta_in, double e_in, double ampmag, double f0_in, double fend_in, double df_in){
	msun = 4.925502303934785*pow(10, -6);
	e0 = e_in;
	eta = eta_in;
	M = M_in*msun/pow(eta, 3./5.);
	y0 = pow(M*M_PI*2*5, 1./3.);
//	y0 = 0.141421;
	f0 = f0_in; 				//the f0 to begin the waveform
	fend = fend_in;				//the final f at which to compute waveform
	df = df_in;					//the df of the waveform
	C_vec = Cvec(e0, y0, eta);			//computes the constants of integration coming from y(e0) = y0
	t_vec = tvec(y0, eta, C_vec);		//computes the constants related to t(e)
	lam_vec = lamvec(y0, eta, C_vec);	//computes the constants related to lamda(e)
	l_vec = lvec(y0, eta, C_vec); 		//computes the constants related to l(e)
	y_vec = yvec(y0, eta, C_vec);		//computes the constants related to y(e)
	acc_fn_e = gsl_interp_accel_alloc (); 	//alloc the interpolation accelerators for the interpolations of F_n(e), F_w(e), y(e), and e(F_n)
	acc_fw_e = gsl_interp_accel_alloc ();
	acc_y_e = gsl_interp_accel_alloc ();
	acc_e_fn = gsl_interp_accel_alloc ();
	elast = 0; 								//initialize the value which tracks the e corresponding to choice of waveform truncations
	over_amp = ampmag*sqrt(10*M_PI*eta)*M*M;
	e_stat_last = e0;						//this holds the last value of the stationary phase inversion, which is used for an initial guess in the last
	phase_container.resize(3);
}

TaylorF2e::~TaylorF2e() {
	// TODO Auto-generated destructor stub
}

double TaylorF2e::get_p_e(double e){
	double y_e =  yevec(e, y_vec);
	return 1./(y_e*y_e);
}

////////////////////////////////////////////////////////////////
//Below we set up interpolations of y(e), F_n(e), F_w(e), and e(F_n) which are used
// to invert the stationary phase condition and in the amplitude
////////////////////////////////////////////////////////////////
void TaylorF2e::init_interps(int N){
	vector<double> y_in (2*N);     				//make this vector bigger just so I don't run into memory issues
	double *y, *e, *F_n, *F_w, *F_n_rev, *e_rev;	//declare the matricies whichi I'll interpolate with GSL
	int i = 0;
//	double ein = 0.9;		//want to start sampling a bit before the initial eccentricity, just to give the interpolation range a little play
//	double de = (double) ein/N;
//	if (e0 < 0.1) {ein = 0.6; de = (double) 0.6/N;}
//	if (e0 < 0.01) {ein = 0.1; de = (double) ein/N;}
//	if (e0 < 0.001) {ein = 0.01; de = (double) ein/N;}
//	if (e0 < 0.0001) {ein = 0.001; de = (double) ein/N;}
//	if (e0 < 0.00001) {ein = 0.0001; de = (double) ein/N;}
//	if (e0 < 0.000001) {ein = 0.00001; de = (double) ein/N;}
//	if (e0 < 0.0000001) {ein = 0.000001; de = (double) ein/N;}
//	if (e0 < 0.00000001) {ein = 0.0000001; de = (double) ein/N;}

	double ein = 0.9;
	double de = abs(11./N);
	double lne = log(ein) - de*i;
	y_in[i] = yevec(exp(lne), y_vec);
//	y_in[i] = yevec(ein - de*i, y_vec);

	while(y_in[i] < 0.35){ 							//sample y(e) until it is larger than 0.4... a regime no waveform should go to anyway
		i++;
//		y_in[i] = yevec(ein - de*i, y_vec);
		lne = log(ein) - de*i;
		y_in[i] = yevec(exp(lne), y_vec);
//		cout << "y val = " << y_in[i] << " e val = " << ein - de*i << " i = " << i << endl;
//		cout << "y val = " << y_in[i] << " e val = " << exp(lne) << " i = " << i << endl;
	}

//	elast = ein - de*i; 							//the last value of e sampled
	elast = exp(lne);
	e = (double*)malloc(sizeof(double) * (i+1));	//allocate appropriate memory for the different matricies to be passsed to GSL interpolators
	e_rev = (double*)malloc(sizeof(double) * (i+1));
	y = (double*)malloc(sizeof(double) * (i+1));
	F_n = (double*)malloc(sizeof(double) * (i+1));
	F_n_rev = (double*)malloc(sizeof(double) * (i+1));
	F_w = (double*)malloc(sizeof(double) * (i+1));

	for(int s = 0; s < i + 1; s++){					//sample y,e,F_w,F_n such that e is increasing (GSL needs the independent variable to be increasing)
		y[s] = y_in[i - s];
		e[s] = exp(log(elast) + de*s);
		F_w[s] = pow(y[s], 3.)*pow(1 - pow(e[s],2),1.5)/(M*2*M_PI);
		F_n[s] = F_w[s]*(1. - 3.*y[s]*y[s] + (-18 + 28*eta - (51 - 26*eta)*pow(e[s],2))/4.*pow(y[s], 4.)
				+ (192 - 896*pow(eta,2) - pow(e[s],4)*(2496 - 1760*eta + 1040*pow(eta,2)) - pow(e[s],2)*(8544 + 5120*pow(eta,2) - eta*(17856 - 123*pow(M_PI,2))) + eta*(14624 - 492*pow(M_PI,2)) -
					     (1920 - 768*eta + (3840 - 1536*eta)*pow(e[s],2))*pow(1 - pow(e[s],2),0.5))/128.*pow(y[s],6.));
//		cout << "y val = " << y[s] << " e val = " << e[s] << " F_w val = " << F_w[s] << " F_n val = " << F_n[s] << endl;
	}
//	cout << "sampled" << endl;

//	cout << "Max value of orbital frequency (n) = " << F_n[0] << endl;

	spline_y_e = gsl_spline_alloc (gsl_interp_cspline, i+1); 		//set up the GSL interpolations
	spline_fn_e = gsl_spline_alloc (gsl_interp_cspline, i+1);
	spline_fw_e = gsl_spline_alloc (gsl_interp_cspline, i+1);
	spline_e_fn = gsl_spline_alloc (gsl_interp_cspline, i+1);
	gsl_spline_init (spline_fn_e, e, F_n, i+1);
	gsl_spline_init (spline_fw_e, e, F_w, i+1);
	gsl_spline_init (spline_y_e, e, y, i+1);

	//need to reverse the order of F_N and e because gsl only interpolates with increasing independent variable
	for(int s = 0; s < i+1; s++){
		F_n_rev[s] = F_n[i - s];
		e_rev[s] = e[i - s];
//		cout << "fn val = " << F_n_rev[s] << " s = " << s << endl;
	}

	gsl_spline_init (spline_e_fn, F_n_rev, e_rev, i+1);

//	cout << "interps set" << endl;
	//dealloc memory
	free(e);
	free(F_n);
	free(F_w);
	free(y);
	free(e_rev);
	free(F_n_rev);
}
//getters for the interpolators
double TaylorF2e::get_fn_e(double e){
	return gsl_spline_eval(spline_fn_e, e, acc_fn_e);
}
double TaylorF2e::get_fw_e(double e){
	return gsl_spline_eval(spline_fw_e, e, acc_fw_e);
}
double TaylorF2e::get_y_e(double e){
	return gsl_spline_eval(spline_y_e, e, acc_y_e);
}
double TaylorF2e::get_e_fn(double fn){
	return gsl_spline_eval(spline_e_fn, fn, acc_e_fn);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// We are going to need to solve for the e_fin which we define to be located when y - (3(1+e))^(-1) < 0 (Section IVB)
// in order to know to what frequencies to evaluate the model
//////////////////////////////////////////////////////////////////////////////////////////////////

double TaylorF2e::fin_cond(double e){
	return get_y_e(e) - 1./(3*(1+e));
}

double TaylorF2e::get_e_fin(){
	// an initial guess
	double e_1 = elast;
	double e_0 = e_1 + e_1/200;
//	cout << "e1 = " << e_1 << " e_0 = " << e_0 << endl;
	// a first iteration of the secant method
	double cond_val_e_1 = fin_cond(e_1);
	double cond_val_e_0 = fin_cond(e_0);
	double e_2 = e_1 - cond_val_e_1*(e_1-e_0)/(cond_val_e_1-cond_val_e_0);
//	cout << "e_2 = " << e_2 << endl;
	while(abs(cond_val_e_1) > 0.01){
		e_0 = e_1;
		cond_val_e_0 = cond_val_e_1;
		e_1 = e_2;
		cond_val_e_1 = fin_cond(e_1);
		e_2 = e_1 - cond_val_e_1*(e_1-e_0)/(cond_val_e_1-cond_val_e_0);
//		cout << "e_2 = " << e_2 << endl;
	}
	return e_2;
}

void TaylorF2e::set_e_fin(){
	e_fin = get_e_fin();
}

////////////////////////////////////////////////////////////
// Code to handle the inversion of the different stationary phase conditions (the _plus and _minus
// and indicies are consistent with the notation of the paper)
////////////////////////////////////////////////////////////

//The inversion of the stationary phase condition for "s" indicies

double TaylorF2e::stat_e_s(double f, int s){
	cout << "evaluating at orbital frequency = " << f/s << endl;
	return get_e_fn(f/s);
}

//For the _plus and _minus indicies these are the functions to be = 0 for the root finder

double TaylorF2e::cond_minus(double& e, double& f, int& j){
	return j*get_fn_e(e) + 2*get_fw_e(e) - f;
}
double TaylorF2e::cond_plus(double e, double f, int j){
	return -j*get_fn_e(e) - 2*get_fw_e(e) - f;
}

//The root finding scheme below inverts the stationary phase conditions via the secant method for the minus index

double TaylorF2e::stat_e_j_minus(double& f, int& j){
	// an initial guess
	calls++;			//to be removed in a final version
	double e_1;
	double e_0;

	e_1 = e_stat_last;			//the initial guess is the last stationary eccentricity (if its the first guess it is e0... implemented elsewhere)
	e_0 = e_1 - e_1/10.;		//the next point for the secant method

	// a first iteration of the secant method
	double cond_val_e_1 = cond_minus(e_1, f, j);
//	cout << "cond val 1 " << endl;
	double cond_val_e_0 = cond_minus(e_0, f, j);
//	cout << "cond val 2 " << endl;
	double e_2 = e_1 - cond_val_e_1*(e_1-e_0)/(cond_val_e_1-cond_val_e_0);

	//iterate method until the condition is satisfied to tolerance
	while(abs(cond_val_e_1) > 0.01){
		e_0 = e_1;
		cond_val_e_0 = cond_val_e_1;
		e_1 = e_2;
//		cout << "e_1 = " << e_1 << endl;
		cond_val_e_1 = cond_minus(e_1,f,j);
		e_2 = e_1 - cond_val_e_1*(e_1-e_0)/(cond_val_e_1-cond_val_e_0);
		count++;
	}
	e_stat_last = e_2;		//store last stationary eccentricity
	return e_2;
}

// a root finding scheme to get e_10hz for a given harmonic j
double TaylorF2e::e_j_10hz(int& j){
	// an initial guess
	double f = 10;
	double e_1;
	double e_0;
//	cout << " Calculating for the j = " << j << " harmonics " << endl;
	e_1 = e0;			//the initial guess is the last stationary eccentricity (if its the first guess it is e0... implemented elsewhere)
	e_0 = e_1 + e_1/10.;		//the next point for the secant method

	// a first iteration of the secant method
	double cond_val_e_1 = cond_minus(e_1, f, j);
//	cout << "cond val 1 " << endl;
	double cond_val_e_0 = cond_minus(e_0, f, j);
//	cout << "cond val 2 " << endl;
	double e_2 = e_1 - cond_val_e_1*(e_1-e_0)/(cond_val_e_1-cond_val_e_0);
	if(e_2 < elast) {e_2 = elast + elast/10;}
	//iterate method until the condition is satisfied to tolerance
	while(abs(cond_val_e_1) > 0.01){
		e_0 = e_1;
		cond_val_e_0 = cond_val_e_1;
		e_1 = e_2;
//		cout << "e_1 = " << e_1 << endl;
		cond_val_e_1 = cond_minus(e_1, f,j);
		e_2 = e_1 - cond_val_e_1*(e_1-e_0)/(cond_val_e_1-cond_val_e_0);
		if(e_2 < elast) {e_2 = elast + elast/10;}
		count++;
	}
	return e_2;
}


/////////////////////////////////////////////////
// Scheme how to sample the different harmonics
// This pretty much relies on drawing a line on figure 6 (now j(e) = 28*e + 2), and not sampling any harmonics at eccentricities below this line
// as well as not sampling any of the harmonics of just l or any with j < -1
////////////////////////////////////////////////

//f(e,j)
double TaylorF2e::fourier_f_minus_e(double e, int j){
	return j*get_fn_e(e) + 2*get_fw_e(e);
}
// e(j), and j(e) at which to cut off the harmonics

// use a line to cut off harmonics

double e_j_cutoff(int j){
	return (j-1)/28.;
}
double j_e_cutoff(double e){
	return 28*e+1;
}

// more conservative

//double e_j_cutoff(int j){
//	return (j)/22.;
//}
//double j_e_cutoff(double e){
//	return 22*e;
//}

//This one appears to be too conservative a scheme

//double e_j_cutoff(int j){
//	if(j-1 <= 0){
//		return 0;
//	} else {
//		return sqrt((double)(j-1)/26.);}
//}
//double j_e_cutoff(double e){
//	return 26*e*e+1;
//}

// convert a frequency to index in final array of {{h_j(f), f}}
double TaylorF2e::raw_freq_to_disc_ind(double f){
	return (f-f0)/df;
}

void TaylorF2e::make_scheme(){
	e_10hz.resize(17);
	for (int i = -1; i < 16; i++){
		e_10hz[i + 1] = e_j_10hz(i);
//		cout << "e 10hz = " << e_10hz[i + 1]  << " j = " << i << endl;
	}
	int k = 1;
	while(j_e_cutoff(e_10hz[k]) > k){
//		cout << "k = " << k << " cutoff val = " << j_e_cutoff(e_10hz[k]) << endl;
		k++;
		if (k == 16){
			break;
		}
	}

//	cout << "k = " << k << endl;


	j_min_max = k-1;				//maximum harmonic to sample
	if (j_min_max > 15) {j_min_max = 15;}			//The highest index j I have computed is j=15

	double e_fin = get_e_fin();						//get the final eccentricity
	double e_hold = 0;								//just holds intermediate values
	j_min_range.resize(j_min_max + 2, vector<int> (2));		//a vector that holds the initial index and last at which to sample any of the j harmonics


	for(int j = -1; j < j_min_max + 1; j++){									//fill this vector
//		cout << " e 10hz = " << e_10hz[j + 1] << endl;
		j_min_range[j + 1][0] = floor(raw_freq_to_disc_ind(fourier_f_minus_e(e_10hz[j + 1], j)));
		if(e_fin > e_j_cutoff(j)){e_hold = e_fin;} else {e_hold = e_j_cutoff(j);}
		if(e_hold > e_10hz[j + 1]) {e_hold = e_10hz[j + 1];}
//		cout << "e hold = " << e_hold << endl;
		j_min_range[j + 1][1] = floor(raw_freq_to_disc_ind(fourier_f_minus_e(e_hold, j)));
//		cout << "j = " << j  << " lower f = " << j_min_range[j + 1][0] << " (" << fourier_f_minus_e(e_10hz[j + 1], j) << ") "<< " higher f = " << j_min_range[j + 1][1] << " (" << fourier_f_minus_e(e_hold, j) << ") " << endl;
	}
}

/////////////////////////////////////////////
// Lookups for the amplitudes
/////////////////////////////////////////////

double TaylorF2e::amplookup_j(double& e, double& y, int& j){
	switch (j){
	case -1: return N_m_1_e_20(y, e, eta);
			 break;
	case 0: return N_0_e_20(y, e, eta);
			 break;
	case 1: return N_p_1_e_20(y, e, eta);
			 break;
	case 2: return N_p_2_e_20(y, e, eta);
			 break;
	case 3: return N_p_3_e_20(y, e, eta);
			 break;
	case 4: return N_p_4_e_20(y, e, eta);
			 break;
	case 5: return N_p_5_e_20(y, e, eta);
			 break;
	case 6: return N_p_6_e_20(y, e, eta);
			 break;
	case 7: return N_p_7_e_20(y, e, eta);
			 break;
	case 8: return N_p_8_e_20(y, e, eta);
			 break;
	case 9: return N_p_9_e_20(y, e, eta);
			 break;
	case 10: return N_p_10_e_20(y, e, eta);
			 break;
	case 11: return N_p_11_e_20(y, e, eta);
			 break;
	case 12: return N_p_12_e_20(y, e, eta);
			 break;
	case 13: return N_p_13_e_20(y, e, eta);
			 break;
	case 14: return N_p_14_e_20(y, e, eta);
			 break;
	case 15: return N_p_15_e_20(y, e, eta);
			 break;
	}
}
double TaylorF2e::amplookup_s(double e, int n){
	if (n == 1){
		return G_1_e_20(e);
	} else if (n == 2){
		return G_2_e_20(e);
	} else if (n == 3){
		return G_3_e_20(e);
	} else if (n == 4){
		return G_4_e_20(e);
	} else if (n == 5){
		return G_5_e_20(e);
	} else if (n == 6){
		return G_6_e_20(e);
	} else if (n == 7){
		return G_7_e_20(e);
	} else if (n == 8){
		return G_8_e_20(e);
	} else if (n == 9){
		return G_9_e_20(e);
	} else if (n == 10){
		return G_10_e_20(e);
	} else if (n == 11){
		return G_11_e_20(e);
	} else if (n == 12){
		return G_12_e_20(e);
	} else if (n == 13){
		return G_13_e_20(e);
	} else if (n == 14){
		return G_14_e_20(e);
	} else if (n == 15){
		return G_15_e_20(e);
	} else {
		cout << "Index error looking up harmonics n" << endl;
		return 0;
	}
}

////////////////////////////////////////////////////////////
// This function samples h_j(f)
////////////////////////////////////////////////////////////


complex<double> TaylorF2e::h_j_minus(double& f, int& j){
	double e = stat_e_j_minus(f, j);						//invert stat phase condition
//	cout << "estat = " << e << endl;
	double y = get_y_e(e);
	double nj = amplookup_j(e, y, j);						//compute amplitude N_j
	complex<double> amp = over_amp*sqrt(1/((j+2)*(96+292*e*e+37*pow(e,4))))*pow(y,-7./2.)*nj;			//compute overall amplitude
	phasevec(phase_container, e, lam_vec, l_vec, t_vec);												//compute phase functions
	double phasefact = 2*M_PI*f*M*phase_container[0] - M_PI/4 - (j*phase_container[2] + 2*phase_container[1]); //compute phase

	return amp*(cos(phasefact) + 1i*sin(phasefact));			// assemble h(f)
}

////////////////////////////////////////////////////////////
// This function generates the waveform using the scheme developed above
////////////////////////////////////////////////////////////

void TaylorF2e::make_F2e_min(){
	int N = (fend-f0)/df + 1;										//amount of frequency samples
	F2_min.resize(j_min_max + 2, vector<complex<double>> (N));		//size up the array which stores the different harmonics
	double f;														//hold values of f (doesn't like to cast otherwise)

	for(int j = -1; j < j_min_max + 1; j++){ //iterate over j
//		cout << "j = " << j << endl;
		e_stat_last = e_10hz[j+1];
		for(int i = j_min_range[j + 1][0]; i < j_min_range[j + 1][1] + 1; i++){ //iterate frequencies which are given by the scheme
			if(i > N - 1) {break;}
//			cout << "f = " << i*df << " i = " << i << " N = " << N <<  endl;
			f = i*df;
			F2_min[j + 1][i] = h_j_minus(f, j);									//store h_j(f)
		}														//reset the last stationary eccentricity each time so that it starts at e0 for each harmonic
	}
//	cout << "done" << endl;
	gsl_spline_free(spline_y_e);
	gsl_spline_free(spline_fn_e);
	gsl_spline_free(spline_fw_e);
	gsl_spline_free(spline_e_fn);
	gsl_interp_accel_free(acc_y_e);
	gsl_interp_accel_free(acc_fn_e);
	gsl_interp_accel_free(acc_fw_e);
	gsl_interp_accel_free(acc_e_fn);
}

vector<vector<complex<double>>> TaylorF2e::get_F2e_min(){
	make_F2e_min();
	return F2_min;
}

//Finish this if we want the full shit

void TaylorF2e::make_F2e_summed(){
	make_F2e_min();
	int N = (fend-f0)/df + 1;
	int jsize = F2_min.size();

	for(int i = 0; i < N; i++){
		for(int j = 0; j < jsize; j++){
			F2_summed[i] += F2_min[j][i]*(cos(j*loff) + 1i*sin(j*loff));
		}
	}
}

} /* namespace std */
