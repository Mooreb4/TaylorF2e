#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include "TaylorF2e.hpp"

using namespace std;

int main (int argc, const char * argv[]){
    double M = stod(argv[1]); // Chirp mass in solar masses
    double eta = stod(argv[2]); // symmetric dimensionless mass ratio
    double A = exp(stod(argv[3])); // overall amplitude given in log
    double eref = stod(argv[4]); // the reference eccentricity
    double f0 = 0; //starting frequency
    double fend = 1000; // ending frequency
    double df = 0.015; // frequency resolution
    clock_t t = clock();
    TaylorF2e F2(M, eta, eref, A, f0, fend, df);
    F2.init_interps(1000); //interpolate needed things with 1000 (ish :) ) points
    F2.make_scheme();     //this is a scheme to sample the different harmonics in an efficient manner
    vector<vector<complex<double>>> vect;
    
    
    vect = F2.get_F2e_min();
    // "vect" holds the frequency response e.g. vect[j][i].real() is the real part of the ith frequency sample of the jth harmonic;
    t = clock() - t;
    cout << setprecision(16) << "Time to get waveform is  " << (float) t/CLOCKS_PER_SEC << " seconds" << endl;
    return 0;
}
