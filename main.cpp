#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include "TaylorF2e.hpp"

using namespace std;

int main (int argc, const char * argv[]){
    double M = stod(argv[1]); // Chirp mass in solar masses
    double eta = stod(argv[2]); // symmetric dimensionless mass ratio
    double eref = stod(argv[3]); // the reference eccentricity
    double theta = stod(argv[4]); // sky angles
    double phi = stod(argv[5]);
    double psi = stod(argv[6]);
    double iota = stod(argv[7]);
    double tc = stod(argv[8]); //time and phase offsets
    double lc = stod(argv[9]);
    double lamc = stod(argv[10]);
    double Dl = stod(argv[11]); //luminosity distance
    double f0 = 0; //starting frequency
    double fend = 1000; // ending frequency
    double df = 0.015; // frequency resolution
    clock_t t = clock();
    TaylorF2e F2(M, eta, eref, theta, phi, psi, iota, lamc, lc, tc, Dl, f0, fend, df);
    F2.init_interps(1000); //interpolate needed things with 1000 (ish :) ) points
    F2.make_scheme();     //this is a scheme to sample the different harmonics in an efficient manner
    F2.make_F2e_min_plus_cross();
    
    vector<complex<double>> plus;
    vector<complex<double>> cross;
    
    plus = F2.get_F2e_min_plus(); // holds h_plus
    cross = F2.get_F2e_min_cross(); // holds h_cross
    
    t = clock() - t;
    cout << setprecision(16) << "Time to get waveform is  " << (float) t/CLOCKS_PER_SEC << " seconds" << endl;
    return 0;
}
