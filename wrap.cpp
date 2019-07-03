#include <iostream>
#include <iomanip>
#include <fstream>
#include "TaylorF2e.hpp"

using namespace std;

extern "C" {
    int generate (complex<double>* hp, complex<double>* hc,
                  double M, double eta, double iota, double eref,
                  double lamc, double lc, double Dl,
                  double fend, double df){

        TaylorF2e F2(M, eta, eref, 0, 0, 0, iota, lamc, lc, 0, Dl, 0, fend, df);
        F2.init_interps(1000); //interpolate needed things with 1000 (ish :) ) points
        F2.make_scheme();     //this is a scheme to sample the different harmonics in an efficient manner
        F2.make_F2e_min_plus_cross();
        vector<complex<double>> plus;
        vector<complex<double>> cross;

        plus = F2.get_F2e_min_plus(); // holds h_plus
        cross = F2.get_F2e_min_cross(); // holds h_cross

        copy(plus.begin(), plus.end(), hp);
        copy(cross.begin(), cross.end(), hc);

        return 0;
    }
}
