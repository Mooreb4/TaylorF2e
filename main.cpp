#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include "TaylorF2e.hpp"

using namespace std;

int main (){
    clock_t t = clock();
    TaylorF2e tstF2(0.2, 80, 20, 0.25, 3./7., 3./7., 3./7., 3./7., 3./7., 0, 500, 0.05);
    // This ^ takes arguements tstF2(e0, p0, M, eta, sky angles type shit..., f0, fend, df)
    
    tstF2.init_interps(100); //interpolate needed things with 100 (ish :) ) points
    tstF2.make_scheme();     //this is a scheme to sample the different harmonics.. you have to do this
    vector<vector<complex<double>>> vect;
    vect = tstF2.get_F2e_min(); //this returns the TaylorF2e sampled
    t = clock() - t;
    cout << setprecision(16) << "Time to get waveform is  " << (float) t/CLOCKS_PER_SEC << " seconds" << endl;
    return 0;
}
