#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <random>

#include "funcs.h"

using namespace std;

int main(int argc, char **argv) {

    string output_file = "../resutls/results.txt";
    if (argc > 1) {
        output_file = argv[1];
    }

    // parameter values
        // masses:
    double m_1 = 10.0;
    if (argc > 2) {
        m_1 = stod(argv[2]);
    }
    
    double m_2 = 1.5;
    if (argc > 3) {
        m_2 = stod(argv[3]);
    }

        // spring constants:
    double k_1 = 8.0;
    if (argc > 4) {
        k_1 = stod(argv[4]);
    }

    double k_2 = 40.0;
    if (argc > 5) {
        k_2 = stod(argv[5]);
    }

        // friction coefficients:
    double fr_1 = 0.2;
        if (argc > 6) {
        fr_1 = stod(argv[6]);
    }
    
    double fr_2 = 0.8;
        if (argc > 7) {
        fr_2 = stod(argv[7]);
    }

        // natural spring lengths:
    double L_1 = 0.5;
    double L_2 = 1.0;

    // Initial conditions
        // x1 and x2 are the initial displacements; y1 and y2 are the initial velocities
    double x_1 = 0.5;
    double v_1 = 0.0;
    double x_2 = 2.25;
    double v_2 = 0.0;

    // ODE solver parameters
    double start_time = 0.0;
    double stop_time = 10.0;
    double step_size = 0.001;

    if (argc > 8) {
        stop_time = stod(argv[4]);
    }

    // setup inital conditions
    vector<double> state_vars_init {x_1, v_1, x_2, v_2};
    vector<double> spring_pars {m_1, m_2, k_1, k_2, L_1, L_2, fr_1, fr_2};

    // run the simulation
    vector<vector<double>> star_var_final = rk4(state_vars_init, spring_pars, start_time, stop_time, step_size);

    // write to file the results of the program
    ofstream myfile;
    myfile.open(output_file);
    for (auto x : star_var_final) {
        for (auto y : x) {
            myfile << y << ", ";
        }
        myfile << "\n";
    }

    myfile.close();
}