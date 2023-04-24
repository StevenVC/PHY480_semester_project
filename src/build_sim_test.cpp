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

    string init_cond_fname = "../init_conds/optim_test.txt";
    if (argc > 1) {
        init_cond_fname = argv[1];
    }

    string seis_data_fname = "../init_conds/elcentro_NS.dat.txt";
    if (argc > 2) {
        seis_data_fname = argv[2];
    }

    string output_fname = "../results/build_sim_res.txt";
    if (argc > 3) {
        output_fname = argv[3];
    }

    vector<double> init_conds {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    double mass = 0.326; // Mkg
    vector<vector<double>> pars {
        {mass, mass, (mass*2)*0.03}, // [metric tonnes]
        {650, 650, 3.7}, // [MN/m]
        {6.2, 6.2, 0.197} // [MN-s/m]
    };

    double h = 0.012;

    vector<double> x_g_data;

    ifstream x_g_file(seis_data_fname);
    string line;

    while (getline(x_g_file, line)) {
        x_g_data.push_back(atof(line.substr(line.find(" ")+1).c_str()));
    }

    x_g_file.close();

    vector<vector<double>> state_hist = rk4_build(init_conds, pars, x_g_data, h);

    ofstream outfile;
    outfile.open(output_fname);

    outfile << "x1, v1, x2, v2, x3, v3" << "\n";

    for (auto x : state_hist) {
        for (auto y : x) {
            outfile << y << ",";
        }
        outfile << "\n";
    }

    outfile.close();
}