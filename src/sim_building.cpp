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

    string output_fname = "../results/sim_build_res.txt";
    if (argc > 3) {
        output_fname = argv[3];
    }

    map<string, double> init_cond_map;
    string init_cond_line;
    string prev_string;
    string init_cond_map_index;

    ifstream init_cond_f(init_cond_fname);

    int count=0;
    while (getline(init_cond_f, init_cond_line)) {
        if (init_cond_line=="") {
            continue;
        }

        switch (init_cond_line[0])
        {
        case 'x':
            prev_string = init_cond_line[0];
            break;
        
        case 'v':
            prev_string = init_cond_line[0];
            break;

        case 'm':
            prev_string = init_cond_line[0];
            break;

        case 'k':
            prev_string = init_cond_line[0];
            break;

        case 'c':
            prev_string = init_cond_line[0];
            break;

        default:
            prev_string = "";
            break;
        }


        // if (count==0) {
        //     init_cond_map_index = init_cond_line;
        //     init_cond_map[init_cond_map_index] = 0.0;
        //     count += 1;
        // }
        // else {
        //     init_cond_map[init_cond_map_index] = stod(init_cond_line);
        //     count = 0;
        // }
        
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