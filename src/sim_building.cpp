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
    // load in program arguments
    string init_cond_fname = "../init_conds/init_cond_build.txt";
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

    double md = 0.03;
    if (argc > 4) {
        md = stod(argv[4]);
    }

    double kd = 3.7;
    if (argc > 5) {
        md = stod(argv[5]);
    }

    double cd = 0.197;
    if (argc > 6) {
        md = stod(argv[6]);
    }

    // load in inital condition arguments that define 
    // the inital system 
    map<string, double> init_cond_map;
    string init_cond_line;
    string init_cond_map_index;

    ifstream init_cond_f(init_cond_fname);

    int count=0;
    while (getline(init_cond_f, init_cond_line)) {
        if (init_cond_line=="") {
            continue;
        }

        if (count==0) {
            init_cond_map_index = init_cond_line;
            count += 1;
        }
        else {
            init_cond_map[init_cond_map_index] = stod(init_cond_line);
            count = 0;
        }
    }

    init_cond_f.close();

    // format the inital conditons
    double h;
    vector<double> init_conds;
    vector<vector<double>> pars(3);

    for ( const auto &p : init_cond_map ) {
        switch (p.first[0]) {
        case 'h':
            h=p.second;
            break;

        case 'x':
            init_conds.push_back(p.second);
            break;

        case 'v':
            init_conds.push_back(p.second);
            break;

        case 'm':
            pars[0].push_back(p.second);
            break;
        
        case 'k':
            pars[1].push_back(p.second);
            break;
        
        case 'c':
            pars[2].push_back(p.second);
            break;

        default:
            break;
        }
    }

    // get the number of floors plus 1 for the dampener
    int pars_width = pars[0].size(); 

    // get dampener parameters if supplied from the cmd
    // line
    if (argc > 4) {
        pars[0][pars_width-1] = stod(argv[4]); // set the dampener mass fraction
    }

    if (argc > 5) {
        pars[1][pars_width-1] = stod(argv[5]); // set the dampener k value
    }

    if (argc > 6) {
        pars[2][pars_width-1] = stod(argv[6]); // set the dampener c value
    }

    // set the mass of the dampener so it is the correct
    // fraction of the building mass
    double m_total = 0;
    for (int i=0; i<pars_width-1; i++) {
        m_total += pars[0][i];
    }

    pars[0][pars_width-1] = pars[0][pars_width-1]*m_total;

    // load in the data representing x_g
    vector<double> x_g_data;

    ifstream x_g_file(seis_data_fname);
    string line;

    while (getline(x_g_file, line)) {
        x_g_data.push_back(atof(line.substr(line.find(" ")+1).c_str()));
    }

    x_g_file.close();

    // run the simulation using rk4
    vector<vector<double>> state_hist = rk4_build(init_conds, pars, x_g_data, h);

    // write the results out to a txt file
    ofstream outfile;
    outfile.open(output_fname);

        // assign the column labels
    count = 0;
    for (int i=0; i<pars_width*2; i++) {
        if (i%2 == 0) {
            outfile << "x" << count << ", ";
        }
        else {
            outfile << "v" << count  << ", ";
            count += 1;
        }
    }
    outfile << "xd, vd," << "\n";

    for (auto x : state_hist) {
        for (auto y : x) {
            outfile << y << ",";
        }
        outfile << "\n";
    }

    outfile.close();
}