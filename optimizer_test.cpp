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

    string init_cond_fname = "init_conds/optim_test.txt";
    if (argc > 1) {
        init_cond_fname = argv[1];
    }

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
            init_cond_map[init_cond_map_index] = 0.0;
            count += 1;
        }
        else {
            init_cond_map[init_cond_map_index] = stod(init_cond_line);
            count = 0;
        }
    }

    init_cond_f.close();

    // load in object 1 parameters "building"
    double m_1 = init_cond_map["m_1"];
    double k_1 = init_cond_map["k_1"];
    double fr_1 = init_cond_map["fr_1"];
    double L_1 = init_cond_map["L_1"];
    double x_1 = init_cond_map["x_1"];
    double v_1 = init_cond_map["v_1"];

    // load in object 2 parameters "dampener"
    double m_2_prior = init_cond_map["m_2_prior"];
    double m_2_guess_old = init_cond_map["m_2_guess"];
    double k_2_prior = init_cond_map["k_2_prior"];
    double k_2_guess_old = init_cond_map["k_2_guess"];
    double fr_2 = init_cond_map["fr_2"];
    double L_2 = init_cond_map["L_2"];
    double x_2 = init_cond_map["x_2"];
    double v_2 = init_cond_map["v_2"];

    // load in rk4 params
    double start_time = init_cond_map["rk4_start_t"];
    double end_time = init_cond_map["rk4_end_t"];
    double rk4_h = init_cond_map["rk4_h"];

    // load mcmc parameters
        // estimate of the errors in the desired building position
    double est_data_error = 0.2;
    if (argc > 2) {
        est_data_error = stod(argv[2]);
    }    

    double d_m_2 = 0.01;
    if (argc > 3) {
        d_m_2 = stod(argv[3]);
    }
    double d_k_2 = d_m_2; // assuming step size is equal for d_m & d_k

    double n_points = 1E5;
    if (argc > 4) {
        n_points = atoi(argv[4]);
    }

    double n_burn = n_points*0.10;
    if (argc > 5) {
        n_burn = atoi(argv[5]);
    }

    // adjust "dampener" parameters that are being optimize
    // if desired
    if (argc > 6) {
        m_2_prior = stod(argv[6]);
    }
    if (argc > 7) {
        m_2_guess_old = stod(argv[7]);
    }
    if (argc > 8) {
        k_2_prior = stod(argv[8]);
    }
    if (argc > 9) {
        k_2_guess_old = stod(argv[9]);
    }

    vector<double> obj_1 {m_1, k_1, L_1, fr_1};
    vector<double> obj_2 {m_2_guess_old, k_2_guess_old, L_2, fr_2};
    vector<double> init_conds {x_1, v_1, x_2, v_2};
    vector<double> rk4_params {start_time, end_time, rk4_h};

    // generate inital data and error on data
        // get the positions of the building based on the initial
        // conditions
    vector<double> build_pos_sim_old = est_sim_vals(
        obj_1,
        obj_2,
        init_conds,
        rk4_params
    )[0];

        // we want the buildign to have the dampening be as effective as 
        // possible ideally the building does not move
    vector<double> build_pos_want (build_pos_sim_old.size(), 0);

    double build_pos_err_old = est_error(
        build_pos_want, 
        build_pos_sim_old, 
        est_data_error, 
        2
    );

    // initialize variables for use in mcmc loop
    vector<double> m_2_g_hist;
    vector<double> k_2_g_hist;
    vector<double> errors;

    double m_2_guess_new;
    double k_2_guess_new;

    vector<double> build_pos_sim_new;
    double build_pos_err_new;

    default_random_engine rand_gen;
    normal_distribution<double> m_2_norm_dist(0.0,d_m_2);
    normal_distribution<double> k_2_norm_dist(0.0,d_k_2);

    // begin mcmc loop
    int iter_count = 0;
    for (int i=0; i<n_points; i++) {
        // update optimization parameters
        m_2_guess_new = m_2_guess_old * m_2_norm_dist(rand_gen);
        k_2_guess_new = k_2_guess_old * m_2_norm_dist(rand_gen);
        
        // update obj_2 vector
        obj_2[0] = m_2_guess_new;
        obj_2[1] = k_2_guess_new;

        // model building position values based on the new
        // mass and spring coefficent values of the
        // dampener
        build_pos_sim_new = est_sim_vals(
            obj_1,
            obj_2,
            init_conds,
            rk4_params
        )[0];

        // calculate the sum of squares value comparing the
        // data and model given our estimate of the errors
        // in the data
        build_pos_err_new = est_error(
            build_pos_want, 
            build_pos_sim_new, 
            est_data_error, 
            2
        );

        // calculate the probability of acceptance for both
        // the mass and spring coeficient parameters, this uses
        // bayesian stats to include the assumption that the
        // paramter space for each parameter is gaussian
        double p_m_2 = 1;


    }
}