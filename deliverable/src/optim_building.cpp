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

    string init_sim_cond_fname = "../init_conds/init_cond_build.txt";
    if (argc > 1) {
        init_sim_cond_fname = argv[1];
    }

    string xg_data_fname = "../init_conds/elcentro_NS.dat.txt";
    if (argc > 2) {
        xg_data_fname = argv[2];
    }

    string init_optim_cond_fname = "../init_conds/init_cond_optim.txt";
    if (argc > 3) {
        init_optim_cond_fname = argv[3];
    }

    string output_fname = "../resutls/optim_build_mcmc.txt";
    if (argc > 4) {
        output_fname = argv[4];
    }

    // load in parameter values from txt files
    map<string, double> init_optim_map = load_init_txt(init_optim_cond_fname);

    double est_des_err = init_optim_map["est_data_error"];
    double mcmc_points = init_optim_map["n_points"];
    double mcmc_burn = init_optim_map["n_burn"] * mcmc_points;

    vector<double> p_g_old {
        init_optim_map["md_guess"],
        init_optim_map["kd_guess"],
        init_optim_map["cd_guess"]
    };

    vector<double> p_prior {
        init_optim_map["md_prior"],
        init_optim_map["kd_prior"],
        init_optim_map["cd_prior"]
    };

    vector<double> p_sig {
        init_optim_map["md_sig"],
        init_optim_map["kd_sig"],
        init_optim_map["cd_sig"]       
    };

    vector<double> p_d {
        init_optim_map["md_d"],
        init_optim_map["kd_d"],
        init_optim_map["cd_d"]  
    };

    // get the simulation results usign the inital 
    // parameter guess'
    vector<vector<double>> sim_res = sim_building(
        xg_data_fname,
        init_sim_cond_fname,
        p_g_old[0],
        p_g_old[1],
        p_g_old[2]
    );

    int len_data = sim_res.size();
    int top_floor_pos_index = sim_res[0].size() - 4; // get index to top floor position

    vector<double> sim_vals_old(len_data, 0);
    vector<double> des_vals(len_data, 0);

    for (int i=0; i<len_data; i++) {
        sim_vals_old[i] = sim_res[i][top_floor_pos_index]; // save top floor position
    }

    // get the err on the simulation results usign the inital
    // parameter guess'
    double sim_err_old = est_error(des_vals, sim_vals_old, est_des_err, 3);

    // p_prior_old = p_factor * prior_build(p_g_old, p_prior, p_sig);

    // initialize variables for use in mcmc loop
    vector<vector<double>> p_g_hist;
    vector<double> sim_err_hist; // history of the error in the top floor position

    vector<double> p_g_new (p_g_old);

    double p_factor;
    double p_prior_old;
    double p_prior_new;
    double p_accept;

    vector<double> sim_vals_new(sim_vals_old);
    double sim_err_new;

    default_random_engine rand_gen;

    normal_distribution<double> md_norm_dist(0.0,p_d[0]);
    normal_distribution<double> kd_norm_dist(0.0,p_d[1]);
    normal_distribution<double> cd_norm_dist(0.0,p_d[2]);

    uniform_real_distribution<double> r_accept(0, 1);

    ofstream log;
    log.open("../results/mcmc_log.txt");

    // begin mcmc loop
    int iter_count = 0;
    for (int i=0; i<mcmc_points; i++) {
        // log << "\n";
        // log << "Iteration: " << i << "\n";
        // log << "mass old:" << p_g_old[0] << "\n";

        // update optimization parameters
        p_g_new[0] = p_g_old[0] +  md_norm_dist(rand_gen);
        p_g_new[1] = p_g_old[1] +  kd_norm_dist(rand_gen);
        p_g_new[2] = p_g_old[2] +  cd_norm_dist(rand_gen);

        // log << "mass new:" << p_g_new[0] << "\n";
        // log << "P_accept: " << p_accept << "\n";
        // log << "\n";

        // calculate new simulation results and error
        // based on new parameter values

            // model building position values based on the new
            // mass and spring coefficent values of the
            // dampener
        update_top_floor_pos(
            xg_data_fname,
            init_sim_cond_fname,
            p_g_new,
            sim_res,
            sim_vals_new
        );

            // calculate the "sum of squares" value comparing the
            // data and model given our estimate of the errors
            // in the data
        update_error_est(
            des_vals, 
            sim_vals_new,
            est_des_err, 
            3,
            sim_err_new
        );

        // calculate the probability of acceptance for both
        // the mass and spring coeficient parameters, this uses
        // bayesian stats to include the assumption that the
        // paramter space for each parameter is gaussian
        p_factor = exp(-sim_err_new+sim_err_old);
        p_prior_old = p_factor * prior_build(p_g_old, p_prior, p_sig);
        p_prior_new = p_factor * prior_build(p_g_new, p_prior, p_sig);

            // calculate current p_accept
        p_accept = p_prior_new/p_prior_old;

        if (r_accept(rand_gen) < p_accept) {
            p_g_old[0] = p_g_new[0];
            p_g_old[1] = p_g_new[1];
            p_g_old[2] = p_g_new[2];
            sim_err_old = sim_err_new;

            if (iter_count > mcmc_burn) {
                p_g_hist.push_back(p_g_new);
                sim_err_hist.push_back(sim_err_new);
            }
        } 

        log << "\n";
        log << "Iteration: " << i << "\n";
        log << "sim error: " << sim_err_new << "\n";
        log << "mass: " << p_g_old[0] << "\n";
        log << "k: " << p_g_old[1] << "\n";
        log << "c: " << p_g_old[2] << "\n";
        log << "P_accept: " << p_accept << "\n";
        log << "\n";

        iter_count += 1;
    }

    log.close();

    ofstream mcmc_output;
    mcmc_output.open(output_fname);

    for (int i=0; i<p_g_hist.size(); i++) {
        mcmc_output << p_g_hist[i][0] << ", " 
                    << p_g_hist[i][1] << ", " 
                    << p_g_hist[i][2] << ", "
                    << sim_err_hist[i] << "\n";
    }
    mcmc_output.close();
}