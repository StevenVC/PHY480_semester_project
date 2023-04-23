#include <vector>
#include <cmath>
#include <string>

using namespace std;

vector<double> spring_f(vector<double>& state_vars, double t, vector<double>& spring_pars) {
    /*
    Defines the differential equations for the coupled spring-mass system.

    inputs:
        state_vars : vector of the state variables; vector<double> [x1,y1,xv_2,yv_2]

        t :  time

        spring_pars : vector of the spring parameters; vector<double> [m1,m2,k1,k2,L1,L2,b1,b2]
    
    return:
        f : solution of the differential equation for the 
            coupled spring-mass system at time t; vector<double> [x1,y1,x2,y2] -> [vel for mass 1, accel for mass 1, vel for mass 2, accel for mass 2]
    */
    double x_1 = state_vars[0];
    double y_1 = state_vars[1];
    double x_2 = state_vars[2];
    double y_2 = state_vars[3];

    double m_1 = spring_pars[0];
    double m_2 = spring_pars[1];
    double k_1 = spring_pars[2];
    double k_2 = spring_pars[3];
    double L_1 = spring_pars[4];
    double L_2 = spring_pars[5];
    double fr_1 = spring_pars[6];
    double fr_2 = spring_pars[7];

    vector<double> f {
        y_1,
        (-fr_1*y_1 - k_1*(x_1 - L_1) + k_2*(x_2 - x_1 - L_2)) / m_1,
        y_2,
        (-fr_2*y_2 - k_2*(x_2- x_1 - L_2)) / m_2
    };

    return f;
}

vector<double> build_f(vector<double>& state_vars, vector<vector<double>>& pars, double x_g_t) {
    /*
    Defines the differential equations for the building system

    inputs:
        state_vars : vector of the state variables; vector<double> [x1, v1, x2, v2, x3, v3, ...] -> [pos of floor 1, vel of foor 1, ...]

        pars : vector of the building parameters; vector<vector<double>> [[m1, m2, m3, ...], [k1, k2, k3, ...], [c1, c2, c3, ...]]

        x_g_t : current ground acceleration; double
    
    returns:
        f : solution of the differential equation for the 
            building system at time t; vector<double> [v1, a1, v2, a2, v3, a3, ...]

    */

    int n_v = state_vars.size();

    vector<double> f (n_v, 0);

    f[0] = state_vars[1];
    f[1] = -(state_vars[0]*(pars[1][0] + pars[1][1])/pars[0][0] - state_vars[2]*pars[1][1]/pars[0][0]) - \
            (state_vars[1]*(pars[2][0] + pars[2][1])/pars[0][0] - state_vars[3]*pars[2][1]/pars[0][0]) - \
            x_g_t;

    f[-2] = state_vars[-1];
    f[n_v-1] = -(-state_vars[n_v-4]*(pars[1][n_v-1])/pars[0][n_v-1] + state_vars[n_v-2]*pars[1][n_v-1]/pars[0][n_v-1]) - \
                (-state_vars[n_v-3]*(pars[2][n_v-1])/pars[0][n_v-1] + state_vars[n_v-1]*pars[2][n_v-1]/pars[0][n_v-1]) -  \
                x_g_t;

    double temp_1;
    double temp_2;

    for (int i=2; i<(n_v-2); i++) {
        // check if updating velocity of position
        if (i%2 == 0) {
            f[i] = state_vars[i+1];
        }
        else {
            
            // adjust to get the positions from state_vars and compute with "k" parameter
            temp_1 = -state_vars[i-3]*(pars[1][i-1]/pars[0][i-1]) + \
                      state_vars[i-1]*(pars[1][i-1] + pars[1][i])/pars[0][i-1] - \
                      state_vars[i+1]*(pars[1][i]/pars[0][1]);

            // adjust to get the velocities from state_vars and compute with "c" parameter
            temp_2 = -state_vars[i-2]*(pars[2][i-1]/pars[0][i-1]) + \
                      state_vars[i]*(pars[2][i-1] + pars[2][i])/pars[0][i-1] - \
                      state_vars[i+2]*(pars[2][i]/pars[0][1]);

            f[i] = -(temp_1) - (temp_2) - x_g_t;
        }
    }

    return f;
}

vector<double> vec_add(vector<double>& a, vector<double>& b) {
    /*
    compute the sum of vector a & b
    */
    int vec_size = a.size();

    vector<double> c(vec_size, 0);

    for (int i=0; i<vec_size; i++) {
        c[i] = a[i] + b[i];
    }

    return c;
}

vector<double> vec_mult_scal(vector<double>& a, double b) {
    /*
    compute the sum of vector a and scalar b
    */
    int vec_size = a.size();

    vector<double> c(vec_size, 0);

    for (int i=0; i<vec_size; i++) {
        c[i] = a[i] * b;
    }

    return c;
}

vector<vector<double>> rk4(vector<double> state_vars, vector<double> spring_pars, double a, double b, double h) {
    /*
    rk4 algorithm

    inputs:
        state_vars : vector of the state variables; vector<double> [x1,y1,xv_2,yv_2]

        spring_pars : vector of the spring parameters; vector<double> [m1,m2,k1,k2,L1,L2,b1,b2]

        a : initial time; double

        b : final time; double

        h : step size; double

    returns:
        p_hist: 
    */
    
    // initialize the time array [a,b] with stepsize h
    int n_times = (b-a)/h;
    vector<double> times(n_times, 0);
    for (double i=0; i<n_times; i++) {
        times[i] = i+h;
    }

    int n_vars = state_vars.size();

    vector<vector<double>> p_hist(n_times, vector<double>(n_vars, 0));

    vector<double> k_1;
    vector<double> k_2;
    vector<double> k_3;
    vector<double> k_4;

    vector<double> temp;
    vector<double> temp_o;

    p_hist[0] = state_vars;

    for (int i=0; i<n_times-1; i++) {
        // calc k_1
        k_1 = spring_f(p_hist[i], times[i], spring_pars);

        // calc k_2
        temp = vec_mult_scal(k_1, h/2);
        temp = vec_add(p_hist[i], temp);

        k_2 = spring_f(temp, times[i], spring_pars);

        // calc k_3
        temp = vec_mult_scal(k_2, h/2);
        temp = vec_add(p_hist[i], temp);

        k_3 = spring_f(temp, times[i], spring_pars);

        // calc k_4
        temp = vec_mult_scal(k_3, h);
        temp = vec_add(p_hist[i], temp);

        k_4 = spring_f(temp, times[i], spring_pars);

        // calc x+1
        temp_o = vec_mult_scal(k_2, 2);
        temp = vec_add(k_1, temp_o);

        temp_o = vec_mult_scal(k_3, 2);
        temp = vec_add(temp, temp_o);

        temp = vec_add(k_4, temp);

        temp = vec_mult_scal(temp, h/6);

        p_hist[i+1] = vec_add(p_hist[i], temp);
    }

    return p_hist;
}

vector<vector<double>> rk4_build(vector<double> state_vars, vector<vector<double>> pars, vector<double> t, double h) {
    /*
    rk4 algorithm

    inputs:
        state_vars : vector of the state variables; vector<double> [x1,y1,xv_2,yv_2]

        spring_pars : vector of the spring parameters; vector<double> [m1,m2,k1,k2,L1,L2,b1,b2]

        t : vector contaning timing data; vector<double> [...]
        
    returns:
        p_hist: 
    */
    
    // initialize the time array [a,b] with stepsize h

    int n_times = t.size();
    int n_vars = state_vars.size();

    vector<vector<double>> p_hist(n_times, vector<double>(n_vars, 0));

    vector<double> k_1;
    vector<double> k_2;
    vector<double> k_3;
    vector<double> k_4;

    vector<double> temp;
    vector<double> temp_o;

    p_hist[0] = state_vars;

    for (int i=0; i<n_times-1; i++) {
        // calc k_1
        k_1 = build_f(p_hist[i], pars, t[i]);

        // calc k_2
        temp = vec_mult_scal(k_1, h/2);
        temp = vec_add(p_hist[i], temp);

        k_2 = build_f(temp, pars, t[i]);

        // calc k_3
        temp = vec_mult_scal(k_2, h/2);
        temp = vec_add(p_hist[i], temp);

        k_3 = build_f(temp, pars, t[i]);

        // calc k_4
        temp = vec_mult_scal(k_3, h);
        temp = vec_add(p_hist[i], temp);

        k_4 = build_f(temp, pars, t[i]);

        // calc x+1
        temp_o = vec_mult_scal(k_2, 2);
        temp = vec_add(k_1, temp_o);

        temp_o = vec_mult_scal(k_3, 2);
        temp = vec_add(temp, temp_o);

        temp = vec_add(k_4, temp);

        temp = vec_mult_scal(temp, h/6);
        
        p_hist[i+1] = vec_add(p_hist[i], temp);
    }

    return p_hist;
}

/*
The following for parameter optimization using mcmc
*/
double est_error(vector<double>& y_real, vector<double>& y_est, double data_error, int n_param) {
    int y_size=y_real.size();

    double error = 0;
    
    for (int i=0; i<y_size; i++) {
        error += pow(y_real[i] - y_est[i], 2) / (n_param*pow(data_error,2));
    }

    error = error/(y_size-n_param);

    return error;
}

vector<vector<double>> est_sim_vals(
                  vector<double>& obj_1, 
                  vector<double>& obj_2,
                  vector<double>& init_conds,
                  vector<double>& rk4_params) {
    /*
    Runs the coupled spring simulation
    
    inputs:
        obj_1: vector of parameters describing the "building";
               {m, k, L, fr}
        
        obj_2: vector of parameters describing the "dampener";
               {m, k, L, fr}
        
        init_conds: vector describing the inital position and 
                    velocity of obj_1 & obj_2;
                    {x_1, v_1, x_2, v_2}

        rk4_params: vector of parameters for the rk4 solver;
                    {start_time, stop_time, step_size}

    returns:
        obj_pos: 2d vector containing the position history for 
                 the "building" and the "dampner";
                 {{x_1_i},{x_2_i}}
    */

    vector<double> spring_pars {obj_1[0], obj_2[0], obj_1[1], obj_2[1],
                                obj_1[2], obj_2[2], obj_1[3], obj_2[3]};

    vector<vector<double>> state_hist = rk4(init_conds, spring_pars, rk4_params[0], rk4_params[1], rk4_params[2]);
    
    vector<vector<double>> obj_pos (2, vector<double> (state_hist.size(), 0));

    for (int i=0; i<state_hist.size(); i++) {
        obj_pos[0][i] = state_hist[i][0];
        obj_pos[1][i] = state_hist[i][2];
    }

    return obj_pos;
}

double prior(double m_2, double k_2, 
             double m_2_g, double k_2_g, 
             double m_2_sig, double k_2_sig) {
    /*
    Given values of m_2 & k_2, return a probability based on prior knoweldge
    that the user supplies in (m_2_g, k_2_g, m_2_sig, k_2_sig). The parameter
    PDF's are assumed to be gaussian
    */

    // double p_m_2 = pow(pow(2.0*M_PI*m_2_sig, 2),-0.5) * exp(-pow(m_2-m_2_g,2)/(2.0*pow(m_2_sig,2)));
        
    // double p_k_2 = pow(pow(2.0*M_PI*k_2_sig, 2),-0.5) * exp(-pow(k_2-k_2_g,2)/(2.0*pow(k_2_sig,2)));

    double p_m_2 = exp((2.0 * pow(m_2, 2) + pow(m_2_g, 2) - (2.0 * m_2 * m_2_g)) / (2 * pow(m_2_sig, 2)));

    double p_k_2 = exp((2.0 * pow(k_2, 2) + pow(k_2_g, 2) - (2.0 * k_2 * k_2_g)) / (2 * pow(k_2_sig, 2)));
    
    double p_prior = p_m_2 * p_k_2;

    return p_prior;
}
