#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <random>

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
            coupled spring-mass system at time t; vector<double> [x1,y1,x2,y2]
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

int main(int argc, char **argv) {

    string output_file = "results.txt";
    if (argc > 1) {
        output_file = argv[1];
    }

    // parameter values
        // masses:
    double m_1 = 1.0;
    if (argc > 2) {
        m_1 = stod(argv[2]);
    }
    
    double m_2 = 1.5;
    if (argc > 3) {
        m_2 = stod(argv[3]);
    }

        // spring constants:
    double k_1 = 8.0;
    double k_2 = 40.0;
    
        // natural spring lengths:
    double L_1 = 0.5;
    double L_2 = 1.0;

        // friction coefficients:
    double fr_1 = 0.2;
    double fr_2 = 0.8;

    // Initial conditions
        // x1 and x2 are the initial displacements; y1 and y2 are the initial velocities
    double x_1 = 0.5;
    double y_1 = 0.0;
    double x_2 = 2.25;
    double y_2 = 0.0;

    // ODE solver parameters
    double abserr = 1.0e-8;
    double relerr = 1.0e-6;
    double start_time = 0.0;
    double stop_time = 10.0;
    double step_size = 0.001;

    if (argc > 4) {
        stop_time = stod(argv[4]);
    } 

    vector<double> state_vars_init {x_1, y_1, x_2, y_2};
    vector<double> spring_pars {m_1, m_2, k_1, k_2, L_1, L_2, fr_1, fr_2};

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