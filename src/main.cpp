#include "placeholder.hpp"
#include "poisson1DCart.hpp"
#include "poisson2DCyl.hpp"

#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <fstream>   // For std::ofstream
#include <iostream>  // For std::cerr

double epsi = 8.85e-12;

double no_epsi = 1;

double ZERO(double x){
    x=1;
    return 0;
}

double constant(double x){
    x=1;
    return 1 * epsi;
}

double ZERO2D(double r, double z){
    r=1;
    z = 2;
    return 0;
}

double CONST2D(double r, double z){
    r=1;
    z = 2;
    return 0.0001 * epsi;
}

double jan(double r, double z){
    return (6 - 4 * r* r - 4*z*z)*exp(-r*r -z*z);
}


int main() {

    int size_r = 85e-12;
    int size_z = 20;

    std::map<int, double> eps_map;
    std::map<int, double> sig_map;

    Eigen::MatrixXd ne = Eigen::MatrixXd::Constant(size_r, size_z, 1e21);
    //Eigen::MatrixXd ne = Eigen::MatrixXd::Constant(size_r, size_z, 0);
    Eigen::MatrixXd ni = Eigen::MatrixXd::Constant(size_r, size_z, 0);

    double fronteira[] = {0,0,1,0};
    double fronteira_livre[] = {0,0,1,1}; // zmin, zmax, r0, rmax

    sig_map[0] = 0;

    //eps_map[0] = no_epsi;
    //eps_map[66] = no_epsi*2;
    //eps_map[132] = no_epsi*10;

    std::vector<std::pair<int, double>> eps_vec = {{0, epsi}};//, {66, 2.0}, {132, 10.0}};
    //std::vector<std::pair<int, double>> sig_vec = {{49,8e-5}};
    std::vector<std::pair<int, double>> sig_vec = {{0,0}};

    Eigen::VectorXd eps = Eigen::VectorXd::Constant(size_z, epsi);
    Eigen::VectorXd sig = Eigen::VectorXd::Constant(size_z, 0);

    //Poisson1DCart testPoisson1D(12, 1, 0.);

    //testPoisson1D.dirichlet(10, false, false, 0, 10.0, &ZERO);

    //testPoisson1D.dirichlet(5, false, false, 0, 10.0, &constant);

    //testPoisson1D.dirichlet(5, true, false, 0, 10.0, &constant);

    //testPoisson1D.dirichlet(5, true, true, 0, 10, &constant);


    // =======
    // Convection schemes Testing
    // =======
    /*
    double height = 10 * 1.6e19;
    int startRange = 25;   // Starting index for the non-zero range
    int endRange = 75;     // Ending index for the non-zero range

    // Set values within the specified range to the desired height
    for (int i = startRange; i < endRange; ++i) {
        for (int j = 0; j < size_r; j++){
            ne(j, i) = height;
        }
        
    }
    */
    // Create the triangular distribution pattern
    /*
    for (int i = 0; i < size_r; ++i) {
        // Determine the peak value for the current row
        int peak = (size_z - 1) / 2; // For 9 columns, peak is at index 4
        for (int j = 0; j <= peak; ++j) {
            if (j > 20 && j <= peak) {
                ne(i, j) = (j-20)*1.6e19;             // Increase to peak
                ne(i, size_z - j - 1) = (j - 20)*1.6e19; // Mirror to create symmetry
            }
        }
    }   

    std::cout << "Triangular Matrix:\n" << ne << std::endl;
    */
    Poisson2DCyl testPoisson2D(size_r,size_z,0.0001,0.0001,eps, sig);

    //testPoisson2D.solve(0,0,0,0, &ZERO2D, fronteira_livre, "zero");

    testPoisson2D.solve(0,0,0,0, ne,ni, fronteira_livre);
    //testPoisson2D.solve(10,0,0,0, ne,ni, fronteira_livre, "zero");
    testPoisson2D.solve_Poisson();
    testPoisson2D.write_fields("i");
    //testPoisson2D.solve_Poisson();    

    std::cout<<"test"<<std::endl;
    
    std::ofstream file("rho_data.txt");
    if (!file.is_open()) {
        std::cerr << "Error opening file!\n";
        return 1;
    }
    for (double time = 0.0; time <=100; time += 1.0) {

        testPoisson2D.push_time(time,8e-14,file);
    }
    
    testPoisson2D.write_fields("f");

    return 0;
}