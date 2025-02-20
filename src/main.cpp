#include "placeholder.hpp"
#include "poisson1DCart.hpp"
#include "poisson2DCyl.hpp"
#include "specie.hpp"
#include "chemistry.hpp"
#include "simulation.hpp"


#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <fstream>   // For std::ofstream
#include <iostream>  // For std::cerr
#include <chrono>

double epsi =8.85418781762e-12;

double no_epsi = 1;

// Function to perform polynomial fitting
Eigen::VectorXd polynomialFit(const std::vector<double>& x, const std::vector<double>& y, int degree, bool type) {
    int n = x.size();
    Eigen::MatrixXd X(n, degree + 1);
    Eigen::VectorXd Y(n);

    // Construct the Vandermonde matrix and the Y vector
    for (int i = 0; i < n; ++i) {
        Y(i) = y[i];
        for (int j = 0; j <= degree; ++j) {
            if (type){
                X(i, j) = pow(x[i], j);
            }else{
                X(i, j) = pow(log(x[i]), j);
            }   
            
        }
    }

    // Solve for the polynomial coefficients using the normal equation
    Eigen::VectorXd coeffs = (X.transpose() * X).ldlt().solve(X.transpose() * Y);
    return coeffs;
}

int main() {

    int size_r = 1;
    int size_z = 501; // 501

    std::map<int, double> eps_map;
    std::map<int, double> sig_map;

    Eigen::MatrixXd ne = Eigen::MatrixXd::Constant(size_r, size_z, 0);
    //Eigen::MatrixXd ne = Eigen::MatrixXd::Constant(size_r, size_z, 0);
    Eigen::MatrixXd ni = Eigen::MatrixXd::Constant(size_r, size_z, 0);

    //double fronteira[] = {0,0,1,0};
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

    // =======
    // Convection schemes Testing
    // =======
    
    double height = 1e20;
    int startRange = 200; //200   // Starting index for the non-zero range
    int endRange = 300;   //300  // Ending index for the non-zero range

    // Set values within the specified range to the desired height
    for (int i = startRange; i < endRange; ++i) {
        for (int j = 0; j < size_r; j++){
            ne(j, i) = height;
            ni(j, i) = height;
        }   
    }

    // ====================
    //
    //  Polynomial Fitting
    //
    // ====================

    std::ifstream mobfile("../bolsig/N2mob.txt"); // Replace "data.txt" with your file name
    if (!mobfile.is_open()) {
        std::cerr << "Error: Unable to open the file!" << std::endl;
        return 1;
    }

    // Read data from the file
    std::vector<double> x, y;
    double xi, yi;
    while (mobfile >> xi >> yi) {
        x.push_back(xi);
        y.push_back(yi);
    }
    mobfile.close();

    // Perform polynomial fitting
    

    std::ifstream tempfile("../bolsig/N2temp.txt"); // Replace "data.txt" with your file name
    if (!tempfile.is_open()) {
        std::cerr << "Error: Unable to open the file!" << std::endl;
        return 1;
    }

    // Read data from the file
    std::vector<double> x2, y2;
    double xi2, yi2;
    while (tempfile >> xi2 >> yi2) {
        x2.push_back(xi2);
        y2.push_back(yi2);
    }
    tempfile.close();
    //std::cout<<y2[0]<<std::endl;
    int degree = 2; // Adjust the degree of the polynomial as needed
    Eigen::VectorXd mob_coeffs = polynomialFit(x, y, degree,1);
    Eigen::VectorXd temp_coeffs = polynomialFit(x2, y2, degree,0);

    std::cout<<mob_coeffs<<std::endl;
    std::cout<<temp_coeffs<<std::endl;


    // =======================
    //
    //
    // =======================

    Poisson2DCyl testPoisson2D(size_r,size_z,20.0e-6,20.0e-6,eps, sig);

    testPoisson2D.set_pol_coeffs(mob_coeffs, temp_coeffs, degree);

    //testPoisson2D.solve(0,0,0,0, &ZERO2D, fronteira_livre, "zero"); 
    std::ofstream file("rho_data.txt");
    std::ofstream dt_file("dt_data.txt");
    
    testPoisson2D.solve(40e3,0,0,0, ne,ni, fronteira_livre);

    testPoisson2D.solve_Poisson();

    testPoisson2D.write_fields("i");
    testPoisson2D.write_dens(file); 
    //testPoisson2D.solve_Poisson();    

    std::cout<<"test"<<std::endl;
    
    
    if (!file.is_open()) {
        std::cerr << "Error opening file!\n";
        return 1;
    }
    int i = 0;
    //while (testPoisson2D.t <= 0) {     //
    auto start = std::chrono::high_resolution_clock::now();
    while (testPoisson2D.t <= 5e-9) {
        testPoisson2D.push_time(i ,8e-14,dt_file, 0);
        i++;

    }
    auto t1 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> d1 = t1 - start;
    
    std::cout << "Time to Run Code: " << d1.count() << " ms" << std::endl;

    testPoisson2D.write_fields("f");
    testPoisson2D.write_dens(file);

    return 0;
}