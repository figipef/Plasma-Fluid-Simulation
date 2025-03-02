//#include "placeholder.hpp"
//#include "poisson1DCart.hpp"
#include "poisson2DCyl.hpp"
#include "specie.hpp"
#include "chemistry.hpp"
#include "simulation.hpp"
#include "poissonsolver2d.hpp"
#include "convection.hpp"


#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <fstream>   // For std::ofstream
#include <iostream>  // For std::cerr
#include <chrono>
#include <sstream>
#include <unordered_map>

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
    std::ifstream input_file("input.txt"); // Open the input file
    if (!input_file) {
        std::cerr << "Error: Could not open the file!" << std::endl;
        return 1;
    }
    
    std::unordered_map<std::string, std::string> values;
    std::string line;
    
    while (std::getline(input_file, line)) {
        std::istringstream iss(line);
        std::string key, value, unit;
        
        if (std::getline(iss, key, '=')) {
            iss >> value >> unit; // Extract value and unit
            values[key] = value;  // Store value in the map
        }
    }
    input_file.close();
    
    // Convert to appropriate types
    double gas_temp = std::stod(values["GAS_TEMP"]);
    double gas_pressure = std::stod(values["GAS_PRESSURE"]);
    double electron_density = std::stod(values["ELECTRON_DENSITY"]);
    double electron_mean_energy = std::stod(values["ELECTRON_MEAN_ENERGY"]);
    double secondary_electron_mean_energy = std::stod(values["SECONDARY_ELECTRON_MEAN_ENERGY"]);
    int grid_size = std::stoi(values["GRID_SIZE"]);
    int grid_init = std::stoi(values["GRID_INIT"]);
    int grid_end = std::stoi(values["GRID_END"]);
    double rel_permitivity = std::stod(values["PERMITIVITY"]);
    double length = std::stod(values["LENGTH"]);
    double left_potential = std::stod(values["LEFT_POTENTIAL"]);
    double right_potential = std::stod(values["RIGHT_POTENTIAL"]);
    std::string potential_type = values["POTENTIAL_TYPE"];
    double frequency = std::stod(values["FREQ"]);
    int number_species = std::stoi(values["NUMBER_SPECIES"]);
    int number_reactions = std::stoi(values["NUMBER_REACTIONS"]);
    
    // Print values to verify
    std::cout << "GAS_TEMP: " << gas_temp << " K" << std::endl;
    std::cout << "GAS_PRESSURE: " << gas_pressure << " TORR" << std::endl;
    std::cout << "ELECTRON_DENSITY: " << electron_density << " M-3" << std::endl;
    std::cout << "ELECTRON_MEAN_ENERGY: " << electron_mean_energy << " EV" << std::endl;
    std::cout << "SECONDARY_ELECTRON_MEAN_ENERGY: " << secondary_electron_mean_energy << " EV" << std::endl;
    std::cout << "GRID_SIZE: " << grid_size << std::endl;
    std::cout << "GRID_INIT: " << grid_init << std::endl;
    std::cout << "GRID_END: " << grid_end << std::endl;
    std::cout << "PERMITIVITY: " << rel_permitivity << std::endl;
    std::cout << "LENGTH: " << length << " M" << std::endl;
    std::cout << "LEFT_POTENTIAL: " << left_potential << " V" << std::endl;
    std::cout << "RIGHT_POTENTIAL: " << right_potential << " V" << std::endl;
    std::cout << "POTENTIAL_TYPE: " << potential_type << std::endl;
    std::cout << "OSCILLATION FREQUENCY: " << frequency << std::endl;
    std::cout << "NUMBER_SPECIES: " << number_species << std::endl;
    std::cout << "NUMBER_REACTIONS: " << number_reactions << std::endl;

    int size_r = 1;
    //int size_z = 501; // 501

    Eigen::MatrixXd ne = Eigen::MatrixXd::Constant(size_r, grid_size, 0);
    Eigen::MatrixXd gas_dens = Eigen::MatrixXd::Constant(size_r, grid_size, 0);
    Eigen::MatrixXd null = Eigen::MatrixXd::Constant(size_r, grid_size, 0);

    //double fronteira[] = {0,0,1,0};
    double fronteira_livre[] = {0,0,1,1}; // zmin, zmax, r0, rmax

    Eigen::VectorXd eps = Eigen::VectorXd::Constant(grid_size, epsi);
    Eigen::VectorXd sig = Eigen::VectorXd::Constant(grid_size, 0);

    double gas_density = (gas_pressure * 133.322368) / (gas_temp * 1.380649e-23); // Calculates the gas density

    double grid_step = length / grid_size;

    for (int i = 0; i < grid_size; i++){
        if (i < grid_init || i >= grid_end){
            eps(i) = epsi * rel_permitivity; // for a e_r|| e_0 || e_r type system
        }
    }

    // Set values within the specified range to the desired height
    for (int i = grid_init; i < grid_end; ++i) {
        for (int j = 0; j < size_r; j++){
            ne(j, i) = electron_density;
            gas_dens(j,i) = gas_density;

        }   
    }

    int react_e[] = {0,0,0,1,1,1,0};
    int react_Ar[] = {0,-1,1,-1,0,1,1};
    int react_Ar_star[] = {0,1,-1,0,-1,-2,-1};
    int react_Ar_plus[] = {0,0,0,1,1,1,0};

    Specie electron("e", -1, 511, react_e, ne);
    Specie argon("Ar", 0, 511, react_Ar, gas_dens);
    Specie argon_star("Ar_star", 0, 511, react_Ar_star, null);
    Specie argon_plus("Ar_plus", 1, 511, react_Ar_plus, null);

    //Specie ion("i", 1, 511, react_e, ni);

    std::vector<Specie> species;
    species.push_back(electron);
    species.push_back(argon);
    species.push_back(argon_star);
    species.push_back(argon_plus);

    // NEW MAIN FROM HERE

    Simulation simul(size_r, grid_size, 20.0e-6, grid_step, "cartesian", eps, species, grid_init, grid_end, electron_mean_energy, secondary_electron_mean_energy);
    std::ofstream file("rho_data.txt");
    
    int a = 0;

    while (simul.get_t() <= 5e-8) {

        PoissonSolver2D solver(left_potential*cos(frequency * simul.get_t()),right_potential,0,0,fronteira_livre,sig,simul);
        std::cout <<left_potential*cos(frequency * simul.get_t())<<std::endl;
        solver.solve();
        std::cout << simul.get_Ez1()<<std::endl;
        simul.write_dens(file);

        //std::cout <<a<<std::endl;
        a++;
        if (a%100 == 0){
            std::cout<< "time "<<simul.get_t()<<std::endl;
        }
        simul.push_time(0);
        solver.solve();
    }
    std::cout<< "iter "<<a<<std::endl;
    simul.write_dens(file);
    
    
    //Poisson2DCyl testPoisson2D(size_r,size_z,20.0e-6,20.0e-6,eps, sig);
    //std::ofstream file2("rho_data2.txt");
    //testPoisson2D.solve(10e3,0,0,0, ne,ni, fronteira_livre);
    //testPoisson2D.solve_Poisson();
    //std::cout <<" RIGHT" <<std::endl;
    //std::cout << testPoisson2D.Ez1<<std::endl;
    //testPoisson2D.write_dens(file2);

    /*
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
    */
    return 0;
}