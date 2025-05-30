//#include "placeholder.hpp"
//#include "poisson1DCart.hpp"
//#include "poisson2DCyl.hpp"
#define _USE_MATH_DEFINES

#include "specie.hpp"
#include "chemistry.hpp"
#include "simulation.hpp"
#include "poissonsolver2d.hpp"
#include "convection.hpp"

#include <memory>
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

double epsi = 8.85418781762e-12;

double charge = 1.6e-19;

double no_epsi = 1;

std::vector<double> computeCellCenters(const std::vector<double>& distances, const std::vector<int>& cellCounts) {
    std::vector<double> centers;

    if (distances.size() != cellCounts.size()) {
        throw std::invalid_argument("distances and cellCounts must have the same size");
    }

    double x_start = 0.0;
    for (size_t i = 0; i < distances.size(); ++i) {
        double dx = distances[i] / cellCounts[i];
        for (int j = 0; j < cellCounts[i]; ++j) {
            double x_center = x_start + (j + 0.5) * dx;
            centers.push_back(x_center);
        }
        x_start += distances[i];
    }

    return centers;
}

std::vector<double> computeCellEdges(const std::vector<double>& centers) {
    std::vector<double> edges;
    size_t N = centers.size();

    if (N == 0) return edges;

    edges.reserve(N + 1);

    // First edge extrapolated on the left
    edges.push_back(centers[0] - 0.5 * (centers[1] - centers[0]));

    // Interior edges: average of neighboring centers
    for (size_t i = 0; i < N - 1; ++i) {
        edges.push_back(0.5 * (centers[i] + centers[i + 1]));
    }

    // Last edge extrapolated on the right
    edges.push_back(centers[N - 1] + 0.5 * (centers[N - 1] - centers[N - 2]));

    return edges;
}

std::vector<double> computeCellSizes(const std::vector<double>& edges) {

    std::vector<double> cellsizes;
    size_t N = edges.size();

    for (size_t i = 0; i < N - 1; ++i) {
        cellsizes.push_back(edges[i+1] - edges[i]);
    }

    return cellsizes;
}

std::vector<double> computeEastPoissonCoefficients(const std::vector<double>& permitivity, const std::vector<double>& centers, const std::vector<double>& edges) {
    
    // Function to calculate the East Poisson coeffecients

    std::vector<double> coeffecients;
    size_t N = permitivity.size();
    
    if (N == 0) return coeffecients;
    
    coeffecients.reserve(N);
    
    for (size_t i = 0; i < N - 1; ++i) {

        // Calculate the east coefficient for all interior cells, except the last one
        double value = - 1.0 * (permitivity[i] * permitivity[i+1])/((permitivity[i+1]*(edges[i+1] - centers[i]))+(permitivity[i] * (centers[i+1] - edges[i+1])));

        coeffecients.push_back(value);
    }
     
    // The last value is for the ghost cell, so it is used the previous obtained coeffecient
    coeffecients.push_back(coeffecients[N-2]);

    return coeffecients;
}   

std::vector<double> computeWestPoissonCoefficients(const std::vector<double>& east_coeffs){
    // Function to calculate the West Poisson coeffecients

    std::vector<double> west_coeffs;
    size_t N = east_coeffs.size();

    // The first value is for the ghost cell, so it is used the next obtained coeffecient 
    west_coeffs.push_back(east_coeffs[0]);

    for (size_t i = 0; i < N-1; ++i) {

        west_coeffs.push_back(east_coeffs[i]);
    }

    return west_coeffs;
}   

std::vector<double> computeCenterPoissonCoefficients(const std::vector<double>& east_coeffs, const std::vector<double>& west_coeffs){
    // Calculates the Center Poisson Coefficients

    std::vector<double> center_coeffs;
    size_t N = east_coeffs.size();


    for (size_t i = 0; i < N; ++i) {

        center_coeffs.push_back(-1.0*(east_coeffs[i] + west_coeffs[i]));
    }

    return center_coeffs;
}

std::vector<int> computeSigmaMask(const std::vector<int>& positions, const int& N_grid){

    // Return a Mask with 0's and 1's for the Interfaces with surface charges
    size_t N = positions.size();
    std::vector<int> sigma_mask(N_grid - 1, 0);

    for (size_t i = 0; i < N; ++i) {

        sigma_mask[positions[i]] = 1;
    } 

    return sigma_mask;   
}

std::vector<double> computeSigmaCoefficients(const std::vector<int>& sigma_mask, const std::vector<double>& permitivity, const std::vector<double>& centers, const std::vector<double>& edges){
    // Calculates the Right Hand Side Surface Charge coeffecients (before and after the surface)

    std::vector<double> sigma_coeffs;
    size_t N = permitivity.size();


    for (size_t i = 0; i < N; ++i) {

        if (sigma_mask[i] && i < N-1){
            double value = permitivity[i] * (centers[i+1] - edges[i + 1]) / ( permitivity[i] * (centers[i+1] - edges[i + 1]) + permitivity[i+1] * (edges[i+1] - centers[i]));
            sigma_coeffs.push_back(value);
        }
        if (sigma_mask[i-1] && i > 0){
            double value = permitivity[i] * (edges[i] - centers[i]) / ( permitivity[i - 1] * (centers[i] - edges[i]) + permitivity[i] * (edges[i] - centers[i - 1]));
            sigma_coeffs.push_back(value);
        }

        if (i > 0 && !sigma_mask[i] && !sigma_mask[i-1]){
            sigma_coeffs.push_back(0);
        } else if ( i == 0  && !sigma_mask[i]){
            sigma_coeffs.push_back(0);
        }

    }

    return sigma_coeffs;
}

double calculate_g_c(double g_dc, double g_cu, double u, double dx, double dt){
    
    if (std::abs(g_dc - g_cu) < 0.6*std::abs(g_dc + g_cu)){
        //std::cout <<"a"<<std::endl;
        //double xmf = std::copysign(1.0, u) * (dx - std::abs(u) * dt)/2.0;
        return g_dc - (dx*copysign(1.0,-u) + std::abs(u)*dt)/(1.5*std::copysign(1.0,u))*(g_dc - g_cu)/(2.0 *dx*copysign(1.0,u));
        //return g_dc - 4/3*(std::copysign(1.0, u)*dx - xmf)*(g_dc - g_cu)/(2.0 *dx* std::copysign(1.0, u));

    } else if (g_dc * g_cu > 0.0 ){
        //std::cout <<"b"<<std::endl;
        return std::copysign(1.0, g_dc)*2.0*std::min(std::abs(g_dc),std::abs(g_cu));
    
    }else{
        //std::cout <<"c"<<std::endl;
        return std::copysign(1.0, g_dc)*std::min(std::abs(g_dc),std::abs(g_cu));
    }
    
    //return std::copysign(1.0, g_dc)*2*(std::abs(g_dc * g_cu)/(std::abs(g_dc) + std::abs(g_cu)+ 1e-308));
}

std::vector<double> computeFluxUNO3(const std::vector<double>& centers, const std::vector<double>& density, const std::vector<double>& e_field, const std::vector<double>& mob_coef, const std::vector<double>& dif_coef,  const int& plasma_init, const int& plasma_end, const double& dt, const int& charge){
    std::vector<double> fluxes;

    size_t N = centers.size();

    for (size_t i = 0; i < N - 1; ++i) {

        // Fill the flux vector for the boundaries and outside them
        if (i < plasma_init || i >= plasma_end - 1){
            fluxes.push_back(0);

        } else {

            double flux = 0;

            // density derivatives respective to the donator cell (center)
            double g_cu = 0;
            double g_dc = 0;
            double g_c = 0; 

            double n_surface = 0; // density at the cell boundary (i | i + 1)

            double v = charge * e_field[i]; // velocity in the directional sense, still needs to multiply by the mobility
            double dx = centers[i+1] - centers[i]; // common x difference

            if (v > 0 ){ // celula i -> i+1

                v = v * mob_coef[i]; // true sense of velocity

                g_dc = (density[i+1] - density[i])/(dx); // downstream - center cells

                if (i-1 < 0){
                    g_cu = 0; // No density difference between wall and first cell
                } else {
                    g_cu = (density[i] - density[i - 1])/(centers[i] - centers[i - 1]); // center - upstream cells
                }
                

                g_c = calculate_g_c(g_dc, g_cu, v, dx, dt);

                n_surface = density[i] + 0.5 * std::copysign(1.0, v) * (dx - std::abs(v) * dt) * g_c;  

            } else if (v < 0) { // célula i <- i+1

                v = v * mob_coef[i+1];

                g_dc = -1.0 * (density[i] - density[i + 1])/(dx);

                if (i > N - 3){

                    g_cu = 0;

                } else {

                    g_cu = -1.0 *(density[i + 1] - density[i + 2])/(centers[i+2] - centers[i+1]);
                }
                

                g_c = calculate_g_c(g_dc, g_cu, v, dx, dt);

                n_surface = density[i + 1] + 0.5 * std::copysign(1.0, v) * (dx - std::abs(v) * dt) * g_c;

            }

            flux += v * n_surface;

            double dens_div = (density[i+1] - density[i]) / (dx);

            if (dens_div < 0){
                flux -= dif_coef[i] * dens_div;
            } else if (dens_div > 0){
                flux -= dif_coef[i + 1] * dens_div;
            }

            fluxes.push_back(flux);
        }
    }   

    return fluxes;
}

int main() {

    // Relative permitivity (For DBD case)
    double outer_rel_permitivity = 5;
    double inner_rel_permitivity = 1;

    // Distances between sections SI (m)
    double d1 = 1e-4;
    double d2 = 1e-4;
    double d3 = 1e-4;

    // Number of cells in sections
    int c1 = 5;
    int c2 = 5;
    int c3 = 5;
    
    std::vector<double> distances = {d1, d2, d3}; // Distances between sections SI (m)
    
    std::vector<int> cellCounts = {c1, c2, c3}; // Number of cells in sections

    std::vector<double> centers = computeCellCenters(distances, cellCounts);

    int N = 0; // Number of internal Cells

    for (size_t i = 0; i < cellCounts.size(); ++i) {
        N += cellCounts[i];
    }

    std::cout << "Number of cells "<< N <<"\n";

    // Print the first few centers for verification
    for (size_t i = 0; i < centers.size(); ++i) {
        std::cout << "x[" << i << "] = " << centers[i] << "\n";
    }
    std::vector<double> edges = computeCellEdges(centers);

    // Print the first few edges for verification
    for (size_t i = 0; i < edges.size(); ++i) {
        std::cout << "edge[" << i << "] = " << edges[i] << "\n";
    }

    std::vector<double> cellsizes = computeCellSizes(edges);

    // Print the first few edges for verification
    for (size_t i = 0; i < cellsizes.size(); ++i) {
        std::cout << "cellsizes[" << i << "] = " << cellsizes[i] << "\n";
    }

    std::vector<double> permitivity(N);

    for (int i = 0; i < N; i++){
        if (i < cellCounts[0] || i >= N - cellCounts[2]){
            permitivity[i] = epsi * outer_rel_permitivity; // for a e_r|| e_0 || e_r type system
        } else {
            permitivity[i] = epsi * inner_rel_permitivity;
        }
    }

    for (size_t i = 0; i < N; ++i) {
        std::cout << "Permitivity[" << i << "] = " << permitivity[i] << "\n";
    }

    std::vector<double> east_poisson_coeffecients = computeEastPoissonCoefficients(permitivity,centers, edges);


    for (size_t i = 0; i < N; ++i) {
        std::cout << "East Poisson Coeffecient[" << i << "] = " << east_poisson_coeffecients[i] << "\n";
    }

    std::vector<double> west_poisson_coeffecients = computeWestPoissonCoefficients(east_poisson_coeffecients);


    for (size_t i = 0; i < N; ++i) {
        std::cout << "West Poisson Coeffecient[" << i << "] = " << west_poisson_coeffecients[i] << "\n";
    }

    std::vector<double> center_poisson_coeffecients = computeCenterPoissonCoefficients(east_poisson_coeffecients, west_poisson_coeffecients);


    for (size_t i = 0; i < N; ++i) {
        std::cout << "Center Poisson Coeffecient[" << i << "] = " << center_poisson_coeffecients[i] << "\n";
    }

    // Positions for the charged surfaces
    // Positions are +1, so if I want to place a suface between cell 5 and 6, the index will be 4!!!
    int p1 = 4;
    int p2 = 9;

    std::vector<int> surface_positions = {p1,p2};

    std::vector<int> sigma_mask = computeSigmaMask(surface_positions, N); 

    for (size_t i = 0; i < N - 1; ++i) {
        std::cout << "Surface[" << i << "] = " << sigma_mask[i] << "\n";
    }

    std::vector<double> sigma_coeffs = computeSigmaCoefficients(sigma_mask, permitivity, centers, edges);

    for (size_t i = 0; i < N ; ++i) {
        std::cout << "Sigma Coeff[" << i << "] = " << sigma_coeffs[i] << "\n";
    }

    



    
    return 0;
}

/*
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
    std::cout <<"gas dens: " <<gas_density<<"\n";
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

    //int react_e[] = {0,0,0,1,1,1,0}; // Real reactions
    //int react_Ar[] = {0,-1,1,-1,0,1,1};
    //int react_Ar_star[] = {0,1,-1,0,-1,-2,-1};
    //int react_Ar_plus[] = {0,0,0,1,1,1,0};

    int react_e[] =       {0,0,1,1,1,0}; // Reactions being considered rn are 2,3,4,5,6,7
    int react_Ar[] =      {-1,1,-1,0,1,1};
    int react_Ar_star[] = {1,-1,0,-1,-2,-1};
    int react_Ar_plus[] = {0,0,1,1,1,0};

    Specie electron("e", -1, 9.1093837e-31, react_e, ne);
    Specie argon("Ar", 0, 6.6335209e-26, react_Ar, gas_dens);
    Specie argon_star("Ar_star", 0, 6.6335209e-26, react_Ar_star, null);
    Specie argon_plus("Ar_plus", 1, 6.6335209e-26, react_Ar_plus, null);

    Specie r1_reag[2] = {electron, argon};
    Specie r1_prod[2] = {electron, argon_star};

    Chemistry r1 = Chemistry(2,2, r1_reag, r1_prod, 0, 11.5, 3.55613441e-14, 1.70881680e-02, 9.42791699e-01); // Excitation

    Specie r2_reag[2] = {electron, argon_star};
    Specie r2_prod[2] = {electron, argon};

    Chemistry r2 = Chemistry(2,2, r2_reag, r2_prod, 0, -11.5, 5.00178366e-14, 2.31631264e-02, 1.63097482e+00); // Superelastic

    Specie r3_reag[2] = {electron, argon};
    Specie r3_prod[2] = {electron, argon_plus};

    Chemistry r3 = Chemistry(2,2, r3_reag, r3_prod, 0, 15.8, 1.03817139e-13, 1.04243404e-02, 9.48850558e-01); // Ionization

    Specie r4_reag[2] = {electron, argon_star};
    Specie r4_prod[2] = {electron, argon_plus};

    Chemistry r4 = Chemistry(2,2, r4_reag, r4_prod, 0, 4.3, 1.16698343e-13, 9.08521436e-03, 1.21470021e+00); // StepWize

    double constant_rate_r5 = 3.4e8 / 6.022e-23;

    Specie r5_reag[2] = {argon_star, argon_star};
    Specie r5_prod[3] = {electron, argon, argon_plus};

    Chemistry r5 = Chemistry(2,2, r5_reag, r5_prod, 1, 0, constant_rate_r5, 0, 0);

    double constant_rate_r6 = 1807 / 6.022e-23;

    Specie r6_reag[2] = {argon_star, argon};
    Specie r6_prod[2] = {argon, argon};

    Chemistry r6 = Chemistry(2,2, r6_reag, r6_prod, 1, 0, constant_rate_r6, 0, 0);

    //Specie ion("i", 1, 511, react_e, ni);

    // ALWAYS ADD THE ELECTRONS FIRST
    std::vector<Specie> species;
    species.push_back(electron);
    species.push_back(argon);
    species.push_back(argon_star);
    species.push_back(argon_plus);

    std::vector<Chemistry> chemistries; // Stepwise Ionization
    chemistries.push_back(r1);
    chemistries.push_back(r2);
    chemistries.push_back(r3);
    chemistries.push_back(r4);
    chemistries.push_back(r5);
    chemistries.push_back(r6);

    // NEW MAIN FROM HERE

    Simulation simul(size_r, grid_size, 20.0e0, grid_step, "cartesian", eps, species, chemistries, grid_init, grid_end, electron_mean_energy, secondary_electron_mean_energy, gas_temp, gas_density, gas_pressure);
    
    std::ofstream file_e_dens("../output/e_dens.txt");
    std::ofstream file_ar_dens("../output/ar_dens.txt");
    std::ofstream file_arplus_dens("../output/arplus_dens.txt");
    std::ofstream file_arstar_dens("../output/arstar_dens.txt");
    std::ofstream file_time_steps("../output/time_steps.txt");
    std::ofstream file_current_dens("../output/current_dens.txt");
    std::ofstream file_e_field("../output/e_field.txt");
    std::ofstream file_e_energy("../output/e_energy.txt");
    std::ofstream file_bflux("../output/bflux.txt");
    std::ofstream file_pot("../output/pot.txt");

    int a = 0;

    Eigen::MatrixXd electron_fluxes = Eigen::MatrixXd::Zero(size_r, grid_size); // To be used by electron energy after an iteration

    PoissonSolver2D solver(fronteira_livre, sig, simul);

    auto start = std::chrono::high_resolution_clock::now();

    while (simul.get_t() <= 100e-6) {
    //while (a < 1) {    
        //if (a%1000 == 0){
        //    solver = std::make_unique<PoissonSolver2D>(fronteira_livre, sig, simul);
        //}

        //PoissonSolver2D solver(left_potential*sin(2*3.1416*frequency * simul.get_t()),right_potential,0,0,fronteira_livre,sig,simul);
        
        //std::cout <<left_potential*sin(frequency * simul.get_t())<<std::endl;
        double j_left = 0;
        double j_right = 0;
        
        solver.update_boundary_voltage(fronteira_livre, left_potential*sin(2*3.1416*frequency * (simul.get_t())),right_potential,0,0);

        auto start1 = std::chrono::high_resolution_clock::now();

        solver.solve();
        //std::cout << simul.get_Ez1()<<std::endl;

        auto start2 = std::chrono::high_resolution_clock::now();

        //std::cout <<a<<std::endl;
        
        double old_t = simul.get_t();
        simul.push_time(0, j_left, j_right, electron_fluxes);
        
        auto start3 = std::chrono::high_resolution_clock::now();

        if (a%10000 == 0){
            //std::cout << "voltage :" <<left_potential*sin(2*3.1416*frequency * (simul.get_t()))<<"\n";
            std::cout<< "time "<<simul.get_t()<<std::endl;
            double dt = simul.get_t() - old_t; // Get the time step for the sigma calculation
            //std::cout<< "dt "<<dt<<std::endl;

            // Saving values

            simul.write_dens(file_e_dens, 0);
            simul.write_dens(file_ar_dens, 1);
            simul.write_dens(file_arstar_dens, 2);
            simul.write_dens(file_arplus_dens, 3);

            simul.write_efield(file_e_field);
            simul.write_potential(file_pot);

            simul.write_e_energy(file_e_energy);

            file_time_steps << "Potencial on left: "<<left_potential*sin(2*3.1416*frequency * (simul.get_t())) << " dt: "<< dt << "\n";
            file_current_dens <<  sig(grid_init - 1) << " | "<< sig(grid_end - 1) << "\n";
            file_bflux << j_left << " | " << j_right<< "\n";
            //std::cout << simul.get_Ez1()<<"\n";
            //std::cout<< j_left<<"   "<< j_right<<"\n";
        }
        
        a++;

        double dt = simul.get_t() - old_t; // Get the time step for the sigma calculation

        sig(grid_init - 1) = sig(grid_init - 1) + dt * j_left; // Check com o professor se faz sentido
        sig(grid_end - 1) = sig(grid_end - 1) + dt * j_right;

        auto start4 = std::chrono::high_resolution_clock::now();
        //for (Specie& s : species){
        //    std::cout <<"especies "<<s.get_density()<<"\n";
        //}
        auto duration1 = std::chrono::duration_cast<std::chrono::nanoseconds>(start1 - start);
        auto duration2 = std::chrono::duration_cast<std::chrono::nanoseconds>(start2 - start1);
        auto duration3 = std::chrono::duration_cast<std::chrono::nanoseconds>(start3 - start2);
        auto duration4 = std::chrono::duration_cast<std::chrono::nanoseconds>(start4 - start3);
        //std::cout << "Time taken1: " << duration1.count() << " seconds" << std::endl;
        //std::cout << "Time taken2: " << duration2.count() << " seconds" << std::endl;
        //std::cout << "Time taken3: " << duration3.count() << " seconds" << std::endl;
        //std::cout << "Time taken4: " << duration4.count() << " seconds" << std::endl;
    }

     // End timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate duration in microseconds
    std::chrono::duration<double> duration = end - start;

    // Print or save the duration
    std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;



    std::cout<< "iter "<<a<<std::endl;
    simul.write_dens(file_e_dens, 0);
    
    return 0;
}
*/