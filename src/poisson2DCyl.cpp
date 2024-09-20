#define _USE_MATH_DEFINES

#include "poisson2DCyl.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <chrono>

Poisson2DCyl::Poisson2DCyl(int n, int m, double dr, double dz, std::vector<std::pair<int, double>> eps, std::map<int, double> sig) {

	r_size = n;
	z_size = m;
	r_step = dr;
	z_step = dz;

	double hdz = dz/2;
	double hdr = dr/2;

	r_grid = new double[r_size];  // Dynamically allocate an array of size n
	z_grid = new double[z_size];

	phic = new double*[n];
	phiw = new double*[n];
	phie = new double*[n];
	phin = new double*[n];
	phis = new double*[n];

	eps_vec = eps;
	sig_map = sig;

	for (int i = 0; i < n; ++i) {

        phic[i] = new double[m];
		phiw[i] = new double[m];
		phie[i] = new double[m];
		phin[i] = new double[m];
		phis[i] = new double[m];

		r_grid[i] = (i + 0.5) * dr;
    }

    for (int j = 0; j < m; j++){

    	z_grid[j] = (j+0.5)  * dz;

	}


	for (int i = 0; i < n; i++) {
		
		double dr2 = ((r_grid[i] + hdr) * (r_grid[i] + hdr) - (r_grid[i] - hdr) * (r_grid[i] - hdr));

		for (int j = 0; j < m; j++){

			phin[i][j] = 0;
			phis[i][j] = 0;
			phie[i][j] = 0;
			phiw[i][j] = 0;
			phic[i][j] = 0;

			double epsj = 0;// eps_map[(--eps_map.lower_bound(j)) -> first];
			double epsj1 = 0;//eps_map[(--eps_map.lower_bound(j+1)) -> first];

			auto it = std::upper_bound(eps_vec.begin(), eps_vec.end(), std::make_pair(j, 0.0), [](const auto& a, const auto& b) { return a.first < b.first; });
		    // Check if the iterator is valid and decrement it to get the value
		    
		    if (it != eps_vec.begin()) {
		        --it;  // Move to the largest key less than or equal to j
		        epsj = it->second;
		    }

		    auto it1 = std::upper_bound(eps_vec.begin(), eps_vec.end(), std::make_pair(j+1, 0.0), [](const auto& a, const auto& b) { return a.first < b.first; });

		    if (it1 != eps_vec.begin()) {
		        --it1;  // Move to the largest key less than or equal to j
		        epsj1 = it1->second;
		    }

			if (j < m){
				phie[i][j] =  -M_PI * dr2 * epsj * epsj1 / ((epsj1 + epsj) * hdz);
			}

			if (i < n){
				phin[i][j] = -2 * M_PI * (r_grid[i] + hdr) * dz * epsj / dr;
			}

			if (i > 0){
				phis[i][j] = phin[i-1][j];
			} else {
				phis[i][j] = phin[i][j];
			}

			if (j > 0){
				phiw[i][j] = phie[i][j-1];
			} else{
				phiw[i][j] = phie[i][j];
			}

			
			phic[i][j] = -1*(phiw[i][j] + phie[i][j] + phis[i][j] + phin[i][j]);

			//std::cout << "phie " << phie[i][j] << "phiw " << phiw[i][j]<< "phin " << phin[i][j]<< "phis " << phis[i][j]<< std::endl;

		}	
	}
	/*
    std::cout<< "grid r" <<std::endl;

    for (int i = 0; i < n; i++) {
    	std::cout<< r_grid[i]<< "; ";
    }
    
    std::cout<< std::endl<< "grid z" <<std::endl;
    for (int j = 0; j < m; j++) {
    	std::cout<< z_grid[j]<< "; ";
    }
    */
}

void Poisson2DCyl::solve(double V0, double VMAX, double VWALL, double VIN ,double (*func)(double, double), double fronteira[4], std::string str){

	auto start = std::chrono::high_resolution_clock::now();

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RHS = Eigen::MatrixXd::Zero(r_size, z_size);

	double hdr = r_step/2;
	double hdz = z_step/2;

	for (int i = 0; i < r_size; i++){

		double dr2 = ((r_grid[i] + hdr) * (r_grid[i] + hdr) - (r_grid[i] - hdr) * (r_grid[i] - hdr));

		//std::cout<<r_grid[i] + hdr<<"  "<<r_grid[i] - hdr <<"  "<<dr2<<std::endl;

		for (int j = 0; j < z_size; j++){

			RHS(i,j) = func(r_grid[i], z_grid[j]) * dr2 * M_PI * z_step;

			if (i == 0){
				if (fronteira[2]){
					phic[i][j] = phic[i][j] + phis[i][j]; //CONDICAO FRONTEIRA R = 0 NEUMANN

				}
				else{
					RHS(i,j) = RHS(i,j) - VIN * phis[i][j]; // voltage of the wall	
			
				}
			}

			else if (i == r_size-1){

				if (fronteira[3]){
					phic[i][j] = phic[i][j] + phin[i][j]; //CONDICAO FRONTEIRA R = MAX NEUMANN
				}
				else{
					RHS(i,j) = RHS(i,j) - VWALL * phin[i][j]; // voltage of the wall

				}
			}

			if (j == 0){

				if(fronteira[0]){
					phic[i][j] = phic[i][j] + phiw[i][j]; //CONDICAO FRONTEIRA z = min NEUMANN
				}
				else{
					RHS(i,j) = RHS(i,j) - V0 * phiw[i][j];
				}		 
			}

			else if (j == z_size-1){

				if (fronteira[1]){

					phic[i][j] = phic[i][j] + phie[i][j]; //CONDICAO FRONTEIRA z = MAX NEUMANN

				}else{
					RHS(i,j) = RHS(i,j) - VMAX * phie[i][j];
				}
			}
		}
	}

	auto t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> d1 = t1 - start;
    
    std::cout << "Time taken to Load RHS: " << d1.count() << " ms" << std::endl;

	if (z_size < 15 && r_size <= 7){
		std::cout << "DENSE Matrix:\n" << RHS << std::endl;
	}

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(r_size * z_size * 5);

	for (int i = 0; i < r_size; i++){

		for (int j = 0; j < z_size; j++){

			int place = i * z_size + j;

			triplets.emplace_back( place, place, phic[i][j]);

			if (  j < z_size - 1){
				triplets.emplace_back(place, place + 1, phie[i][j]);
			}

			if ( j > 0){
				triplets.emplace_back(place, place - 1, phiw[i][j]);
			}

			if (i < r_size - 1){
				triplets.emplace_back(place , place + z_size, phin[i][j]);
			}

			if (i > 0){
				triplets.emplace_back(place,  place - z_size, phis[i][j]);
			}
		}

	}

	Eigen::SparseMatrix<double> PHI(r_size * z_size, r_size * z_size);
	PHI.setFromTriplets(triplets.begin(), triplets.end());

	auto t2 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> d2 = t2 - t1;
    
    std::cout << "Time taken to Load PHI's (coeffecients): " << d2.count() << " ms" << std::endl;

	delete &triplets;

	if (z_size < 15 && r_size <= 7){

		std::cout << "Sparse Matrix:\n" << Eigen::MatrixXd(PHI) << std::endl;
	}
	Eigen::VectorXd b(RHS.size());

    b = Eigen::Map<const Eigen::VectorXd>(RHS.data(), RHS.size());

    if (z_size < 15 && r_size <= 7){

		std::cout << "Sparse Matrix column:\n" << b << std::endl;
	}

    //Eigen::SparseLU<Eigen::SparseMatrix<double>> solver; // no memory for 1000x1000
    //Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver; //VERY FAST !!
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver; // ALSO FAST (faster? 500ms for 1000x1000)

    //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver; // Very Slow

    solver.analyzePattern(PHI);

    solver.factorize(PHI);

    //solver.compute(PHI);

    auto t3 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> d3 = t3 - t2;
    
    std::cout << "Time taken to compute PHI (coeffecients): " << d3.count() << " ms" << std::endl;

    Eigen::VectorXd v = solver.solve(b);

    if (z_size < 15 && r_size <= 7){

		std::cout << "Sparse Matrix column:\n" << v << std::endl;
	}

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> V = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(v.data(), r_size, z_size);

    auto end = std::chrono::high_resolution_clock::now();
    
    // Calculate the duration in milliseconds

    std::chrono::duration<double, std::milli> d4 = end - t3;

    std::chrono::duration<double, std::milli> duration = end - start;
    
    std::cout << "Time taken to solve: " << d4.count() << " ms" << std::endl;

    std::cout << "Time taken: " << duration.count() << " ms" << std::endl;

    std::ofstream file("../" + str + ".txt");

    if (file.is_open()) {
        // Write the matrix to the file
        for (int i = 0; i < V.rows(); ++i) {
            for (int j = 0; j < V.cols(); ++j) {
                file << V(i, j);
                if (j != V.cols() - 1) {
                    file << " ";  // Space-separated values
                }
            }
            file << "\n";  // Newline after each row
        }
        file.close();
        std::cout << "Matrix written to " + str + ".txt\n";
    } else {
        std::cerr << "Unable to open file\n";
    }


    // Calculo de Campo Elétrico

    Eigen::MatrixXd E_r(r_size-1, z_size);
    Eigen::MatrixXd E1_z(r_size, z_size - 1);
    Eigen::MatrixXd E2_z(r_size, z_size - 1);

    for (int i = 0; i < r_size - 1; i++) { // Calculate the radial field (same dieletric constant)
		
		for (int j = 0; j <  z_size; j++){

			E_r(i,j) = (V(i,j) - V(i+1,j)) / (r_step);

		}		
	} 

	for (int i = 0; i < r_size; i++) {
		
		for (int j = 0; j <  z_size - 1; j++){

			double epsj1 = 0;// eps_map[(--eps_map.lower_bound(j)) -> first];
			double epsj2 = 0;//eps_map[(--eps_map.lower_bound(j+1)) -> first];

			auto it = std::upper_bound(eps_vec.begin(), eps_vec.end(), std::make_pair(j, 0.0), [](const auto& a, const auto& b) { return a.first < b.first; });
		    // Check if the iterator is valid and decrement it to get the value
		    
		    if (it != eps_vec.begin()) {
		        --it;  // Move to the largest key less than or equal to j
		        epsj1 = it->second;
		    }

		    auto it1 = std::upper_bound(eps_vec.begin(), eps_vec.end(), std::make_pair(j+1, 0.0), [](const auto& a, const auto& b) { return a.first < b.first; });

		    if (it1 != eps_vec.begin()) {
		        --it1;  // Move to the largest key less than or equal to j
		        epsj2 = it1->second;
		    }

			E1_z(i,j) = (epsj2*V(i,j) - epsj2*V(i,j+1)) / (epsj1*hdz + epsj2*hdz);
			E2_z(i,j) = (epsj1*V(i,j) - epsj1*V(i,j+1)) / (epsj1*hdz + epsj2*hdz);
			
		}
	}

	std::ofstream fileEr("../Er.txt");

    if (fileEr.is_open()) {
        // Write the matrix to the file
        for (int i = 0; i < E_r.rows(); ++i) {
            for (int j = 0; j < E_r.cols(); ++j) {
                fileEr << E_r(i, j);
                if (j != E_r.cols() - 1) {
                    fileEr << " ";  // Space-separated values
                }
            }
            fileEr << "\n";  // Newline after each row
        }
        fileEr.close();
        std::cout << "Matrix written to Er.txt\n";
    } else {
        std::cerr << "Unable to open file\n";
    }

    std::ofstream fileE1z("../E1z.txt");

    if (fileE1z.is_open()) {
        // Write the matrix to the file
        for (int i = 0; i < E1_z.rows(); ++i) {
            for (int j = 0; j < E1_z.cols(); ++j) {
                fileE1z << E1_z(i, j);
                if (j != E1_z.cols() - 1) {
                    fileE1z << " ";  // Space-separated values
                }
            }
            fileE1z << "\n";  // Newline after each row
        }
        fileE1z.close();
        std::cout << "Matrix written to E1_z.txt\n";
    } else {
        std::cerr << "Unable to open file\n";
    }

};
