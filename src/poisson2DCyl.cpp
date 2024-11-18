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


Poisson2DCyl::Poisson2DCyl(int n, int m, double dr, double dz, Eigen::VectorXd eps_initial, Eigen::VectorXd sig_initial) {

	r_size = n;
	z_size = m;
	r_step = dr;
	z_step = dz;

	t = 0;

	double hdz = dz/2;
	double hdr = dr/2;

	r_grid = new double[r_size];  // Dynamically allocate an array of size n
	z_grid = new double[z_size];

	phic = new double*[n];
	phiw = new double*[n];
	phie = new double*[n];
	phin = new double*[n];
	phis = new double*[n];

	eps = eps_initial;
	sig = sig_initial;

	//std::cout <<sig<<std::endl;
	//std::cout <<eps<<std::endl;

	S_hori = Eigen::MatrixXd::Zero(r_size, z_size);
	S_vert = Eigen::MatrixXd::Zero(r_size, z_size);
	vols = Eigen::MatrixXd::Zero(r_size, z_size);
	rho = Eigen::MatrixXd::Zero(r_size, z_size);

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

			// Calcular volumes e superficies (superiores/N e á direita/E) de celulas i,j

			S_hori(i,j) = M_PI * dr2;
			S_vert(i,j) = 2 * M_PI * (r_grid[i] + hdr) * dz;
			vols(i,j) = M_PI * dr2 * dz;

			phin[i][j] = 0;
			phis[i][j] = 0;
			phie[i][j] = 0;
			phiw[i][j] = 0;
			phic[i][j] = 0;

			double epsj = eps(j);// eps_map[(--eps_map.lower_bound(j)) -> first];
			double epsj1 = epsj;

			if (j+1 < m){
				double epsj1 = eps(j+1);//eps_map[(--eps_map.lower_bound(j+1)) -> first];
			}

			if (j < m){
				phie[i][j] =  - S_hori(i,j) * epsj * epsj1 / ((epsj1 + epsj) * hdz);
			}

			if (i < n){
				phin[i][j] = -S_vert(i,j) * epsj / dr;
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

}

void Poisson2DCyl::solve(double V0, double VMAX, double VWALL, double VIN ,Eigen::MatrixXd n_e, Eigen::MatrixXd n_i, double fronteira[4]){

	ne = n_e;
	ni = n_i;

	v0 = V0;
	vmax = VMAX;
	vwall = VWALL;
	vin = VIN;
	std::copy(fronteira, fronteira +4 , front);

	auto start = std::chrono::high_resolution_clock::now();

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RHS = Eigen::MatrixXd::Zero(r_size, z_size);

	double hdr = r_step/2;
	double hdz = z_step/2;

	for (int i = 0; i < r_size; i++){

		//double dr2 = ((r_grid[i] + hdr) * (r_grid[i] + hdr) - (r_grid[i] - hdr) * (r_grid[i] - hdr)); // ativar caso voltare a usar para calculo de rhos

		//std::cout<<r_grid[i] + hdr<<"  "<<r_grid[i] - hdr <<"  "<<dr2<<std::endl;

		for (int j = 0; j < z_size; j++){

			double epsj = eps(j);// eps_map[(--eps_map.lower_bound(j)) -> first];
			double epsj1 = epsj;

			if ( j + 1 < z_size){
				double epsj1 = eps(j + 1);//eps_map[(--eps_map.lower_bound(j+1)) -> first];
			}

			if (sig(j) != 0){

				RHS(i,j) = RHS(i,j) + (epsj * sig(j) * S_hori(i,j))/(epsj + epsj1); // Alterado o RHS pela matrize de superficie por verificar !!
		        RHS(i,j+1) = RHS(i,j+1) + (epsj1 * sig(j) * S_hori(i,j))/(epsj + epsj1); // Alterado o RHS pela matrize de superficie por verificar !!
			}

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

	rhs_i = RHS;

	auto t1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> d1 = t1 - start;
    
    std::cout << "Time taken to Load RHS: " << d1.count() << " ms" << std::endl;

	if (z_size < 15 && r_size <= 7){
		//std::cout << "DENSE Matrix:\n" << RHS << std::endl;
	}

	//std::cout << "DENSE Matrix:\n" << RHS << std::endl;

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

	phi = PHI;

	solver.analyzePattern(PHI);

    solver.factorize(PHI);

	auto t2 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> d2 = t2 - t1;
    
    std::cout << "Time taken to Load PHI's (coeffecients): " << d2.count() << " ms" << std::endl;

	//delete &triplets;

	if (z_size < 15 && r_size <= 7){

		std::cout << "Sparse Matrix:\n" << Eigen::MatrixXd(PHI) << std::endl;
	}
	
}

void Poisson2DCyl::solve_Poisson(){

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RHS = Eigen::MatrixXd::Zero(r_size, z_size);

	calculate_charge_density();

	//std::cout << vols<<std::endl;

	for (int i = 0; i < r_size; i++){
		for (int j = 0; j < z_size; j++){
			//std::cout << rho(i,j) * vols(i,j) << std::endl;
			RHS(i,j) = rhs_i(i,j) + rho(i,j) * vols(i,j);// ;// * vols(i,j); // Alterado o RHS pela matrize de volume por verificar !!
		}
	}

	//std::cout << ne << std::endl;
	//std::cout << rhs_i << std::endl;

	Eigen::VectorXd b(RHS.size());
    b = Eigen::Map<const Eigen::VectorXd>(RHS.data(), RHS.size());

    
    //std::cout << b << std::endl;

    Eigen::VectorXd v = solver.solve(b);

    //std::cout << v << std::endl;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> V = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(v.data(), r_size, z_size);    

    // Calculo de Campo Elétrico

    Eigen::MatrixXd E_r(r_size-1, z_size);
    Eigen::MatrixXd E1_z(r_size, z_size - 1);
    Eigen::MatrixXd E2_z(r_size, z_size - 1);

	double hdz = z_step/2;

	for (int i = 0; i < r_size; i++) {
		
		for (int j = 0; j <  z_size; j++){

			if (j < z_size - 1){

				double epsj1 = eps(j);
				double epsj2 = eps(j+1);

				E1_z(i,j) = (- sig(j) * hdz +epsj2*V(i,j) - epsj2*V(i,j+1)) / (epsj1*hdz + epsj2*hdz);
				E2_z(i,j) = (sig(j) * hdz +epsj1*V(i,j) - epsj1*V(i,j+1)) / (epsj1*hdz + epsj2*hdz);
			}

			if ( i < r_size -1){
				E_r(i,j) = (V(i,j) - V(i+1,j)) / (r_step);
			}

		}
	}

	Pot = V;
	Er = E_r;
	Ez1 = E1_z;
	Ez2 = E2_z;

	std::cout<<"finished Poisson"<<std::endl;
}

void Poisson2DCyl::push_time(double ti, double dt, std::ofstream& file){

	double mu = 0.001; // To be changed later !!!
	double De = 0.1;

	// double ne_N_flux; // Used to minimize calculations of fluxes | To be implemented has to be with a vector to save it and reuse them !!
	double ne_E_flux; // Used to minimize calculations of fluxes

	Eigen::MatrixXd Se = Eigen::MatrixXd::Zero(r_size, z_size); // Source term TO BE CHANGED

	//std::cout <<Se<<std::endl;

	// Calculate Fluxes at each cell
	
	Eigen::MatrixXd f_W = Eigen::MatrixXd::Zero(r_size, z_size);
	Eigen::MatrixXd f_E = Eigen::MatrixXd::Zero(r_size, z_size);
	Eigen::MatrixXd f_N = Eigen::MatrixXd::Zero(r_size, z_size);
	Eigen::MatrixXd f_S = Eigen::MatrixXd::Zero(r_size, z_size);

	Eigen::MatrixXd new_ne = Eigen::MatrixXd::Zero(r_size, z_size);
	//double f_W;
	//double f_E;
	//double f_N;
	//double f_S;

	// for testing purposes
	/*
	Er = Eigen::MatrixXd::Zero(r_size - 1, z_size);
	Ez1 = Eigen::MatrixXd::Constant(r_size, z_size - 1, 1);
	Ez2 = Eigen::MatrixXd::Constant(r_size, z_size - 1, 1);

	De = 0;
	mu = 0.625;
	*/
	double vr_max = std::max({std::fabs(mu * Ez1.maxCoeff()), std::fabs(mu * Ez2.maxCoeff())});
	double vz_max = std::max({std::fabs(mu * Er.maxCoeff()), std::fabs(mu * Er.maxCoeff())});

	dt = std::min({0.5 * z_step/vz_max, 0.5 * r_step/vr_max, 0.125* z_step*z_step/De, 0.125 * r_step*r_step, 0.5 * 8.85e-12/(mu * ne.maxCoeff()* 1.6e-19)});

	//std::cout <<Ez1<<std::endl;
	//std::cout <<ne<<std::endl;
	//std::cout<<" "<<std::endl;
	//std::cout<<Ez2<<std::endl;
	//dt = z_step;

	for (int i = 0; i < r_size; i++){

		double g_cu = 0;
		double g_dc = 0;
		double g_c = 0;
		double mf = 0;

		for (int j = 0; j < z_size; ++j){ // ter atencao n(i,j) * mu

			
			if (i == 0){

				f_S(i,j) = 0;

			} else {

				//f_S(i,j) = -(ne(i,j) * mu * Er(i - 1,j) - De * (ne(i,j) - ne(i - 1, j)) / r_step) * S_vert(i - 1,j);
				f_S(i,j) = 0;
			}

			if (i == r_size - 1){

				f_N(i,j) = 0;

			} else {

				//f_N(i,j) = (-ne(i,j) * mu * Er(i,j) - De * (ne(i + 1,j) - ne(i, j)) / r_step) * S_vert(i,j); 
				f_N(i,j) = 0;
			}

			// Handle f_E calculations using UNO3

			if (j == z_size - 1){

				f_E(i,j) = 0;

			} else {

				if (-Ez1(i, j) > 0) {

					g_dc = (ne(i, j + 1) - ne(i, j))/z_step;
				    g_cu = (ne(i, j) - ne(i, j - 1))/z_step;

				    g_c = calculate_g_c(g_dc, g_cu, -mu * Ez1(i,j), z_step, dt);
					mf = midWayFlux(ne(i,j), -mu * Ez1(i,j), z_step, dt, g_c);

				} else {

					if (j == z_size - 2) {

						g_dc = (ne(i, j) - ne(i, j + 1))/z_step;
			    		g_cu = (ne(i, j + 1) - 0)/z_step;

					} else {

						g_dc = (ne(i, j) - ne(i, j + 1))/z_step;
			    		g_cu = (ne(i, j + 1) - ne(i, j + 2))/z_step;
					
					}

					g_c = calculate_g_c(g_dc, g_cu, -mu * Ez1(i,j), z_step, dt);
					mf = midWayFlux(ne(i,j + 1), - mu * Ez1(i,j), z_step, dt, g_c);

				}

				f_E(i, j) = (-mf * mu * Ez1(i, j) - De * (ne(i, j + 1) - ne(i, j)) / z_step) * S_hori(i, j);
			}

			// Handle f_W calculations

			if (j == 0){

				f_W(i,j) = 0;

			} else {

				if (-Ez1(i, j - 1) > 0) {

					if (j == 1) {

						g_dc = (ne(i, j) - ne(i, j - 1))/z_step;
			    		g_cu = (ne(i, j - 1) - 0)/z_step;

					} else {

						g_dc = (ne(i, j) - ne(i, j - 1))/z_step;
			    		g_cu = (ne(i, j - 1) - ne(i, j - 2))/z_step;
					
					}

				    g_c = calculate_g_c(g_dc, g_cu, -mu * Ez1(i,j - 1), z_step, dt);
					mf = midWayFlux(ne(i,j - 1), -mu * Ez1(i,j - 1), z_step, dt, g_c);

				} else {

					g_dc = (ne(i, j - 1) - ne(i, j))/z_step;
			    	g_cu = (ne(i, j) - ne(i, j + 1))/z_step;

					g_c = calculate_g_c(g_dc, g_cu, -mu * Ez1(i,j - 1), z_step, dt);
					mf = midWayFlux(ne(i,j), - mu * Ez1(i,j - 1), z_step, dt, g_c);

				}

				f_W(i, j) = (-mf * mu * Ez1(i, j - 1) - De * (ne(i, j) - ne(i, j - 1)) / z_step) * S_hori(i, j - 1);
			}

			new_ne(i, j) = ne(i, j) + dt * (Se(i, j) + (f_W(i, j) - f_E(i, j) + f_S(i, j) - f_N(i, j)) / vols(i, j));
			
		}
	}

	ne = new_ne;

	//std::cout << "Fw : "<< f_W<<std::endl;
	//std::cout << "FE : "<< f_E<<std::endl;
	//std::cout<<"---"<<std::endl;
	//std::cout<<f_W<<std::endl;
	//std::cout<<"---"<<std::endl;
	//std::cout<<f_E<<std::endl;

	// 8.85e-12/(mu * ne.maxCoeff())});
	//dt = dt;
	//std::cout<< 0.5 * z_step/vz_max << " "<< 0.5 * r_step/vr_max <<" "<< 0.125* z_step*z_step/De <<" "<< 0.125 * r_step*r_step <<" "<< 0.5 * 8.85e-12/(mu * ne.maxCoeff() * 1.6e-19)<<std::endl;

	//Eigen::MatrixXd result = (f_W - f_E + f_S - f_N).array() / vols.array();

	//ne = ne + dt * (Se + result);//+(f_W - f_E + f_S - f_N)/vols(i,j));

	t = t+dt;

	std::cout<< "dt = "<<dt <<std::endl;

	solve_Poisson();

	std::cout<<"Pushed time to "<<t<<std::endl;
	//std::cout<<ne<<std::endl;

	file << "Time: " << t << "\n";
    file << rho << "\n";
    file << "----\n";  // Separator for the next matrix 
}

void Poisson2DCyl::calculate_charge_density(){
	rho = (ni - ne)* 1.6e-19;
}

void Poisson2DCyl::write_fields(std::string str){
	
	std::ofstream file("../"+str+"V.txt");

    if (file.is_open()) {
        // Write the matrix to the file
        for (int i = 0; i < Pot.rows(); ++i) {
            for (int j = 0; j < Pot.cols(); ++j) {
                file << Pot(i, j);
                if (j != Pot.cols() - 1) {
                    file << " ";  // Space-separated values
                }
            }
            file << "\n";  // Newline after each row
        }
        file.close();
        std::cout << "Matrix written to " + str + "V.txt\n";
    } else {
        std::cerr << "Unable to open file\n";
    }

	std::ofstream fileEr("../"+ str+"Er.txt");

    if (fileEr.is_open()) {
        // Write the matrix to the file
        for (int i = 0; i < Er.rows(); ++i) {
            for (int j = 0; j < Er.cols(); ++j) {
                fileEr << Er(i, j);
                if (j != Er.cols() - 1) {
                    fileEr << " ";  // Space-separated values
                }
            }
            fileEr << "\n";  // Newline after each row
        }
        fileEr.close();
        std::cout << "Matrix written to "+ str+"Er.txt\n";
    } else {
        std::cerr << "Unable to open file\n";
    }

    std::ofstream fileE1z("../" + str+"Ez1.txt");

    if (fileE1z.is_open()) {
        // Write the matrix to the file
        for (int i = 0; i < Ez1.rows(); ++i) {
            for (int j = 0; j < Ez1.cols(); ++j) {
                fileE1z << Ez1(i, j);
                if (j != Ez1.cols() - 1) {
                    fileE1z << " ";  // Space-separated values
                }
            }
            fileE1z << "\n";  // Newline after each row
        }
        fileE1z.close();
        std::cout << "Matrix written to "+ str+"Ez1.txt\n";
    } else {
        std::cerr << "Unable to open file\n";
    }
}


// private functions

int Poisson2DCyl::getSign(double a, double b) {

    if (a > b) {
        return 1;  // Positive
    } else if (a < b) {
        return -1; // Negative
    } else {
        return 0;  // Zero
    }
}

double Poisson2DCyl::midWayFlux(double donor_cell, double u, double dx, double dt, double g_c){
	return donor_cell + 0.5 * getSign(u,0)*(dx - std::abs(u) * dt) * g_c;
}

double Poisson2DCyl::calculate_g_c(double g_dc, double g_cu, double u, double dx, double dt){

	if (std::abs(g_dc - g_cu) <= 0.6*(g_dc + g_cu)){
		
		return g_dc - 4/3 * (dx - 0.5 * getSign(u,0)*(dx - std::abs(u)*dt))*(g_dc - g_cu)/(2 *dx);
	
	} else if (g_dc * g_cu > 0 ){
		
		return getSign(g_dc,0)*2*std::min(g_dc,g_cu);
	
	}else{

		return getSign(g_dc,0)*std::min(g_dc,g_cu);
	}
}