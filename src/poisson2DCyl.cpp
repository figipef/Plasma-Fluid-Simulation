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
#include <iomanip> // For std::setprecision and std::fixed
 
Poisson2DCyl::Poisson2DCyl(int n, int m, double dr, double dz, Eigen::VectorXd eps_initial, Eigen::VectorXd sig_initial) {

	r_size = n;
	z_size = m;
	r_step = dr;
	z_step = dz;

	t = 0;

	double hdz = dz*0.5;
	double hdr = dr*0.5;

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

			//S_hori(i,j) = dr*dh;
			//S_vert(i,j) = dz*dh;
			//vols(i,j) = dz*dh*dr;

			phin[i][j] = 0;
			phis[i][j] = 0;
			phie[i][j] = 0;
			phiw[i][j] = 0;
			phic[i][j] = 0;

			double epsj = eps(j);// eps_map[(--eps_map.lower_bound(j)) -> first];
			double epsj1 = epsj;

			if (j+1 < m){
				epsj1 = eps(j+1);//eps_map[(--eps_map.lower_bound(j+1)) -> first];
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

			if (r_size == 1){
				phic[i][j] = -1*(phiw[i][j] + phie[i][j]);
			} else {
				phic[i][j] = -1*(phiw[i][j] + phie[i][j] + phis[i][j] + phin[i][j]);
			}
		}	
	}
}

void Poisson2DCyl::set_pol_coeffs(Eigen::VectorXd mob_coeffs, Eigen::VectorXd temp_coeffs, int degree){
	pol_degree = degree;
	mob_pol_vec = mob_coeffs; 
	temp_pol_vec = temp_coeffs;
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

	for (int i = 0; i < r_size; i++){

		for (int j = 0; j < z_size; j++){

			double epsj = eps(j);// eps_map[(--eps_map.lower_bound(j)) -> first];
			double epsj1 = epsj;

			if ( j + 1 < z_size){
				epsj1 = eps(j + 1);//eps_map[(--eps_map.lower_bound(j+1)) -> first];
			}

			if (sig(j) != 0){

				RHS(i,j) = RHS(i,j) + (epsj * sig(j) * S_hori(i,j))/(epsj + epsj1); // Alterado o RHS pela matrize de superficie por verificar !!
		        RHS(i,j+1) = RHS(i,j+1) + (epsj1 * sig(j) * S_hori(i,j))/(epsj + epsj1); // Alterado o RHS pela matrize de superficie por verificar !!
			}

			if (i == 0){
				if (fronteira[2]){
					if (r_size > 1){
						phic[i][j] = phic[i][j] + phis[i][j]; //CONDICAO FRONTEIRA R = 0 NEUMANN
					}
					
				}else{
					RHS(i,j) = RHS(i,j) - VIN * phis[i][j]; // voltage of the wall	
			
				}
			}

			else if (i == r_size-1){

				if (fronteira[3]){
					if (r_size > 1){
						phic[i][j] = phic[i][j] + phin[i][j]; //CONDICAO FRONTEIRA R = MAX NEUMANN
					}
					
				}else{
					RHS(i,j) = RHS(i,j) - VWALL * phin[i][j]; // voltage of the wall

				}
			}

			if (j == 0){

				if(fronteira[0]){
					phic[i][j] = phic[i][j] + phiw[i][j]; //CONDICAO FRONTEIRA z = min NEUMANN
				}else{
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
}

void Poisson2DCyl::solve_Poisson(){

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RHS = Eigen::MatrixXd::Zero(r_size, z_size);
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> HELP = Eigen::MatrixXd::Zero(r_size, z_size);

	calculate_charge_density();

	for (int i = 0; i < r_size; i++){
		for (int j = 0; j < z_size; j++){

			RHS(i,j) = rhs_i(i,j) + rho(i,j) * vols(i,j);
		}
	}

	Eigen::VectorXd b(RHS.size());
    b = Eigen::Map<const Eigen::VectorXd>(RHS.data(), RHS.size());

    Eigen::VectorXd v = solver.solve(b);

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

				E1_z(i,j) = (- sig(j) * hdz + epsj2*V(i,j) - epsj2*V(i,j+1)) / (epsj1*hdz + epsj2*hdz); 
				E2_z(i,j) = (sig(j) * hdz +epsj1*V(i,j) - epsj1*V(i,j+1)) / (epsj1*hdz + epsj2*hdz);

				//E1_z(i,j) = (V(i,j) - V(i,j+1)) / (z_step); 
				//E2_z(i,j) = (V(i,j) - V(i,j+1)) / (z_step);
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
}

void Poisson2DCyl::push_time(int ti, double dt,std::ofstream& dt_file, int  int_mode){

	//double mu = 0.03; // To be changed later !!!
	//double De = 0.1;

	Eigen::MatrixXd mu = Eigen::MatrixXd::Zero(r_size, z_size);
	Eigen::MatrixXd De = Eigen::MatrixXd::Zero(r_size, z_size);
	Eigen::MatrixXd Se = Eigen::MatrixXd::Zero(r_size, z_size);

	//std::cout <<Ez1<<std::endl;

	for (int i = 0; i < r_size; i++){
	
		for (int j = 0; j < z_size; ++j){ 

			double x;
			
			if (j == 0){
				x = std::abs(Ez1(i,j) / (2.43e25) * 1e21);
			}else if (j == z_size - 1){
				x = std::abs(Ez1(i,j-1) / (2.43e25) * 1e21);
			}else{
				x = std::abs((Ez1(i,j) + Ez1(i,j-1))/2 / (2.43e25) * 1e21);
			}

			

			mu(i,j) = Pol(x, 0)/2.43e25;

			De(i,j) = Pol(x, 1)/2.43e25; 

			Se(i,j) = Pol(x, 2) * 2.43e25 * ne(i,j); 
		}
	}

	// Calculate Fluxes at each cell
	
	Eigen::MatrixXd new_ne = Eigen::MatrixXd::Zero(r_size, z_size);
	Eigen::MatrixXd new_ni = Eigen::MatrixXd::Zero(r_size, z_size);

	Eigen::MatrixXd new_new_ne = Eigen::MatrixXd::Zero(r_size, z_size);
	Eigen::MatrixXd new_new_ni = Eigen::MatrixXd::Zero(r_size, z_size);

	Eigen::MatrixXd midfluxne = Eigen::MatrixXd::Zero(r_size, z_size);
	Eigen::MatrixXd midfluxni = Eigen::MatrixXd::Zero(r_size, z_size);

	double vr_max = 1e-308;
	if ( r_size > 1){
		vr_max = std::abs(mu.maxCoeff() * Er.maxCoeff());
	}

	double vz_max = std::max({std::abs(mu.maxCoeff() * Ez1.maxCoeff()), std::abs(mu.maxCoeff() * Ez2.maxCoeff()), 1e-308});

	//dt = std::min({0.125 * z_step/vz_max, 0.25* z_step*z_step/De, 0.5 * 8.85e-12/(mu * ne.maxCoeff()* 1.6e-19)});
	
	double A_cvt = 0.1;
	double A_dif = 0.1;

	if (r_size == 1){
	
		dt = std::min({A_cvt * z_step/(vz_max + 1e-308), A_dif* z_step*z_step/(De.maxCoeff() + 1e-308)});
	
	} else {
	
		dt = std::min({A_cvt * z_step/vz_max, A_cvt * r_step/vr_max, A_dif* z_step*z_step/De.maxCoeff(), A_dif * r_step*r_step/De.maxCoeff()});
	
	}
	//dt = 1e-13;
	dt = 20e-12;  

	if (int_mode == 0){

		for (int i = 0; i < r_size; i++){
	
			for (int j = 0; j < z_size; ++j){ // changed from 0 and z _size to 1 and z_size - 1
	
				midfluxne(i,j) =  calcFlux_UNO3(i, j, mu(i,j), De(i,j), dt, ne);
	
				//midfluxni(i,j) = calcFlux(i, j, 0, 0, dt, ni, alph);
	
				new_ne(i, j) = ne(i, j) + dt * (midfluxne(i, j)/vols(i, j) + Se(i,j));
	
				new_ni(i, j) = ni(i, j) + dt * Se(i,j);//midfluxni(i, j);

			}
		}

	} else if (int_mode == 1){

		for (int i = 0; i < r_size; i++){
	
			for (int j = 0; j < z_size; ++j){ 
	
				midfluxne(i,j) = calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, ne);
	
				//midfluxni(i,j) = calcFlux_superbee(i, j, 0, 0, dt, ni, alph);
	
				new_new_ne(i, j) = ne(i, j) + dt * midfluxne(i, j)/vols(i, j);
	
				//new_new_ni(i, j) = ni(i, j) + dt * midfluxni(i, j);

			}
		}

		for (int i = 0; i < r_size; i++){
	
			for (int j = 0; j < z_size; ++j){ // ter atencao n(i,j) * mu
	
				new_ne(i,j) = ne(i, j) + dt/2*(midfluxne(i, j) + calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, new_new_ne))/vols(i, j);

				//new_ni(i,j) = ni(i, j) + dt/2*(midfluxni(i, j) + calcFlux_superbee(i, j,  0,  0, dt, new_new_ni, alph))/vols(i, j);

				new_ni(i, j) = ni(i, j) + dt * Se(i,j);//midfluxni(i, j);
			}	
		}

	} else if (int_mode == 2){

		Eigen::MatrixXd k1 = Eigen::MatrixXd::Zero(r_size, z_size);
		Eigen::MatrixXd k2 = Eigen::MatrixXd::Zero(r_size, z_size);
		Eigen::MatrixXd k3 = Eigen::MatrixXd::Zero(r_size, z_size);
		Eigen::MatrixXd k4 = Eigen::MatrixXd::Zero(r_size, z_size);

		Eigen::MatrixXd helpk1 = Eigen::MatrixXd::Zero(r_size, z_size);
		Eigen::MatrixXd helpk2 = Eigen::MatrixXd::Zero(r_size, z_size);
		Eigen::MatrixXd helpk3 = Eigen::MatrixXd::Zero(r_size, z_size);
		Eigen::MatrixXd helpk4 = Eigen::MatrixXd::Zero(r_size, z_size);

		for (int i = 0; i < r_size; i++){

			for (int j = 0; j < z_size; ++j){ 
	
				k1(i,j) = calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, ne);
	
			}
		}
	
		helpk1 = ne + k1 * dt/2.;
	
		for (int i = 0; i < r_size; i++){
	
			for (int j = 0; j < z_size; ++j){ 
	
				k2(i,j) = calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, helpk1);
	
			}
		}
	
		helpk2 = ne + k2 * dt/2.;
	
		for (int i = 0; i < r_size; i++){
	
			for (int j = 0; j < z_size; ++j){ 
	
				k3(i,j) = calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, helpk2);
	
			}
		}
	
		helpk3 = ne + k3 * dt; // Atencao, provavelmente adicionar /vols() ????
	
		for (int i = 0; i < r_size; i++){
	
			for (int j = 0; j < z_size; ++j){ 
	
				k4(i,j) = calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, helpk3);
			}
		}
		
	
		for (int i = 0; i < r_size; i++){
	
			for (int j = 0; j < z_size; ++j){ 
	
				new_ne(i,j) = ne(i,j) + dt*(((k1(i,j) + 2*k2(i,j) + 2*k3(i,j) + k4(i,j))/vols(i, j))/6 + Se(i,j));

				new_ni(i, j) = ni(i, j) + dt * Se(i,j);
			}
		}
	
	}

	ne = new_ne;
	ni = new_ni;

	t = t+dt;

	solve_Poisson();

	write_dt(dt_file, dt);

	if ( ti % 999 == 0){
		//std::cout <<De<<std::endl;
		//std::cout <<Se*dt<<std::endl;

		std::cout << " t = "<<t<<std::endl;
		std::cout << " dt =" <<dt<<std::endl;
	}	
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

void Poisson2DCyl::write_dens(std::ofstream& file){
	//std::cout << std::fixed << std::setprecision(21);
	//std::cout<<"Density writen at "<<Pot<<std::endl;
	file << "Time: " << t << "\n";
    file << ne << "\n";
    file << "----\n";  // Separator for the next matrix 
}

void Poisson2DCyl::write_dt(std::ofstream& file, double dt){
    file << dt << " ";
}

//===============================

//      Private Functions

//===============================

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
	return donor_cell + 0.5 * std::copysign(1.0,u)*(dx - std::abs(u) * dt) * g_c;
}

double Poisson2DCyl::calculate_g_c(double g_dc, double g_cu, double u, double dx, double dt){
	
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

double Poisson2DCyl::calculate_g_c_UNO2(double g_dc, double g_cu){

	return std::copysign(1.0, g_dc)*2*(std::abs(g_dc * g_cu)/(std::abs(g_dc) + std::abs(g_cu)+ 1e-308));
}

double Poisson2DCyl::phia(double x){
	//return std::max(0., std::min(std::min(1., (2+x)/6 ), x));
	return std::max({0., std::max(std::min(1., x*2), std::min(x, 2.))});
}

double Poisson2DCyl::korenLimiter(double x){
	return std::max(0., std::min({1.0, (2.0+x)/6.0 , x}));
}

double Poisson2DCyl::calcFlux_superbee(int i , int j ,double mu, double De, double dt, Eigen::MatrixXd& n){
	int currlim = 1;
	double EPS0 = 8.85418781762e-12;
	double ECHARGE = 1.6e-19;

	double flux_left = 0.0;
	double flux_right = 0.0;

	if (j == 0) {

		double re = (n(i,j) - 0) / (n(i,j+1) - n(i,j) + 1e-8); 

        double ce = (-dt * mu * Ez1(i,j)) / z_step;

        double re_prime = (n(i,j+1) - n(i,j+2)) / (n(i,j) - n(i,j+1) + 1e-8);
        double ce_prime = -ce;

        double nem;
        if (mu * Ez1(i,j) <= 0) {
            nem = n(i,j) + (phia(re) / 2) * (1 - ce) * (n(i,j+1) - n(i,j));
        } else {
            nem = n(i,j+1) + (phia(re_prime) / 2) * (1 - ce_prime) * (n(i,j) - n(i,j+1));
        }

        flux_left = 0.0;
        flux_right = -mu * Ez1(i,j) * nem - (De / z_step) * (n(i,j+1) - n(i,j));

        if (currlim == 1){
        	double EField_star_right = std::max(std::abs(Ez1(i,j)), (De * std::abs(n(i,j+1) - n(i,j))) / (mu*z_step * std::max({n(i,j), n(i,j+1), 1e-8})));

        	if (std::abs(flux_right) > (EPS0 * EField_star_right) / (ECHARGE * dt)) {
            	flux_right = (std::copysign(1.0, flux_right) * EPS0 * EField_star_right) / (ECHARGE * dt);
        	}
        }
        

    } else if (j == z_size - 1) {

        double rw = (n(i,j-1) - n(i,j-2)) / (n(i,j) - n(i,j-1) + 1e-8);
        double cw = (-dt * mu * Ez1(i,j-1)) / z_step;

        double rw_prime;
        if ( j == z_size - 2){
        	rw_prime = (n(i,j) - 0) / (n(i,j-1) - n(i,j) + 1e-8);
        } else {
        	rw_prime = (n(i,j) - 0) / (n(i,j-1) - n(i,j) + 1e-8);
        }

        double cw_prime = -cw;

        double nw;

        if (mu * Ez1(i,j-1) <= 0) {
            nw = n(i,j-1) + (phia(rw) / 2) * (1 - cw) * (n(i,j) - n(i,j-1));
        } else {
            nw = n(i,j) + (phia(rw_prime) / 2) * (1 - cw_prime) * (n(i,j-1) - n(i,j));
        }

        flux_left = -mu * Ez1(i,j-1) * nw - (De / (z_step) * (n(i,j) - n(i,j-1)));
        flux_right = 0;
        if (currlim == 1){
        	double EField_star_left = std::max( std::abs(Ez1(i,j-1)), (De * std::abs(n(i,j) - n(i,j-1))) / (mu*z_step * std::max({n(i,j-1), n(i,j), 1e-8})));

	        if (std::abs(flux_left) > (EPS0 * EField_star_left) / (ECHARGE * dt)) {
	            flux_left = (std::copysign(1.0, flux_left) * EPS0 * EField_star_left) / (ECHARGE * dt);
	        }
        }
        

    } else {
        // Initialize variables
    	double re, ce, rw, cw, re_prime, ce_prime, rw_prime, cw_prime;
    	double ne, nw;
    	double EField_star_right, EField_star_left;
	
    	// Calculate re, ce, rw, cw, and their derivatives
    	re = (n(i,j) - n(i,j-1)) / (n(i,j+1) - n(i,j) + 1e-8);
    	ce = (-dt * mu * Ez1(i,j)) / z_step;

		if (j == 1){
			rw = (n(i,j-1) - 0) / (n(i,j) - n(i,j-1) + 1e-8);
		} else {
			rw = (n(i,j-1) - n(i,j-2)) / (n(i,j) - n(i,j-1) + 1e-8);
		}
    	
    	cw = (-dt * mu * Ez1(i,j-1)) / z_step;

		if ( j == z_size - 2){
			re_prime = (n(i,j+1) - 0) / (n(i,j) - n(i,j+1) + 1e-8);
		} else {
			re_prime = (n(i,j+1) - n(i,j+2)) / (n(i,j) - n(i,j+1) + 1e-8);
		}
    	
    	ce_prime = -ce;
	
    	rw_prime = (n(i,j) - n(i,j+1)) / (n(i,j-1) - n(i,j) + 1e-8);
    	cw_prime = -cw;
	
    	// Compute ne
    	if (mu * Ez1(i,j) <= 0) { // Zero or positive velocity
    	    ne = n(i,j) + (phia(re) / 2.0) * (1.0 - ce) * (n(i,j+1) - n(i,j));
    	} else {
    	    ne = n(i,j+1) + (phia(re_prime) / 2.0) * (1.0 - ce_prime) * (n(i,j) - n(i,j+1));
    	}
	
    	// Compute nw
    	if (mu * Ez1(i,j-1) <= 0) {
    	    nw = n(i,j-1) + (phia(rw) / 2.0) * (1.0 - cw) * (n(i,j) - n(i,j-1));
    	} else {
    	    nw = n(i,j) + (phia(rw_prime) / 2.0) * (1.0 - cw_prime) * (n(i,j-1) - n(i,j));
    	}
	
    	// Compute fluxes
    	flux_left = -mu * Ez1(i,j-1) * nw - (De / z_step) * (n(i,j) - n(i,j-1));
	
    	flux_right = -mu * Ez1(i,j) * ne - (De / z_step) * (n(i,j+1) - n(i,j));
		
		if ( currlim == 1){
			// Compute EField_star_right and adjust flux_right if needed
    		EField_star_right = std::max(std::abs(Ez1(i,j)), (De * std::abs(n(i,j+1) - n(i,j))) / (mu*z_step * std::max({n(i,j),n(i,j+1), 1e-8})));
	
			//std::cout << flux_left<< " "<< flux_right <<std::endl;
	
    		if (std::abs(flux_right) > (EPS0 * EField_star_right) / (ECHARGE * dt)) {
    		    flux_right = std::copysign((EPS0 * EField_star_right) / (ECHARGE * dt), flux_right);
    		}
		
    		// Compute EField_star_left and adjust flux_left if needed
    		EField_star_left = std::max(std::abs(Ez1(i,j-1)),(De * std::abs(n(i,j) - n(i,j-1))) / (mu*z_step * std::max({n(i,j-1), n(i,j), 1e-8})));
		
    		if (std::abs(flux_left) > (EPS0 * EField_star_left) / (ECHARGE * dt)) {
    		    flux_left = std::copysign((EPS0 * EField_star_left) / (ECHARGE * dt), flux_left);
    		}
		}
    }
    //std::cout << flux_left<< " "<< flux_right <<std::endl;
    return flux_left - flux_right;    
}

double Poisson2DCyl::calcFlux_UNO3(int i , int j ,double mu, double De, double dt, Eigen::MatrixXd& n){
	int currlim = 1;
	double EPS0 = 8.85418781762e-12;
	double ECHARGE = 1.6e-19;

	double flux_left = 0.0;
	double flux_right = 0.0;
    
	if (j == 0) { // Só precisamos de ne
		double g_c;
		double g_dc;
		double g_cu;
		double nem;
		double ve = -mu * Ez1(i,j);

        if (mu * Ez1(i,j) <= 0) { // Velocidade nula ou positiva

        	g_dc = (n(i, j + 1) - n(i, j))/z_step;
			g_cu = (n(i, j) - 0)/z_step;

			g_c = calculate_g_c(g_dc, g_cu, ve, z_step, dt);
		    nem = n(i,j) + 0.5 * std::copysign(1.0, ve) * (z_step - std::abs(ve) * dt) * g_c;

        } else {

        	g_dc = -1.0 * (n(i, j) - n(i, j + 1))/z_step;
			g_cu = -1.0 *(n(i, j + 1) - n(i, j + 2))/z_step;

			g_c = calculate_g_c(g_dc, g_cu, ve, z_step, dt);
		    nem = n(i,j+1) +  0.5 * std::copysign(1.0, ve) * (z_step - std::abs(ve) * dt) * g_c;
        }

        flux_left = 0.0;
        flux_left = -mu * Ez1(i,j) * nem - (De / z_step) * (n(i,j+1) - n(i,j)); // CHANGED
        flux_right = -mu * Ez1(i,j) * nem - (De / z_step) * (n(i,j+1) - n(i,j));

        if ( currlim == 1){
        	double EField_star_right = std::max(std::abs(Ez1(i,j)), (De * std::abs(n(i,j+1) - n(i,j))) / (mu*z_step * std::max({n(i,j), n(i,j+1), 1e-8})));

        	if (std::abs(flux_right) > (EPS0 * EField_star_right) / (ECHARGE * dt)) {
        	    flux_right = (std::copysign(1.0, flux_right) * EPS0 * EField_star_right) / (ECHARGE * dt);
        	    flux_left = flux_right; // CHANGED
        	}
        }
        
    } else if (j == z_size - 1) { // Só precisamos de nw
    	double g_c;
		double g_dc;
		double g_cu;
		double nw;
		double vw = -mu * Ez1(i,j-1);
        if (mu * Ez1(i,j-1) <= 0) { // Velocidade nula ou positiva

			g_dc = (n(i, j) - n(i, j - 1))/z_step;
			g_cu = (n(i, j - 1) - n(i, j - 2))/z_step;

			g_c = calculate_g_c(g_dc, g_cu, vw, z_step, dt);
		    nw = n(i,j-1) + 0.5 * std::copysign(1.0, vw) * (z_step - std::abs(vw) * dt) * g_c;

        } else {

			g_cu = -1.0 * (n(i, j) - 0)/z_step;
			g_dc = -1.0 * (n(i, j - 1) - n(i, j))/z_step;

		    //G_A = calculate_g_c(G_WP, G_PE, vw, z_step, dt);
		    g_c = calculate_g_c(g_dc, g_cu, vw, z_step, dt);
		    nw = n(i,j) + 0.5 * std::copysign(1.0, vw) * (z_step - std::abs(vw) * dt) * g_c;
        }

        flux_left = -mu * Ez1(i,j-1) * nw - (De / z_step) * (n(i,j) - n(i,j-1));
    	flux_right = 0;
    	flux_right = flux_left; // CHANGED
    	if ( currlim == 1){
    		double EField_star_left = std::max( std::abs(Ez1(i,j-1)), (De * std::abs(n(i,j) - n(i,j-1))) / (mu*z_step * std::max({n(i,j-1), n(i,j), 1e-8})));

        	if (std::abs(flux_left) > (EPS0 * EField_star_left) / (ECHARGE * dt)) {
        	    flux_left = (std::copysign(1.0, flux_left) * EPS0 * EField_star_left) / (ECHARGE * dt);
        	    flux_right = flux_left; // CHANGED
        	}
    	}
	
    } else { // Precisamos de ne e de nw

		double ve = -mu * Ez1(i,j);
		double vw = -mu * Ez1(i,j-1);

		double g_c = 0;
		double g_cu = 0;
		double g_dc = 0;
		double nem = 0.0;
		double nw = 0.0;
		
		// Handle velocity for ne
		if (mu * Ez1(i,j) <= 0) { // Zero or positive velocity

			g_dc = (n(i, j + 1) - n(i, j))/z_step;
			g_cu = (n(i, j) - n(i, j - 1))/z_step;

			g_c = calculate_g_c(g_dc, g_cu, ve, z_step, dt);
		    nem = n(i,j) + 0.5 * std::copysign(1.0, ve) * (z_step - std::abs(ve) * dt) * g_c;



		} else {

			g_dc = -1.0 * (n(i, j) - n(i, j + 1))/z_step;

			if ( j == z_size - 2){
				g_cu = -1.0 *(n(i, j + 1) -0)/z_step;
			}else{
				g_cu = -1.0 *(n(i, j + 1) - n(i, j + 2))/z_step;
			}
	    	
		    g_c = calculate_g_c(g_dc, g_cu, ve, z_step, dt);

		    nem = n(i,j+1) +  0.5 * std::copysign(1.0, ve) * (z_step - std::abs(ve) * dt) * g_c;
		   	
		}
		
		// Handle velocity for nw
		if (mu * Ez1(i,j-1) <= 0) { // Zero or positive velocity
			
			g_dc = (n(i, j) - n(i, j - 1))/z_step;

			if ( j == 1){
				g_cu = (n(i, j - 1) - 0)/z_step;
			}else{
				g_cu = (n(i, j - 1) - n(i, j - 2))/z_step;
			}

		    g_c = calculate_g_c(g_dc, g_cu, vw, z_step, dt);
		    nw = n(i,j-1) + 0.5 * std::copysign(1.0, vw) * (z_step - std::abs(vw) * dt) * g_c;



		} else {
			g_cu = -1.0 * (n(i, j) - n(i, j + 1))/z_step;
			g_dc = -1.0 * (n(i, j - 1) - n(i, j))/z_step;

		    g_c = calculate_g_c(g_dc, g_cu, vw, z_step, dt);
		    nw = n(i,j) + 0.5 * std::copysign(1.0, vw) * (z_step - std::abs(vw) * dt) * g_c;
		}

		// Compute fluxes
    	flux_left = vw * nw - (De / z_step) * (n(i,j) - n(i,j-1));
	
    	flux_right = ve * nem - (De / z_step) * (n(i,j+1) - n(i,j));
		
		if ( currlim == 1){
			// Compute EField_star_right and adjust flux_right if needed
    		double EField_star_right = std::max(std::abs(Ez1(i,j)), (De * std::abs(n(i,j+1) - n(i,j))) / (mu*z_step * std::max({n(i,j),n(i,j+1), 1e-8})));
	
			//std::cout << flux_left<< " "<< flux_right <<std::endl;
	
    		if (std::abs(flux_right) > (EPS0 * EField_star_right) / (ECHARGE * dt)) {
    		    flux_right = std::copysign((EPS0 * EField_star_right) / (ECHARGE * dt), flux_right);
    		}
		
    		// Compute EField_star_left and adjust flux_left if needed
    		double EField_star_left = std::max(std::abs(Ez1(i,j-1)),(De * std::abs(n(i,j) - n(i,j-1))) / (mu*z_step * std::max({n(i,j-1), n(i,j), 1e-8})));
		
    		if (std::abs(flux_left) > (EPS0 * EField_star_left) / (ECHARGE * dt)) {
    		    flux_left = std::copysign((EPS0 * EField_star_left) / (ECHARGE * dt), flux_left);
    		}
		}
    }

    return (flux_left - flux_right)*S_hori(i,j);
}

double Poisson2DCyl::calcFlux_Koren(int i , int j ,double mu, double De, double dt, Eigen::MatrixXd& n){
	int currlim = 1;
	double EPS0 = 8.85418781762e-12;
	double ECHARGE = 1.6e-19;

	double flux_left = 0.0;
	double flux_right = 0.0;

	if (j == 0) { // Só precisamos de ne

		double ve = -mu * Ez1(i,j);

        if (ve >= 0) { // Velocidade nula ou positiva

        	double r_i = (n(i,j) - 0.0 + 1e-308)/(n(i,j+1) - n(i,j) + 1e-308);
        	flux_right = ve * (n(i,j) + korenLimiter(r_i)*(n(i,j+1) - n(i,j)));

        } else {

        	double r_i = (n(i,j+1) - n(i, j) + 1e-308)/(n(i,j+2) - n(i,j+1) + 1e-308);
        	flux_right = ve * (n(i,j+1) - korenLimiter(1/r_i)*(n(i,j+1) - n(i,j)));

        }

        flux_left = 0.0;
        flux_right = flux_right - (De / z_step) * (n(i,j+1) - n(i,j));

        if ( currlim == 1){
        	double EField_star_right = std::max(std::abs(Ez1(i,j)), (De * std::abs(n(i,j+1) - n(i,j))) / (mu*z_step * std::max({n(i,j), n(i,j+1), 1e-8})));

        	if (std::abs(flux_right) > (EPS0 * EField_star_right) / (ECHARGE * dt)) {
        	    flux_right = (std::copysign(1.0, flux_right) * EPS0 * EField_star_right) / (ECHARGE * dt);
        	}
        }
        
    } else if (j == z_size - 1) { // Só precisamos de nw

		double vw = -mu * Ez1(i,j-1);

        if (vw >= 0) { // Velocidade nula ou positiva

        	double r_i = (n(i,j - 1) - n(i, j - 2) + 1e-308)/(n(i,j) - n(i,j - 1) + 1e-308);
        	flux_left = vw * (n(i,j - 1) + korenLimiter(r_i)*(n(i,j) - n(i,j - 1)));

        } else {

        	double r_i = (n(i,j) - n(i, j - 1) + 1e-308)/(0.0 - n(i,j) + 1e-308);
        	flux_left = vw * (n(i,j) - korenLimiter(1/r_i)*(n(i,j) - n(i,j - 1)));

        }

        flux_left = flux_left - (De / z_step) * (n(i,j) - n(i,j-1));
    	flux_right = 0;

    	if ( currlim == 1){
    		double EField_star_left = std::max( std::abs(Ez1(i,j-1)), (De * std::abs(n(i,j) - n(i,j-1))) / (mu*z_step * std::max({n(i,j-1), n(i,j), 1e-8})));

        	if (std::abs(flux_left) > (EPS0 * EField_star_left) / (ECHARGE * dt)) {
        	    flux_left = (std::copysign(1.0, flux_left) * EPS0 * EField_star_left) / (ECHARGE * dt);
        	}
    	}
	
    } else { // Precisamos de ne e de nw

		double ve = -mu * Ez1(i,j);
		double vw = -mu * Ez1(i,j-1);
		
		// Handle velocity for ne
		if (ve >= 0) { // Zero or positive velocity

			double r_i = (n(i,j) - n(i, j -1) + 1e-308)/(n(i,j+1) - n(i,j) + 1e-308);
        	flux_right = ve * (n(i,j) + korenLimiter(r_i)*(n(i,j+1) - n(i,j)));

		} else {

			double r_i;
        	
			if ( j == z_size - 2){
				r_i = (n(i,j+1) - n(i, j)+ 1e-308)/(0.0 - n(i,j+1) + 1e-308);
			}else{
				r_i = (n(i,j+1) - n(i, j)+ 1e-308)/(n(i,j+2) - n(i,j+1) + 1e-308);
			}

			flux_right = ve * (n(i,j+1) - korenLimiter(1/r_i)*(n(i,j+1) - n(i,j)));
		}
		
		// Handle velocity for nw
		if (vw >= 0) { // Zero or positive velocity
			
			double r_i;

			if ( j == 1){
				r_i = (n(i,j - 1) - 0.0 + 1e-308)/(n(i,j) - n(i,j - 1) + 1e-308);
			}else{
				r_i = (n(i,j - 1) - n(i, j - 2) + 1e-308)/(n(i,j) - n(i,j - 1) + 1e-308);
			}

			
        	flux_left = vw * (n(i,j - 1) + korenLimiter(r_i)*(n(i,j) - n(i,j - 1)));

		} else {

			double r_i = (n(i,j) - n(i, j - 1) + 1e-308)/(n(i,j + 1) - n(i,j) + 1e-308);
        	flux_left = vw * (n(i,j) - korenLimiter(1/r_i)*(n(i,j) - n(i,j - 1)));

		}

		// Compute fluxes
    	flux_left = flux_left - (De / z_step) * (n(i,j) - n(i,j-1));
	
    	flux_right = flux_right - (De / z_step) * (n(i,j+1) - n(i,j));
		
		if ( currlim == 1){
			// Compute EField_star_right and adjust flux_right if needed
    		double EField_star_right = std::max(std::abs(Ez1(i,j)), (De * std::abs(n(i,j+1) - n(i,j))) / (mu*z_step * std::max({n(i,j),n(i,j+1), 1e-8})));
	
			//std::cout << flux_left<< " "<< flux_right <<std::endl;
	
    		if (std::abs(flux_right) > (EPS0 * EField_star_right) / (ECHARGE * dt)) {
    		    flux_right = std::copysign((EPS0 * EField_star_right) / (ECHARGE * dt), flux_right);
    		}
		
    		// Compute EField_star_left and adjust flux_left if needed
    		double EField_star_left = std::max(std::abs(Ez1(i,j-1)),(De * std::abs(n(i,j) - n(i,j-1))) / (mu*z_step * std::max({n(i,j-1), n(i,j), 1e-8})));
		
    		if (std::abs(flux_left) > (EPS0 * EField_star_left) / (ECHARGE * dt)) {
    		    flux_left = std::copysign((EPS0 * EField_star_left) / (ECHARGE * dt), flux_left);
    		}
		}
    }

    return (flux_left - flux_right)*S_hori(i,j);
}

double Poisson2DCyl::Pol(double x, int type){
	if (type == 0){
		// Mobility coef

		if (x < 1e-3){
    	    return 0.3743E+26;
		}else{
    	    return -1.3181433835569645e+25 * atan(log10(x) * 2.45634702952387 + 2.6338705981932145) + 2.0018811537266833e+25;
		}
//-1.3181433835569645e+25 2.45634702952387 -2.6338705981932145 2.0018811537266833e+25
	} else if (type == 1) {

		// Diff coef

    	if (x <= 0.1){

        	return 9.798e23;

    	}else if (x > 0.1 && x <= 10*sqrt(10)){

        	return 2.9e23 *atan(2.818*log10(x) - 2.2653e-1) + 1.937e24;

    	} else if ( x > 10*sqrt(10) && x <= 100){

	    	return 3.7709e24*(log10(x) - 1.66159) * (log10(x) - 1.66159) - 100.0*(log10(x) - 1.0)+1.7099e24;
    	}else{
    		return 4.918e24 * (exp(3.711e-1 * log10(x))) -8.059e24;
    	}

	}else{

		// Ioniz rate
		//return 1.922e-22 * exp(5.45887 * log10(x) + 1.345774) - 1.356355e-17;
		if (x < 10){
        	return 0;
		}
		else{
        	return 1e-13 * (exp(-0.5*( log(x*0.0001) * log(x*0.0001))));
		}
	}
}