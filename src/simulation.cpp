#define _USE_MATH_DEFINES

#include "simulation.hpp"
#include "specie.hpp"

Simulation::Simulation(){	
}

Simulation::Simulation(int n, double dz){
	// IN WORK
	z_size = n;
	z_step = dz;
}

Simulation::Simulation(int n, double dz, Eigen::VectorXd _eps){
	// IN WORK

	z_size = n;
	z_step = dz;
	eps = _eps;
}

Simulation::Simulation(int n, int m, double dr, double dz){
	// IN WORK

	r_size = n;
	z_size = m;
	r_step = dr;
	z_step = dz;
}

Simulation::Simulation(int n, int m, double dr, double dz, std::string _geom, Eigen::VectorXd _eps){

	r_size = n;
	z_size = m;
	r_step = dr;
	z_step = dz;
	eps = _eps;
	geometry = _geom;

	// Dynamically allocate an array of size n

	r_grid = new double[r_size];  
	z_grid = new double[z_size];

	phic = new double*[n];
	phiw = new double*[n];
	phie = new double*[n];
	phin = new double*[n];
	phis = new double*[n];

	for (int i = 0; i < n; ++i) {

		// Create new arrays for the 2D phi arrays

        phic[i] = new double[m];
		phiw[i] = new double[m];
		phie[i] = new double[m];
		phin[i] = new double[m];
		phis[i] = new double[m];

		r_grid[i] = (i + 0.5) * dr; // Set up the real grid values
    }

    // Set up the real grid values

    for (int j = 0; j < m; j++){

    	z_grid[j] = (j+0.5)  * dz; 

	}

	// Simple dummy values to help the dr2 calculation
	double hdz = dz*0.5;
	double hdr = dr*0.5;

	for (int i = 0; i < n; i++) {
		
		double dr2 = ((r_grid[i] + hdr) * (r_grid[i] + hdr) - (r_grid[i] - hdr) * (r_grid[i] - hdr));

		for (int j = 0; j < m; j++){

			// Calcular volumes e superficies (superiores/N e á direita/E) de celulas i,j

			if (geometry=="cylindrical"){

				S_hori(i,j) = M_PI * dr2;
				S_vert(i,j) = 2 * M_PI * (r_grid[i] + hdr) * dz;
				vols(i,j) = M_PI * dr2 * dz;

			} else if (geometry == "cartesian"){

				S_hori(i,j) = dr;
				S_vert(i,j) = dz;
				vols(i,j) = dz*dr;

			} else {

				std::cout<<"No defined Geometry"<<std::endl;
				break;
			}
			
			// Set up the phi matrices values

			phin[i][j] = 0;
			phis[i][j] = 0;
			phie[i][j] = 0;
			phiw[i][j] = 0;
			phic[i][j] = 0;

			double epsj = eps(j);
			double epsj1 = epsj;

			if (j+1 < m){
				epsj1 = eps(j+1);
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

void Simulation::update_charge_density(){
	// PLACEHOLDER
	rho = Eigen::MatrixXd::Zero(r_size, z_size);
}

void Simulation::set_grid1D(int n, double dz){
	// IN WORK

	z_size = n;
	z_step = dz;
}

void Simulation::set_grid2D(int n, int m, double dr, double dz){
	// IN WORK

	r_size = n;
	z_size = m;
	r_step = dr;
	z_step = dz;
}

void Simulation::set_dieletric(Eigen::VectorXd _eps){
	eps = _eps;
}

void Simulation::set_geometry(std::string _geom){
	geometry = _geom;
}

void Simulation::set_potential(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a){
	pot = a;
}

void Simulation::set_Er(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a){
	er = a;
}

void Simulation::set_Ez1(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a){
	ez1 = a;
}

void Simulation::set_Ez2(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a){
	ez2 = a;
}

// GETTER FUNTIONS

int Simulation::get_r_size(){
	return r_size;
}

int Simulation::get_z_size(){
	return z_size;
}

double Simulation::get_dr(){
	return r_step;
}

double Simulation::get_dz(){
	return z_step;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Simulation::get_s_hori(){
	return S_hori;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Simulation::get_vols(){
	return vols;
}

Eigen::Vector<double, Eigen::Dynamic> Simulation::get_eps(){

	return eps;
}

double** Simulation::get_phie(){
	return phie;
}
double** Simulation::get_phiw(){
	return phiw;
}
double** Simulation::get_phic(){
	return phic;
}
double** Simulation::get_phin(){
	return phin;
}
double** Simulation::get_phis(){
	return phis;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Simulation::get_rho(){
	return rho;
}