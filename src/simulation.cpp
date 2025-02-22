#define _USE_MATH_DEFINES

#include "simulation.hpp"
#include "convection.hpp"
#include "specie.hpp"

Simulation::Simulation(std::vector<Specie>& _species) : species(_species){	
}

Simulation::Simulation(int n, double dz, std::vector<Specie>& _species) : species(_species){
	// IN WORK
	z_size = n;
	z_step = dz;
}

Simulation::Simulation(int n, double dz, Eigen::VectorXd _eps, std::vector<Specie>& _species) : species(_species){
	// IN WORK

	z_size = n;
	z_step = dz;
	eps = _eps;
}

Simulation::Simulation(int n, int m, double dr, double dz, std::vector<Specie>& _species) : species(_species){
	// IN WORK

	r_size = n;
	z_size = m;
	r_step = dr;
	z_step = dz;
}

Simulation::Simulation(int n, int m, double dr, double dz, std::string _geom, Eigen::VectorXd& _eps, std::vector<Specie>& _species) : species(_species){
 
	t = 0;

	r_size = n;
	z_size = m;
	r_step = dr;
	z_step = dz;
	eps = _eps;
	geometry = _geom;
	//species = _species;

	S_hori = Eigen::MatrixXd::Zero(r_size, z_size);
	S_vert = Eigen::MatrixXd::Zero(r_size, z_size);
	vols = Eigen::MatrixXd::Zero(r_size, z_size);
	rho = Eigen::MatrixXd::Zero(r_size, z_size);

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
	std::cout <<"Initalized Class Simulation"<<std::endl;
}

void Simulation::update_charge_density(){
	// PLACEHOLDER
	rho = Eigen::MatrixXd::Zero(r_size, z_size);

	for (Specie s : species){
		rho = rho + s.get_density() * s.get_charge() * 1.6e-19; // Can be optimized if the real value of charge is saved in Specie
	}
}

void Simulation::push_time(int int_mode){

	//double mu = 0.03; // To be changed later !!!
	//double De = 0.1;
	Convection conv(z_step, z_size, S_hori, ez1);
	double dt;

	Eigen::MatrixXd mu = Eigen::MatrixXd::Zero(r_size, z_size);
	Eigen::MatrixXd De = Eigen::MatrixXd::Zero(r_size, z_size);
	Eigen::MatrixXd Se = Eigen::MatrixXd::Zero(r_size, z_size);

	//std::cout <<Ez1<<std::endl;

	for (int i = 0; i < r_size; i++){
	
		for (int j = 0; j < z_size; ++j){ 
			/*
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
			*/
			mu(i,j) = 0.03;
			De(i,j) = 0.1;
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
		vr_max = std::abs(mu.maxCoeff() * er.maxCoeff());
	}

	double vz_max = std::max({std::abs(mu.maxCoeff() * ez1.maxCoeff()), std::abs(mu.maxCoeff() * ez2.maxCoeff()), 1e-308});

	//dt = std::min({0.125 * z_step/vz_max, 0.25* z_step*z_step/De, 0.5 * 8.85e-12/(mu * ne.maxCoeff()* 1.6e-19)});
	
	double A_cvt = 0.1;
	double A_dif = 0.1;

	if (r_size == 1){
	
		dt = std::min({A_cvt * z_step/(vz_max + 1e-308), A_dif* z_step*z_step/(De.maxCoeff() + 1e-308)});
	
	} else {
	
		dt = std::min({A_cvt * z_step/vz_max, A_cvt * r_step/vr_max, A_dif* z_step*z_step/De.maxCoeff(), A_dif * r_step*r_step/De.maxCoeff()});
	
	}
	dt = 1e-13;
	dt = 20e-12;

	int counter = 0; // FOR TESTING PURPOSES REMOVE WHEN DONE, only applied to the first integration mode

	for (Specie& s : species){

		counter++;

		Eigen::MatrixXd new_n = Eigen::MatrixXd::Zero(r_size, z_size);
		Eigen::MatrixXd new_new_n = Eigen::MatrixXd::Zero(r_size, z_size);
		Eigen::MatrixXd midfluxn = Eigen::MatrixXd::Zero(r_size, z_size);

		if (int_mode == 0){

			for (int i = 0; i < r_size; i++){
		
				for (int j = 0; j < z_size; ++j){ // changed from 0 and z _size to 1 and z_size - 1
					
					if (counter == 1){
						midfluxn(i,j) =  conv.calcFlux_UNO3(i, j, mu(i,j), De(i,j), dt, 1, s.get_density());
						new_n(i, j) = s.get_density()(i, j) + dt * (midfluxn(i, j)/vols(i, j) + Se(i,j));
					} else {
						new_n(i, j) = s.get_density()(i, j) + dt * Se(i,j);
					}
		
					//midfluxni(i,j) = calcFlux(i, j, 0, 0, dt, ni, alph);
		
					//new_ne(i, j) = ne(i, j) + dt * (midfluxn(i, j)/vols(i, j) + Se(i,j));
		
					//new_ni(i, j) = ni(i, j) + dt * Se(i,j);//midfluxni(i, j);

				}
			}

		} else if (int_mode == 1){

			for (int i = 0; i < r_size; i++){
		
				for (int j = 0; j < z_size; ++j){ 
		
					midfluxn(i,j) = conv.calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, 1, s.get_density());
		
					//midfluxni(i,j) = calcFlux_superbee(i, j, 0, 0, dt, ni, alph);
		
					new_new_n(i, j) = s.get_density()(i, j) + dt * midfluxn(i, j)/vols(i, j);
		
					//new_new_ni(i, j) = ni(i, j) + dt * midfluxni(i, j);

				}
			}

			for (int i = 0; i < r_size; i++){
		
				for (int j = 0; j < z_size; ++j){ // ter atencao n(i,j) * mu
		
					new_n(i,j) = s.get_density()(i, j) + dt/2*(midfluxn(i, j) + conv.calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, 1, new_new_n))/vols(i, j);

					//new_ni(i,j) = ni(i, j) + dt/2*(midfluxni(i, j) + calcFlux_superbee(i, j,  0,  0, dt, new_new_ni, alph))/vols(i, j);

					//new_ni(i, j) = ni(i, j) + dt * Se(i,j);//midfluxni(i, j); NEWLY CHANGED
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
		
					k1(i,j) = conv.calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, 1, s.get_density());
		
				}
			}
		
			helpk1 = s.get_density() + k1 * dt/2.;
		
			for (int i = 0; i < r_size; i++){
		
				for (int j = 0; j < z_size; ++j){ 
		
					k2(i,j) = conv.calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, 1, helpk1);
		
				}
			}
		
			helpk2 = s.get_density() + k2 * dt/2.;
		
			for (int i = 0; i < r_size; i++){
		
				for (int j = 0; j < z_size; ++j){ 
		
					k3(i,j) = conv.calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, 1, helpk2);
		
				}
			}
		
			helpk3 = s.get_density() + k3 * dt; // Atencao, provavelmente adicionar /vols() ????
		
			for (int i = 0; i < r_size; i++){
		
				for (int j = 0; j < z_size; ++j){ 
		
					k4(i,j) = conv.calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, 1, helpk3);
				}
			}
			
		
			for (int i = 0; i < r_size; i++){
		
				for (int j = 0; j < z_size; ++j){ 
		
					new_n(i,j) = s.get_density()(i,j) + dt*(((k1(i,j) + 2*k2(i,j) + 2*k3(i,j) + k4(i,j))/vols(i, j))/6 + Se(i,j));

					// new_ni(i, j) = ni(i, j) + dt * Se(i,j); NEWLY CHANGED
				}
			}
		
		}
		//std::cout<<new_n<<std::endl;
		s.update_density(new_n);
	}

	//ne = new_n;
	//ni = new_ni;

	t = t+dt; // UPDATE TIME

	//solve_Poisson();

	//write_dt(dt_file, dt);

	//if ( ti % 999 == 0){
		//std::cout <<De<<std::endl;
		//std::cout <<Se*dt<<std::endl;

	//std::cout << " t = "<<t<<std::endl;
	//std::cout << " dt =" <<dt<<std::endl;
	//}	
}

void Simulation::write_dens(std::ofstream& file){
	//std::cout << std::fixed << std::setprecision(21);
	//std::cout<<"Density writen at "<<Pot<<std::endl;
	file << "Time: " << t << "\n";
    file << species[0].get_density() << "\n";
    file << "----\n";  // Separator for the next matrix 
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

void Simulation::set_species(std::vector<Specie> _species){
	species = _species;
}

void Simulation::set_dieletric(Eigen::VectorXd _eps){
	eps = _eps;
}

void Simulation::set_geometry(std::string _geom){
	geometry = _geom;
}

void Simulation::set_potential(const Eigen::MatrixXd& a){
	pot = a; 
}

void Simulation::set_Er(const Eigen::MatrixXd& a){
	er = a;
}

void Simulation::set_Ez1(const Eigen::MatrixXd& a){
	//std::cout<<a<<std::endl;
	ez1 = a;
	//std::cout<<ez1<<std::endl;
}

void Simulation::set_Ez2(const Eigen::MatrixXd& a){
	ez2 = a;
}

// GETTER FUNTIONS

double Simulation::get_t(){
	return t;
}

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

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Simulation::get_Ez1(){
	return ez1;
}