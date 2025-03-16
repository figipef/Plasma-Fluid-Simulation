#define _USE_MATH_DEFINES

#include "simulation.hpp"
#include "convection.hpp"
#include "specie.hpp"

Simulation::Simulation(std::vector<Specie>& _species, std::vector<Chemistry>& _chemistries) : species(_species), chemistries(_chemistries){	
}

Simulation::Simulation(int n, double dz, std::vector<Specie>& _species, std::vector<Chemistry>& _chemistries) : species(_species), chemistries(_chemistries){
	// IN WORK
	z_size = n;
	z_step = dz;
}

Simulation::Simulation(int n, double dz, Eigen::VectorXd _eps, std::vector<Specie>& _species, std::vector<Chemistry>& _chemistries) : species(_species), chemistries(_chemistries){
	// IN WORK

	z_size = n;
	z_step = dz;
	eps = _eps;
}

Simulation::Simulation(int n, int m, double dr, double dz, std::vector<Specie>& _species, std::vector<Chemistry>& _chemistries) : species(_species), chemistries(_chemistries){
	// IN WORK

	r_size = n;
	z_size = m;
	r_step = dr;
	z_step = dz;
}

Simulation::Simulation(int n, int m, double dr, double dz, std::string _geom, Eigen::VectorXd& _eps, std::vector<Specie>& _species, std::vector<Chemistry>& _chemistries, int _grid_init, int _grind_end, double electron_energy, double sec_e_em_energy, double _gas_temp, double _gas_dens, double _gas_pres) : species(_species), chemistries(_chemistries), gas_dens(_gas_dens), gas_pres(_gas_pres){
 
	t = 0;
	gas_temp = _gas_temp;

	secondary_emission_energy = sec_e_em_energy;

	r_size = n;
	z_size = m;
	r_step = dr;
	z_step = dz;
	eps = _eps;
	geometry = _geom;

	grid_init = _grid_init;
	grid_end = _grind_end;

	for (Specie s : species){
		int fill[] = {0};
		if (s.get_name() == "e"){
			Eigen::MatrixXd temp_matrix = s.get_density() * electron_energy;
			Specie electron_ener("e_energy", -1, 511, fill, temp_matrix);
			species.push_back(electron_ener);
		}
	}

	
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

void Simulation::push_time(int int_mode, double& j_l, double& j_r){

	Convection conv(z_step, z_size, S_hori, ez1);
	double dt;

	// To calculate the electron temperature
	const double k_B = 8.617333e-5;  // Boltzmann constant in eV/K
    const double e = 1.602176634e-19; // Elementary charge in C

	Eigen::MatrixXd mu = Eigen::MatrixXd::Zero(r_size, z_size);
	Eigen::MatrixXd De = Eigen::MatrixXd::Zero(r_size, z_size);

	Eigen::MatrixXd mu_energy = Eigen::MatrixXd::Zero(r_size, z_size);
	Eigen::MatrixXd De_energy = Eigen::MatrixXd::Zero(r_size, z_size);

	Eigen::MatrixXd mu_gas = Eigen::MatrixXd::Zero(r_size, z_size);
	Eigen::MatrixXd De_gas = Eigen::MatrixXd::Zero(r_size, z_size);

	Eigen::MatrixXd Se = Eigen::MatrixXd::Zero(r_size, z_size);
	Eigen::MatrixXd e_ener = Eigen::MatrixXd::Zero(r_size, z_size);

	Eigen::MatrixXd electron_fluxes = Eigen::MatrixXd::Zero(r_size, z_size); // To be used by electron energy after

	std::vector<Eigen::MatrixXd> reaction_matrix; // Reaction matrix list

	//std::vector<DataPoint> dif_data = { // Difusion data
    //    {0.001678, 1.092e+25}, {0.002221, 2.196e+25}, {0.003451, 3.298e+25}, 
    //    {0.007177, 3.861e+25}, {0.01671, 3.777e+25}, {0.03987, 3.272e+25}, 
    //    {0.09591, 2.646e+25}, {0.2593, 2.074e+25}, {1.902, 1.607e+25}, 
    //    {14.64, 1.23e+25}, {74.43, 8.949e+24}, {260.3, 6.922e+24}, 
    //    {647.5, 6.551e+24}, {1240.0, 7.629e+24}, {2032.0, 9.983e+24}, 
    //    {3075.0, 1.362e+25}, {4442.0, 1.893e+25}, {6182.0, 2.686e+25}
    //};

	// Fill the reactions matrices
	for (Chemistry& c : chemistries){ 
		Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(r_size, z_size);
		reaction_matrix.push_back(mat);
	}

	for (int i = 0; i < r_size; i++){
	
		for (int j = grid_init; j < grid_end; ++j){

			e_ener(i,j) = species.back().get_density()(i,j) / (species[0].get_density()(i,j) + 1e-308);

			// Set the values of the reaction matrix using the electron energy
			int counter = 0;
			for (Chemistry& c : chemistries){ // Cycle through the reactions to calculate (8)

				reaction_matrix[counter](i,j) = c.calc_reaction_rate(e_ener(i,j)) * c.reagents[0].get_density()(i,j) * c.reagents[1].get_density()(i,j); // Calculate R_i for each reaction
				counter++;
			}

			double x;
			
			if (j == grid_init){
				x = std::abs(ez2(i,j) / (gas_dens) * 1e21);
			}else if (j == grid_end - 1){
				x = std::abs(ez1(i,j-1) / (gas_dens) * 1e21);
			}else{
				x = std::abs((ez1(i,j) + ez1(i,j-1))/2 / (gas_dens) * 1e21);
			}
			//std::cout <<x<<"\n";
			/* DIFFUSION COEFFECIENT CALCULATION THROUGH LIN INTERPOLATION
			double dif_coef = 0;

			// Find the two data points surrounding x
    		for (size_t i = 0; i < dif_data.size() - 1; ++i) {
    		    if (x >= dif_data[i].x && x <= dif_data[i + 1].x) {
    		        // Interpolate between data[i] and data[i + 1]
    		        double x0 = dif_data[i].x;
    		        double y0 = dif_data[i].y;
    		        double x1 = dif_data[i + 1].x;
    		        double y1 = dif_data[i + 1].y;
            
    		        // Perform the linear interpolation
    		        dif_coef = y0 + (x - x0) * (y1 - y0) / (x1 - x0);
    		    } else if (x <= dif_data[0].x){
    		    	dif_coef = 0;
    		    } else if (x >= dif_data[dif_data.size() - 1].x){
    		    	dif_coef = 0;
    		    }
    		}
			*/ 
			//mu(i,j) = Pol(x, 0)/2.43e25;
			mu(i,j) = (-1.0475e26* std::atan(std::log10(x) * 2.59783 - 4.66191) + 1.5733318e26)/gas_dens;
			De(i,j) = mu(i,j) * e_ener(i,j) * 2/3;

			mu_energy(i,j) = mu(i,j) * 5/3;
			De_energy(i,j) = mu_energy(i,j) * e_ener(i,j) * 2/3;

			mu_gas(i,j) = 0.12/gas_pres;
			//De(i,j) = dif_coef/gas_dens; 

			//Se(i,j) = Pol(x, 2) * 2.43e25 * ne(i,j); 
			
			//mu(i,j) = 0.03;
			//De(i,j) = 0.1;
		}
	}
	//std::cout <<e_ener<< std::endl;
	std::cout <<mu<<"\n";

	std::cout <<De<<"\n";
	// Calculate Fluxes at each cell

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

	//dt = 1e-13;
	//dt = 20e-12;

	int counter = 0; // FOR TESTING PURPOSES REMOVE WHEN DONE, only applied to the first integration mode

	for (Specie& s : species){

		counter++;

		if (s.get_name() == species[0].get_name()){
			s.update_mob_coef(mu);
			s.update_dif_coef(De);
		} else if (s.get_name() == species.back().get_name()) {
			s.update_mob_coef(mu_energy);
			s.update_dif_coef(De_energy);
		} else {
			s.update_mob_coef(mu_gas);
			s.update_dif_coef(De_gas);
		}

		Eigen::MatrixXd new_n = Eigen::MatrixXd::Zero(r_size, z_size);
		Eigen::MatrixXd new_new_n = Eigen::MatrixXd::Zero(r_size, z_size);
		Eigen::MatrixXd midfluxn = Eigen::MatrixXd::Zero(r_size, z_size);

		if (int_mode == 0){

			for (int i = 0; i < r_size; i++){
		
				for (int j = grid_init; j < grid_end; ++j){ // changed from 0 and z _size to 1 and z_size - 1

					// Calculate the Reactions at the cell i,j
					double temp_S = 0; // temporary sum dummy variable
					int counter = 0;
					for (Chemistry& c : chemistries){ // Cycle through the reactions to calculate (8)

						temp_S = temp_S + s.react_net[counter] * reaction_matrix[counter](i,j); // Calculate R_i for each reaction
						counter++;
					}

					Se(i,j) = temp_S;

					double aux_flux = 0;

					if (j == grid_init || j == grid_end - 1){ // Special border cases for fluxes

						if (s.get_name() == species[0].get_name()){

							// For the electrons

							aux_flux = 0.5 * s.get_density()(i,j) * calc_vthermal(s, e_ener(i,j) * 11606. * 2/3); // Base electron wall fluxes

							for (Specie& ss : species){ // Sum the fluxes of p species
								if (ss.get_charge() >= 1){
									aux_flux = aux_flux - 1 * 0.5 * ss.get_density()(i,j) * calc_vthermal(ss, gas_temp);
								}
							}

							aux_flux = aux_flux * S_hori(i,j);

							std::cout<< s.get_name() << j<<std::endl;

						} else if (s.get_name() == species.back().get_name()) {
							
							// For the electron energy density specie

							aux_flux = 2/3 * s.get_density()(i,j) * calc_vthermal(s, e_ener(i,j) * 11606. * 2/3); // Base electron wall fluxes

							for (Specie& ss : species){ // Sum the fluxes of p species
								if (ss.get_charge() >= 1){
									aux_flux = aux_flux - 1 * 0.5 * ss.get_density()(i,j) * calc_vthermal(ss, gas_temp) * secondary_emission_energy;
								}
							}

							aux_flux = aux_flux * S_hori(i,j);

							std::cout<< s.get_name() << j<< " aaaa"<<std::endl; 

						// For the other species

						} else if (s.get_name() == "Ar") {

							aux_flux = -0.5 * s.get_density()(i,j) * calc_vthermal(s, gas_temp) * S_hori(i,j); // Probably needs to be the sum of the excited and charged species

						} else{

							aux_flux = 0.5 * s.get_density()(i,j) * calc_vthermal(s, gas_temp) * S_hori(i,j);
						}

						if (j == grid_init){ // Verificar que nao é feito na energia eletronica
							j_l = j_l + aux_flux * s.get_charge() * 1.6e-19; // Can be optimized if the real value of charge is saved in Specie 
						}

						if (j == grid_end - 1){
							j_r = j_r + aux_flux * s.get_charge() * 1.6e-19; // Can be optimized if the real value of charge is saved in Specie 
						}
					}

					if (s.get_name() == species.back().get_name()) { // Calculate the Creation for the electron energy

						double temp_sum = 0;

						for (Chemistry& c : chemistries){ // Cycle through the reactions to calculate (8)
							temp_sum = temp_sum + c.calc_reaction_rate(e_ener(i,j)) * c.reagents[1].get_density()(i,j) * c.react_energy_delta; 
						}

						Se(i,j) = - electron_fluxes(i,j) * ez1(i,j) * 1.6e-19 - species[0].get_density()(i,j) * temp_sum;
					}

					midfluxn(i,j) =  conv.calcFlux_UNO3(i, j, s.get_mob_coef()(i,j), s.get_dif_coef()(i,j), dt, 1, s.get_density(), s.get_charge());
					new_n(i, j) = s.get_density()(i, j) + dt * ((midfluxn(i, j) + aux_flux)/vols(i, j) + Se(i,j));

					if (s.get_name() == species[0].get_name()){
						electron_fluxes(i,j) = midfluxn(i, j) + aux_flux; // Save the electron fluxes to calculate the electron energy creation 
					}
				}
			}

		} else if (int_mode == 1){

			for (int i = 0; i < r_size; i++){
		
				for (int j = grid_init; j < grid_end; ++j){ 
		
					midfluxn(i,j) = conv.calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, 1, s.get_density(), s.get_charge());
		
					new_new_n(i, j) = s.get_density()(i, j) + dt * midfluxn(i, j)/vols(i, j);


				}
			}

			for (int i = 0; i < r_size; i++){
		
				for (int j = grid_init; j < grid_end; ++j){ // ter atencao n(i,j) * mu
		
					new_n(i,j) = s.get_density()(i, j) + dt/2*(midfluxn(i, j) + conv.calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, 1, new_new_n, s.get_charge()))/vols(i, j);

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

				for (int j = grid_init; j < grid_end; ++j){ 
		
					k1(i,j) = conv.calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, 1, s.get_density(), s.get_charge());
		
				}
			}
		
			helpk1 = s.get_density() + k1 * dt/2.;
		
			for (int i = 0; i < r_size; i++){
		
				for (int j = grid_init; j < grid_end; ++j){ 
		
					k2(i,j) = conv.calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, 1, helpk1, s.get_charge());
		
				}
			}
		
			helpk2 = s.get_density() + k2 * dt/2.;
		
			for (int i = 0; i < r_size; i++){
		
				for (int j = grid_init; j < grid_end; ++j){ 
		
					k3(i,j) = conv.calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, 1, helpk2, s.get_charge());
		
				}
			}
		
			helpk3 = s.get_density() + k3 * dt; // Atencao, provavelmente adicionar /vols() ????
		
			for (int i = 0; i < r_size; i++){
		
				for (int j = grid_init; j < grid_end; ++j){ 
		
					k4(i,j) = conv.calcFlux_Koren(i, j, mu(i,j), De(i,j), dt, 1, helpk3, s.get_charge());
				}
			}
			
		
			for (int i = 0; i < r_size; i++){
		
				for (int j = grid_init; j < grid_end; ++j){ 
		
					new_n(i,j) = s.get_density()(i,j) + dt*(((k1(i,j) + 2*k2(i,j) + 2*k3(i,j) + k4(i,j))/vols(i, j))/6 + Se(i,j));
				}
			}
		
		}
		//std::cout<<new_n<<std::endl;
		s.update_density(new_n);
	}

	//ne = new_n;
	//ni = new_ni;

	t = t+dt; // UPDATE TIME
	std::cout <<dt<<std::endl;
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

double Simulation::calc_vthermal(Specie specie, double temp){
	return sqrt(8 * 1.380649e-23 * temp / (M_PI * specie.get_mass()));
}

