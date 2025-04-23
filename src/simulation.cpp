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
 
	//t = 0.48e-3;
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
			//Specie electron_ener("e_energy", -1, 511, fill, temp_matrix);
			Specie electron_ener("e_energy", -1, 9.1093837e-31, fill, temp_matrix); // Initialize the electron energy as a new specie
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

void Simulation::push_time(int int_mode, double& j_l, double& j_r, Eigen::MatrixXd& last_flux){

	// ------- Variable intialization --------

	Convection conv(z_step, z_size, S_hori, ez1, ez2);
	double dt;

	double A_cvt = 0.25;
	double A_dif = 0.25;

    // Initialize necessary matrices

	Eigen::MatrixXd Se = Eigen::MatrixXd::Zero(r_size, z_size);
	e_ener = Eigen::MatrixXd::Zero(r_size, z_size);

	Eigen::MatrixXd electron_fluxes = Eigen::MatrixXd::Zero(r_size, z_size); // To be used by electron energy after

	std::vector<Eigen::MatrixXd> reaction_matrix; // Reaction matrix list

	// Fill the reactions matrices with blanks
	reaction_matrix.resize(chemistries.size(), Eigen::MatrixXd::Zero(r_size, z_size));


	// ------- ALGORITHM ------

	calc_e_energy(e_ener); // Update the electron energy matrix

	calc_reaction_vector(reaction_matrix, e_ener); // Calculate the R_i for all the reactions

	calc_update_species_coeffs(e_ener); // Calculates & Updates the Difusion and Mobility coeffectients for the species

	calc_time_step(dt, A_cvt, A_dif); // Calculates the timestep, dt, dinamically


	//dt = 1e-13;
	//dt = 20e-12;

	// In class or function-scope storage
	std::vector<Eigen::MatrixXd> new_n_all(species.size(), Eigen::MatrixXd::Zero(r_size, z_size));
	std::vector<Eigen::MatrixXd> midfluxn_all(species.size(), Eigen::MatrixXd::Zero(r_size, z_size));

	for (size_t s_idx = 0; s_idx < species.size(); ++s_idx) {

    	Specie& s = species[s_idx];
		// Booleans for the flux cases
		bool is_electron = (&s == &species[0]);
		bool is_energy = (&s == &species.back());
		bool is_gas = (s.get_name() == "Ar"); // this could be cached/indexed too

		// Matrices are alocated
		Eigen::MatrixXd& new_n = new_n_all[s_idx];
		Eigen::MatrixXd& midfluxn = midfluxn_all[s_idx];
		new_n.setZero();
		midfluxn.setZero();

		// Save the species matrices to optimize the calculations
    	const auto& dens = s.get_density();
    	const auto& mob = s.get_mob_coef();
    	const auto& dif = s.get_dif_coef();

    	const auto& elec_dens = species[0].get_density();
		//Eigen::MatrixXd new_new_n = Eigen::MatrixXd::Zero(r_size, z_size); FOR THE THE OTHER INTEGRATION METHODS

		if (int_mode == 0){

			for (int i = 0; i < r_size; i++){
		
				for (int j = grid_init; j < grid_end; ++j){ // changed from 0 and z _size to 1 and z_size - 1

					// Change to calculate outside i,j loop
					Se(i,j) = calc_net_reaction_source(s,i,j, reaction_matrix); // Calculate the Creation term according to all reactions

					double aux_flux = 0;

					if (j == grid_init || j == grid_end - 1){ // Special border cases for fluxes
						
						if (is_electron){
							//std::cout <<Se<<"\n"<<j<<"\n\n";
							//std::cout << mob<<"\n\n";
							//std::cout << dif<<"\n\n";
							// For the electrons

							aux_flux = -0.5 * dens(i,j) * calc_vthermal(s, e_ener(i,j) * 11606. * 2/3); // Base electron wall fluxes
							//std::cout << "e "<< calc_vthermal(s, e_ener(i,j) * 11606. * 2/3)<<"\n";

							for (Specie& ss : species){ // Sum the fluxes of p species
								if (ss.get_charge() >= 1){
									aux_flux = aux_flux + 1 * 0.5 * ss.get_density()(i,j) * calc_vthermal(ss, gas_temp);
								}
							}

							aux_flux = aux_flux * S_hori(i,j);

							//std::cout<< s.get_name() << j<<std::endl;

						} else if (is_energy) {
							
							// For the electron energy density specie

							aux_flux = -2/3 * dens(i,j) * calc_vthermal(s, e_ener(i,j) * 11606. * 2/3); // Base electron wall fluxes

							//std::cout << "e_ene "<< calc_vthermal(s, e_ener(i,j) * 11606. * 2/3)<<"\n";

							for (Specie& ss : species){ // Sum the fluxes of p species
								if (ss.get_charge() >= 1){
									aux_flux = aux_flux + 1 * 0.5 * ss.get_density()(i,j) * calc_vthermal(ss, gas_temp) * secondary_emission_energy;
								}
							}

							aux_flux = aux_flux * S_hori(i,j);

							//std::cout<< s.get_name() << j<< " aaaa"<<std::endl; 

						// For the other species

						} else if (is_gas) {

							//aux_flux = -0.5 * s.get_density()(i,j) * calc_vthermal(s, gas_temp) * S_hori(i,j); // Probably needs to be the sum of the excited and charged species
							
							for (Specie& ss : species){ // Sum the fluxes of excited/charged species to the wall
 								if (ss.get_name() == "Ar_star" || ss.get_name() == "Ar_plus"){
									aux_flux = aux_flux + 1 * 0.5 * ss.get_density()(i,j) * calc_vthermal(ss, gas_temp);
								}
							}

							aux_flux = aux_flux * S_hori(i,j);
							//std::cout <<" gases "<<calc_vthermal(s, gas_temp) << "\n";

						} else{

							aux_flux = -0.5 * dens(i,j) * calc_vthermal(s, gas_temp) * S_hori(i,j);
							//std::cout <<" gases "<<calc_vthermal(s, gas_temp) << "\n";
						}
						
						if (j == grid_init && !is_energy){ // Verificar que nao é feito na energia eletronica
							j_l = j_l + aux_flux * s.get_charge() * 1.6e-19; // Can be optimized if the real value of charge is saved in Specie 
						}

						if (j == grid_end - 1  && !is_energy){

							j_r = j_r + aux_flux * s.get_charge() * 1.6e-19; // Can be optimized if the real value of charge is saved in Specie 
						}
					}

					if (is_energy) { // Calculate the Creation for the electron energy

						double temp_sum = 0;

						for (Chemistry& c : chemistries){ // Cycle through the reactions to calculate (8)
							temp_sum = temp_sum + c.calc_reaction_rate(e_ener(i,j)) * c.reagents[1].get_density()(i,j) * c.react_energy_delta; 
						}

						if (j == grid_end-1){
							Se(i,j) = - last_flux(i,j) * ez2(i,j) * 1.6e-19 - elec_dens(i,j) * temp_sum;
						} else {
							Se(i,j) = - last_flux(i,j) * ez1(i,j) * 1.6e-19 - elec_dens(i,j) * temp_sum;
						}
						
					}

					midfluxn(i,j) =  conv.calcFlux_UNO3(i, j, mob(i,j), dif(i,j), dt, 1, dens, s.get_charge(), grid_init, grid_end);
					new_n(i, j) = dens(i, j) + dt * ((midfluxn(i, j) + aux_flux)/vols(i, j) + Se(i,j));

					if (is_electron){
						//std::cout << midfluxn(i, j) << " "<< aux_flux<<"\n";
						electron_fluxes(i,j) = midfluxn(i, j) + aux_flux; // Save the electron fluxes to calculate the electron energy creation 
					}
				}
			}

		}/* else if (int_mode == 1){

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
		
		}*/

		new_n_all[s_idx] = new_n;
		//s.update_density(new_n); // Probably will need to only update the densities after, as to not make the error of sing the ipdated densities in the wrong timestep
	}

	for (size_t s_idx = 0; s_idx < species.size(); ++s_idx) {
		species[s_idx].update_density(new_n_all[s_idx]);
	}

	last_flux = electron_fluxes;

	t = t+dt; // UPDATE TIME

}

void Simulation::write_dens(std::ofstream& file, int i){
	//std::cout << std::fixed << std::setprecision(21);
	//std::cout<<"Density writen at "<<Pot<<std::endl;
	file << "Time: " << t << "\n";
    file << species[i].get_density() << "\n";
    file << "----\n";  // Separator for the next matrix 
}

void Simulation::write_efield(std::ofstream& file){
	//std::cout << std::fixed << std::setprecision(21);
	//std::cout<<"Density writen at "<<Pot<<std::endl;
	file << "Time: " << t << "\n";
    file << ez1 << "\n";
    file << "----\n";  // Separator for the next matrix 
}

void Simulation::write_e_energy(std::ofstream& file){
	//std::cout << std::fixed << std::setprecision(21);
	//std::cout<<"Density writen at "<<Pot<<std::endl;
	file << "Time: " << t << "\n";
    file << e_ener << "\n";
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

void Simulation::calc_e_energy(Eigen::MatrixXd& e_ener){

	Eigen::MatrixXd ne = species[0].get_density();
	Eigen::MatrixXd e_ne = species.back().get_density();


	for (int i = 0; i < r_size; i++){
	
		for (int j = grid_init; j < grid_end; ++j){

			e_ener(i,j) = e_ne(i,j) / (ne(i,j) + 1e-308);
		}
	}
}

void Simulation::calc_reaction_vector(std::vector<Eigen::MatrixXd>& reaction_matrix, Eigen::MatrixXd e_ener){

	int idx = 0;

	for (Chemistry& c : chemistries){

		Eigen::MatrixXd n0 = c.reagents[0].get_density();
		Eigen::MatrixXd n1 = c.reagents[1].get_density();

		for (int i = 0; i < r_size; i++){
		
			for (int j = grid_init; j < grid_end; ++j){
	
				reaction_matrix[idx](i,j) = c.calc_reaction_rate(e_ener(i,j)) * n0(i,j) * n1(i,j); // Calculate R_i for each reaction

			}
		}

		idx ++;
	}
}

void Simulation::calc_update_species_coeffs(Eigen::MatrixXd e_ener) {

    Eigen::MatrixXd e_mu = Eigen::MatrixXd::Zero(r_size, z_size); // Electron mobility coefficients
    Eigen::MatrixXd e_De = Eigen::MatrixXd::Zero(r_size, z_size); // Electron diffusion coefficients

    std::string first_name = species.front().get_name();
    std::string last_name  = species.back().get_name();

    // Precompute electron mobility and diffusion
    for (int i = 0; i < r_size; ++i) {
        for (int j = grid_init; j < grid_end; ++j) {
            double ener = e_ener(i, j);
            if (ener <= 0.025) {
                e_mu(i, j) = 6.4e25 / gas_dens;
            } else if (ener <= 0.14) {
                e_mu(i, j) = 2.831166363547715e+27 * std::pow(ener, 0.9628859587344575) / gas_dens;
            } else {
                e_mu(i, j) = 4.289500273955298e+25 * std::pow(ener, -1.2096423778714385) / gas_dens;
            }
            e_De(i, j) = e_mu(i, j) * ener * 2.0 / 3.0;
        }
    }

    // Loop over species
    for (Specie& s : species) {
        Eigen::MatrixXd mu = Eigen::MatrixXd::Zero(r_size, z_size);
        Eigen::MatrixXd De = Eigen::MatrixXd::Zero(r_size, z_size);

        const std::string& name = s.get_name();
        int charge = s.get_charge();

        if (name == first_name) {
            s.update_mob_coef(e_mu);
            s.update_dif_coef(e_De);
        } 
        else if (name == last_name) {
            for (int i = 0; i < r_size; ++i) {
                for (int j = grid_init; j < grid_end; ++j) {
                    mu(i, j) = e_mu(i, j) * 5.0 / 3.0;
                    De(i, j) = mu(i, j) * e_ener(i, j) * 2.0 / 3.0;
                }
            }
            s.update_mob_coef(mu);
            s.update_dif_coef(De);
        } 
        else if (charge != 0) {
            mu.setConstant(0.12 / gas_pres);
            s.update_mob_coef(mu);
            s.update_dif_coef(De); // De is already zero
        } 
        else {
            s.update_mob_coef(mu); // mu is already zero
            s.update_dif_coef(De); // De is already zero
        }
    }
}

void Simulation::calc_time_step(double& dt, double A_cvt, double A_dif){

	
    const Eigen::MatrixXd& mu = species[0].get_mob_coef(); // Get the electron and diffusion coeffecients
    const Eigen::MatrixXd& De = species[0].get_dif_coef();

    const double mu_max = mu.maxCoeff();
    const double De_max = De.maxCoeff();
    const double eps = 1e-308;

    const double vz_max = std::max({
        std::abs(mu_max * ez1.maxCoeff()),
        std::abs(mu_max * ez2.maxCoeff()),
        eps
    });

    // this comment still has the dieletric relaxation consideration 
	//dt = std::min({0.125 * z_step/vz_max, 0.25* z_step*z_step/De, 0.5 * 8.85e-12/(mu * ne.maxCoeff()* 1.6e-19)});
	
	if (r_size == 1) { // for 1 dimension

        dt = std::min({
            A_cvt * z_step / (vz_max + eps),
            A_dif * z_step * z_step / (De_max + eps)
        });

        //std::cout << De_max<<"\n";

    } else { // for 2 dimensions

        const double vr_max = std::max(std::abs(mu_max * er.maxCoeff()), eps);

        dt = std::min({
            A_cvt * z_step / (vz_max + eps),
            A_cvt * r_step / (vr_max + eps),
            A_dif * z_step * z_step / (De_max + eps),
            A_dif * r_step * r_step / (De_max + eps)
        });
    }
}

double Simulation::calc_net_reaction_source(const Specie& s, int i, int j, const std::vector<Eigen::MatrixXd> reaction_matrix) {
    double result = 0;
    for (size_t k = 0; k < chemistries.size(); ++k) {
        result += s.react_net[k] * reaction_matrix[k](i, j);
    }
    return result;
}