#include "poissonsolver2d.hpp"
#include "simulation.hpp"

PoissonSolver2D::PoissonSolver2D(double V0, double VMAX, double VWALL, double VIN ,double fronteira[4], Eigen::VectorXd& _sig, Simulation& _simul) : simul(_simul){

	// Set up all variables necessary from simul

	sig = _sig;

	int r_size = simul.get_r_size();
	int z_size = simul.get_z_size(); 

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S_hori = simul.get_s_hori();

	Eigen::Vector<double, Eigen::Dynamic> eps = simul.get_eps();

	double **phie = _simul.get_phie();
	double **phiw = _simul.get_phiw();
	double **phic = _simul.get_phic();
	double **phin = _simul.get_phin();
	double **phis = _simul.get_phis();

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

	// Create the matrix that's needed for the calcualtion

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

	solver.analyzePattern(PHI);

    solver.factorize(PHI); // Set up the factorization so that solving the system is faster

	auto t2 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> d2 = t2 - t1;
    
    std::cout << "Time taken to Load PHI's (coeffecients): " << d2.count() << " ms" << std::endl;
}

void PoissonSolver2D::solve(){

	simul.update_charge_density(); // Update the charge density values

	// get the necessary parameters from the simulation

	int r_size = simul.get_r_size();
	int z_size = simul.get_z_size();

	double r_step = simul.get_dr();
	double z_step = simul.get_dz();

	Eigen::Vector<double, Eigen::Dynamic> eps = simul.get_eps();
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rho = simul.get_rho();
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vols = simul.get_vols();

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RHS = Eigen::MatrixXd::Zero(r_size, z_size); // Right hand side of the equation

	for (int i = 0; i < r_size; i++){
		for (int j = 0; j < z_size; j++){

			RHS(i,j) = rhs_i(i,j) + rho(i,j) * vols(i,j); // calculate the right hand side
		}
	}

	Eigen::VectorXd b(RHS.size());
    b = Eigen::Map<const Eigen::VectorXd>(RHS.data(), RHS.size());

    Eigen::VectorXd v = solver.solve(b); // Just sovle the matrix system

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

				// Calculate the Eletric field

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

	// Save the values

	simul.set_potential(V);
	simul.set_Er(E_r); 
	simul.set_Ez1(E1_z);
	simul.set_Ez2(E2_z);
}