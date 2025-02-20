#include "poissonsolver2d.hpp"

PoissonSolver2D::PoissonSolver2D(double V0, double VMAX, double VWALL, double VIN ,double fronteira[4], Eigen::VectorXd eps, Eigen::VectorXd sig){

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

