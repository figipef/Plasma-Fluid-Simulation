#pragma once

#include "simulation.hpp"

#include <vector>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

class PoissonSolver2D{

	public:
	
		PoissonSolver2D();

		PoissonSolver2D(double[4], Eigen::VectorXd&, Simulation&); 

		void update_boundary_voltage(double[4], double, double, double, double);

		void solve();

	private:

		// Variables for the initialization of the Solver
		Simulation& simul;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rhs_i;  // Right hand side only considering the PHI's
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rhs_ii; // Right hand-side considering the boundary voltages
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		Eigen::VectorXd sig;

		//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	
};