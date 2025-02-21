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

		PoissonSolver2D(double V0, double VMAX, double VWALL, double VIN ,double[4], Eigen::VectorXd, Simulation);

		void solve();

	private:

		// Variables for the initialization of the Solver
		Simulation simul;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rhs_i;
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		Eigen::VectorXd sig;

		//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	
};