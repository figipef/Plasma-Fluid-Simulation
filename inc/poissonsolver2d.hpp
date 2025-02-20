#pragma once

#include <Eigen/Dense>
#include <chrono>

class PoissonSolver2D{

	public:
	
		PoissonSolver2D();

		PoissonSolver2D(double V0, double VMAX, double VWALL, double VIN ,double[4], Eigen::VectorXd, Eigen::VectorXd);

	private:

		double rhs_i;
	
};