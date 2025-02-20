#pragma once

#include "specie.hpp"

#include <string>
#include <Eigen/Dense>

class Simulation{

	public:

		Simulation(); // Default constructor

		Simulation(int, double);

		Simulation(int, double, Eigen::VectorXd);

		Simulation(int, int, double, double); // Constructor with only grid and Dieletric

		Simulation(int, int, double, double, Eigen::VectorXd); // Constructor with only grid and Dieletric

		void set_grid1D(int, double);
		
		void set_grid2D(int, int, double, double);

		void set_dieletric(Eigen::VectorXd);

		void set_species(Specie);

		void set_geometry(std::string);


	private:

		int r_size;
		int z_size;

		double r_step;
		double z_step;

		Eigen::Vector<double, Eigen::Dynamic> eps;

		std::string geometry;

		Specie species;

};