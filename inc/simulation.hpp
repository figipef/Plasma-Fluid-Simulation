#pragma once

#include "specie.hpp"

#include <string>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>

class Simulation{

	public:

		Simulation(); // Default constructor

		Simulation(int, double);

		Simulation(int, double, Eigen::VectorXd);

		Simulation(int, int, double, double); // Constructor with only grid and Dieletric

		Simulation(int, int, double, double, std::string, Eigen::VectorXd); // Constructor with only grid and Dieletric

		void update_charge_density();

		void set_grid1D(int, double);
		
		void set_grid2D(int, int, double, double);

		void set_dieletric(Eigen::VectorXd);

		void set_species(Specie);

		void set_geometry(std::string);

		void set_potential(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>);

		void set_Er(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>);

		void set_Ez1(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>);

		void set_Ez2(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>);

		int get_r_size();
		int get_z_size();

		double get_dr();
		double get_dz();

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> get_s_hori();
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> get_vols();


		Eigen::Vector<double, Eigen::Dynamic> get_eps();

		double** get_phie();
		double** get_phiw();
		double** get_phic();
		double** get_phin();
		double** get_phis();

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> get_rho();

	private:

		// Grid size
		int r_size;
		int z_size;

		// Grid real step (m)
		double r_step;
		double z_step;

		// Matrices for geometrical stored values 
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S_hori; // Surfaces
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S_vert; // Surfaces
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vols;	// Volumes

		// Dieletric material 
		Eigen::Vector<double, Eigen::Dynamic> eps;

		// Geometry use
		std::string geometry;

		// Species included in the Simulation
		Specie species;

		// Grid real values

		double *r_grid; 
		double *z_grid;

		// Values to aid the Poisson Solver

		double **phie;
		double **phiw;

		double **phic;

		double **phin;
		double **phis;

		// Updated values 
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rho; // Charge density
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> pot; // Eletric Potential
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> er; // Radial Eletric Field
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ez1; // Horizontal Eletric Field on the left
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ez2; // Horizontal Eletric Field on the right

};