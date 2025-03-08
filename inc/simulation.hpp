#pragma once

#include "specie.hpp"
#include "chemistry.hpp"

#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>

class Simulation{

	public:

		Simulation(std::vector<Specie>&, std::vector<Chemistry>&); // Default constructor

		Simulation(int, double, std::vector<Specie>&, std::vector<Chemistry>&);

		Simulation(int, double, Eigen::VectorXd,std::vector<Specie>&, std::vector<Chemistry>&);

		Simulation(int, int, double, double, std::vector<Specie>&, std::vector<Chemistry>&); // Constructor with only grid and Dieletric 

		Simulation(int, int, double, double, std::string, Eigen::VectorXd&, std::vector<Specie>&, std::vector<Chemistry>&, int, int, double, double, double); // Constructor with grid, Dieletric and Species

		void update_charge_density();

		void push_time(int);

		void write_dens(std::ofstream&);

		void set_grid1D(int, double);
		
		void set_grid2D(int, int, double, double);

		void set_dieletric(Eigen::VectorXd);

		void set_species(std::vector<Specie>);

		void set_geometry(std::string);

		void set_potential(const Eigen::MatrixXd& a);

		void set_Er(const Eigen::MatrixXd& a);

		void set_Ez1(const Eigen::MatrixXd& a);

		void set_Ez2(const Eigen::MatrixXd& a);

		double get_t();

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

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> get_Ez1();

		double calc_vthermal(Specie, double);

	private:

		// Time
		double t;

		// Gas temperature;
		double gas_temp;

		// Grid size
		int r_size;
		int z_size;

		// Grid real step (m)
		double r_step;
		double z_step;

		// Grid gas init and end
		int grid_init;
		int grid_end;

		// Matrices for geometrical stored values 
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S_hori; // Surfaces
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S_vert; // Surfaces
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vols;	// Volumes

		// Dieletric material 
		Eigen::Vector<double, Eigen::Dynamic> eps;

		// Geometry use
		std::string geometry;

		// Species included in the Simulation
		std::vector<Specie>& species;

		// Reactions included in the Simulation
		std::vector<Chemistry>& chemistries;

		// Electron emission energies
		double secondary_emission_energy;

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