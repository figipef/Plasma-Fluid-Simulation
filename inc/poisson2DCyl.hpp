#pragma once

#include <map>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

class Poisson2DCyl{

	public:
		double r_step; // dr and dz
		double z_step;

		double v0;
		double vmax;
		double vwall;
		double vin;
		double front[4];
		std::string name;

		int r_size; // n x m
		int z_size;
		
		double *r_grid; // Grid real values
		double *z_grid;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S_hori;
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S_vert;
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vols;


		std::vector<std::pair<int, double>> eps_vec; // Dieletric constant across along z
		std::vector<std::pair<int, double>> sig_vec; // Surface charge densety along z

		//double **v;	
		
		double **phie;
		double **phiw;

		double **phic;

		double **phin;
		double **phis;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rhs_i;
		Eigen::SparseMatrix<double> phi;
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Pot;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Er;
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Ez1;
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Ez2;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ne;
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ni;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rho;
		
		Poisson2DCyl(int, int, double, double, std::vector<std::pair<int, double>>, std::vector<std::pair<int, double>>); // n, m, rstep, zstep eps, sig

		void solve(double, double, double, double ,double (*)(double, double), double[4], std::string);
		void solve(double, double, double, double ,Eigen::MatrixXd, Eigen::MatrixXd, double[4], std::string);

		void solve_Poisson();

		void push_time(double, double, std::ofstream&);

		void calculate_charge_density();

		void write_fields();
};
