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
		double dh = 1;

		double t;

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

		double pol_degree;
		Eigen::Vector<double, Eigen::Dynamic> mob_pol_vec;
		Eigen::Vector<double, Eigen::Dynamic> temp_pol_vec;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S_hori;
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S_vert;
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vols;


		std::vector<std::pair<int, double>> eps_vec; // Dieletric constant across along z
		std::vector<std::pair<int, double>> sig_vec; // Surface charge densety along z

		Eigen::Vector<double, Eigen::Dynamic> eps;
		Eigen::Vector<double, Eigen::Dynamic> sig;

		//double **v;	
		
		double **phie;
		double **phiw;

		double **phic;

		double **phin;
		double **phis;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rhs_i;
		Eigen::SparseMatrix<double> phi;
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Pot;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Er;
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Ez1;
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Ez2;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ne;
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ni;

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rho;
		
		Poisson2DCyl(int, int, double, double, Eigen::VectorXd, Eigen::VectorXd); // n, m, rstep, zstep eps, sig

		//void solve(double, double, double, double ,double (*)(double, double), double[4]);

		void set_pol_coeffs(Eigen::VectorXd, Eigen::VectorXd, int);

		void solve(double, double, double, double ,Eigen::MatrixXd, Eigen::MatrixXd, double[4]);

		void solve_Poisson();

		void push_time(int, double,std::ofstream&, int);

		void calculate_charge_density();

		void write_fields(std::string);

		void write_dens(std::ofstream&);

		void write_dt(std::ofstream&, double);

	private:
		
		int getSign(double, double);

		void uno2ConvectionScheme(Eigen::MatrixXd&, Eigen::MatrixXd&, double, Eigen::MatrixXd, Eigen::MatrixXd);

		double midWayFlux(double, double, double, double, double);

		double calculate_g_c(double, double, double, double, double);
		
		double calculate_g_c_UNO2(double, double);

		double phia(double);

		double korenLimiter(double);

		double calcFlux(int, int, double, double, double, Eigen::MatrixXd&, double);

		double calcFlux_superbee(int, int, double, double, double, Eigen::MatrixXd&);

		double calcFlux_UNO3(int , int ,double, double, double, Eigen::MatrixXd&);

		double calcFlux_Koren(int, int, double, double, double, Eigen::MatrixXd&);

		double Pol(double, int);

};
