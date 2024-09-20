#pragma once

#include <map>
#include <string>
#include <vector>

class Poisson2DCyl{

	public:
		double r_step; // dr and dz
		double z_step;

		int r_size; // n x m
		int z_size;
		
		double *r_grid; // Grid real values
		double *z_grid;


		std::vector<std::pair<int, double>> eps_vec; // Dieletric constant across along z
		std::map<int, double> sig_map; // Surface charge densety along z

		double **v;	
		
		double **phie;
		double **phiw;

		double **phic;

		double **phin;
		double **phis;

		Poisson2DCyl(int, int, double, double, std::vector<std::pair<int, double>>, std::map<int, double>); // n, m, rstep, zstep eps, sig

		void solve(double, double, double, double ,double (*)(double, double), double[4], std::string);

};
