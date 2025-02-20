#include "simulation.hpp"

Simulation::Simulation(){}

Simulation::Simulation(int n, double dz){

	z_size = n;
	z_step = dz;

}

Simulation::Simulation(int n, double dz, Eigen::VectorXd _eps){

	z_size = n;
	z_step = dz;
	eps = _eps;

}

Simulation::Simulation(int n, int m, double dr, double dz){

	r_size = n;
	z_size = m;
	r_step = dr;
	z_step = dz;

}

Simulation::Simulation(int n, int m, double dr, double dz, Eigen::VectorXd _eps){

	r_size = n;
	z_size = m;
	r_step = dr;
	z_step = dz;
	eps = _eps;

}

void Simulation::set_grid1D(int n, double dz){

	z_size = n;
	z_step = dz;

}

void Simulation::set_grid2D(int n, int m, double dr, double dz){

	r_size = n;
	z_size = m;
	r_step = dr;
	z_step = dz;

}

void Simulation::set_dieletric(Eigen::VectorXd _eps){
	eps = _eps;
}

void Simulation::set_geometry(std::string _geom){
	geometry = _geom;
}