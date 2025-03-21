#include "specie.hpp"

Specie::Specie() {
    // Default constructor definition
}

Specie::Specie(std::string _name, int _charge, double _mass, int* _react_net, Eigen::MatrixXd& _n) : n(_n){
	name = _name;
	charge = _charge;
	mass = _mass;
	react_net = _react_net;

	qm_ratio = charge/mass * E_CHARGE;

	//n = _n;
}

//GETTERS

std::string Specie::get_name(){
	return name;
}

double Specie::get_mass(){
	return mass;
}

int Specie::get_charge(){
	return charge;
}

double Specie::get_qm_ratio(){
	return qm_ratio;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Specie::get_density(){
	return n;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Specie::get_mob_coef(){
	return mu;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Specie::get_dif_coef(){
	return De;
}

void Specie::update_density(Eigen::MatrixXd& new_n){
	n = new_n;
}

void Specie::update_mob_coef(Eigen::MatrixXd& new_mu){
	mu = new_mu;
}

void Specie::update_dif_coef(Eigen::MatrixXd& new_De){
	De = new_De;
}
