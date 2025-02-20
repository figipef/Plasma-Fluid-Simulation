#include "specie.hpp"

Specie::Specie(std::string _name, int _charge, double _mass, int* _react_net){
	name = _name;
	charge = _charge;
	mass = _mass;
	react_net = _react_net;

	qm_ratio = charge/mass * E_CHARGE;
}

std::string Specie::get_name(){
	return name;
}

double Specie::get_mass(){
	return mass;
}

double Specie::get_charge(){
	return charge;

}

double Specie::get_qm_ratio(){
	return qm_ratio;
}