#include "specie.hpp"

specie::specie(std::string _name, int _charge, double _mass, int* _react_net){
	name = _name;
	charge = _charge;
	mass = _mass;
	react_net = _react_net;

	qm_ratio = charge/mass * E_CHARGE;
}

std::string specie::get_name(){
	return name;
}

double specie::get_mass(){
	return mass;
}

double specie::get_charge(){
	return charge;

}

double specie::get_qm_ratio(){
	return qm_ratio;
}