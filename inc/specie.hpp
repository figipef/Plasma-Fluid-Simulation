#pragma once

#include <string>

double const E_CHARGE = 1.6e-19;

class specie{

	public:

		specie(std::string, int, double, int*);

		std::string get_name();

		double get_mass();

		double get_charge();

		double get_qm_ratio();

	private:

		std::string name;

		int charge;

		double mass;

		int *react_net; // net value for each reaction in the array form

		double qm_ratio;
};