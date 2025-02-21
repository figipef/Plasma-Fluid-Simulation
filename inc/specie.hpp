#pragma once

#include <string>
#include <Eigen/Dense>

double const E_CHARGE = 1.6e-19;

class Specie{

	public:

		Specie();

		Specie(std::string, int, double, int*);

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

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> n; // density (m^-3)
};