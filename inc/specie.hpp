#pragma once

#include <string>
#include <Eigen/Dense>

double const E_CHARGE = 1.6e-19;

class Specie{

	public:

		Specie();

		Specie(std::string, int, double, int*, Eigen::MatrixXd&);

		std::string get_name();

		double get_mass();

		int get_charge();

		double get_qm_ratio();

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> get_density();

		void update_density(Eigen::MatrixXd&);

	private:

		std::string name;

		int charge; // whole units

		double mass;

		int *react_net; // net value for each reaction in the array form

		double qm_ratio; // charge mass ratio may be needed for some calculations

		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> n; // density (m^-3)
};