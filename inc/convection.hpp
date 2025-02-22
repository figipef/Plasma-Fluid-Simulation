#pragma once

#include <Eigen/Dense>

class Convection {

	public:
		Convection();
		Convection(double, double, Eigen::MatrixXd, Eigen::MatrixXd&);

		double calcFlux_superbee(int, int, double, double, double, int, Eigen::MatrixXd);

		double calcFlux_UNO3(int, int, double, double, double, int, Eigen::MatrixXd);

		double calcFlux_UNO2(int, int, double, double, double, int, Eigen::MatrixXd);

		double calcFlux_Koren(int, int, double, double, double, int, Eigen::MatrixXd);

	private:

		double midWayFlux(double, double, double, double, double);

		double calculate_g_c(double, double, double, double, double);

		double calculate_g_c_UNO2(double, double);

		double phia(double);

		double korenLimiter(double);

		double z_step;
		double z_size;
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S_hori;
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Ez1;

};

