#pragma once

#include "Specie.hpp"

#include <fstream>
#include <iostream>
#include <cmath>

class Chemistry{

	public:

		Chemistry(int, int, Specie*, Specie*, int, double, double, double, double); // n_reagets, n_products, reagents, products, reaction delta energy

		double calc_reaction_rate(double);

		double setup_species_react_net(std::ofstream&);

		Specie* reagents;

		int n_reagents;

		Specie* products;

		int n_products;

		double react_energy_delta; // In eV !!!

		// used for the reaction rates a1  * (exp(-0.5*( log(x*a2) / a3)**2))
		double a1;
		double a2;
		double a3;

		int type; // If it is constant or not!!

};