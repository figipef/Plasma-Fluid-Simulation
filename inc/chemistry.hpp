#pragma once

#include <fstream>
#include <iostream>

class chemistry{

	public:

		chemistry(int n_species, int n_reactions, );

		double calc_reaction_rate();

		double setup_species_react_net(std::ofstream&);

	private:

		double *reaction_rates; // For each reaction (i)

		int **species_react_net; // For each species (i) save the net gain/loss for each reaction (j)

}