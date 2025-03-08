#include "chemistry.hpp"

// n_reagets, n_products, reagents, products, reaction delta energy
Chemistry::Chemistry(int _n_reagents, int _n_products, Specie* _reagents, Specie* _products, int _type,
 	double _react_energy_delta, double _a1, double _a2, double _a3) : n_reagents(_n_reagents), n_products(_n_products),type(_type), react_energy_delta(_react_energy_delta), a1(_a1), a2(_a2), a3(_a3){

	 // Ensure reagents and products arrays are properly allocated
    reagents = new Specie[n_reagents];
    products = new Specie[n_products];

	std::copy(_reagents, _reagents + _n_reagents, reagents);
	std::cout <<"Ad\n";
	std::copy(_products, _products + _n_products, products);
	std::cout <<"Ad\n";
} 

// Calculate the ionzitation rate or reaction rate for an Electric field in townsend
double Chemistry::calc_reaction_rate(double E_field){
	std::cout <<"hi \n";
	if (type == 0){	
		return a1 * (std::exp(-0.5*std::pow( std::log(E_field*a2) / a3 , 2)));
	} else {
		return a1;
	}
	 
}