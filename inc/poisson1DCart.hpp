#pragma once

class Poisson1DCart{

	public:
		double stepx;
		int size;
		
		int *xgrid;
		
		double *v;
		
		double *phie;
		double *phiw;
		double *phic;

		Poisson1DCart(int, double, double);

		void placeholderfunc();
		void dirichlet(int, bool, bool, double, double,  double (*)(double));
		void print_voltages();
		~Poisson1DCart();

};
