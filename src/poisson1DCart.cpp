#include "poisson1DCart.hpp"

#include <Eigen/Dense>
#include <iostream>
#include <vector>

Poisson1DCart::Poisson1DCart(int n, double dx, double xmin) {

	double epsi0 = 8.85e-12;

	size = n;
	stepx = dx;

	xgrid = new int[size];  // Dynamically allocate an array of size n
	v = new double[size];

	phic = new double[size];
	phiw = new double[size];
	phie = new double[size];

	for (int i = 0; i < size; i++) {

		xgrid[i] = xmin + i * dx;  // Initialize array elements (optional)
		v[i] = 0;

		phie[i] = -epsi0/dx;
		phiw[i] = -epsi0/dx;
		phic[i] = -phie[i] - phiw[i];
	}

    std::cout<< "grid" <<std::endl;

    for (int i = 0; i < size; i++) {
    	std::cout<< xgrid[i]<< "; ";
    }
}

void Poisson1DCart::placeholderfunc(){

	std::cout<< "Testing" <<std::endl;

}

void Poisson1DCart::dirichlet(int n, bool nl, bool nr, double VMIN, double VAPP, double (*func)(double)){

	double epsi0 = 8.85e-12;

	Eigen::MatrixXd m = Eigen::MatrixXd::Zero(n,n);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(n);

	if (nl == 1){
		phic[0] = -phiw[0];	
	} else {
		phic[0] = -phie[0] - phiw[0];
	}

	if (nr == 1){
		phic[n-1] = -phie[n];
	} else{
		phic[n-1] = -phie[n] - phiw[n];
	}

	for (int i = 0; i <= n-1; i++){
		
		if (i >= 1){
			m(i,i-1) = phiw[i];
		}



		m(i,i) = phic[i];

		if (i < n-1){
			m(i,i+1) = phie[i];
		}

		b(i) = func(xgrid[i]) * stepx;

	}

	if ( nl == false){
		b(0) = -1 * VMIN * phiw[0] + b(0);
	}
	else {
		b(0) = -epsi0 * VMIN + b(0);
	}

	if (nr == false){
		b(n-1) = -1 * VAPP * phie[n-1] + b(n-1);
	}else{
		b(n-1) = -epsi0 * VAPP + b(0);
	}
	

	std::cout << "Here is the matrix A:\n" << m << std::endl;
	std::cout << "Here is the vector b:\n" << b << std::endl;

	Eigen::VectorXd x = m.colPivHouseholderQr().solve(b);

	std::cout << "The solution is:\n" << x << std::endl;

}

void Poisson1DCart::print_voltages(){
	std::cout<<"Voltage Array"<<std::endl;

	for (int i = 0; i < size; i++) {

		std::cout<<v[i]<<", ";
	}
	std::cout<<std::endl;
}

// Destructor
Poisson1DCart::~Poisson1DCart() {
	delete[] xgrid;  // Deallocate the array memory
}