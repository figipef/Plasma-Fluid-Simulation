#include "placeholder.hpp"
#include "poisson1DCart.hpp"
#include "poisson2DCyl.hpp"

#include <cmath>
#include <vector>

double epsi = 8.85e-12;

double no_epsi = 1;

double ZERO(double x){
    x=1;
    return 0;
}

double constant(double x){
    x=1;
    return 1 * epsi;
}

double ZERO2D(double r, double z){
    r=1;
    z = 2;
    return 0;
}

double CONST2D(double r, double z){
    r=1;
    z = 2;
    return 0.0001 * epsi;
}

double jan(double r, double z){
    return (6 - 4 * r* r - 4*z*z)*exp(-r*r -z*z);
}


int main() {

    std::map<int, double> eps_map;
    std::map<int, double> sig_map;

    double fronteira[] = {0,0,1,0};
    double fronteira_livre[] = {0,0,1,0}; // zmin, zmax, r0, rmax

    sig_map[0] = 0;

    eps_map[0] = no_epsi;
    eps_map[66] = no_epsi*2;
    eps_map[132] = no_epsi*10;

    std::vector<std::pair<int, double>> eps_vec = {{0, 1.0}, {66, 2.0}, {132, 10.0}};

    Poisson1DCart testPoisson1D(12, 1, 0.);

    testPoisson1D.dirichlet(10, false, false, 0, 10.0, &ZERO);

    //testPoisson1D.dirichlet(5, false, false, 0, 10.0, &constant);

    //testPoisson1D.dirichlet(5, true, false, 0, 10.0, &constant);

    //testPoisson1D.dirichlet(5, true, true, 0, 10, &constant);

    Poisson2DCyl testPoisson2D(200,200,0.01,0.01,eps_vec, sig_map);

    testPoisson2D.solve(10,0,0,0, &ZERO2D, fronteira_livre, "zero");


    return 0;
}