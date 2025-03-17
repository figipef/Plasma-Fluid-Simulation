#include "convection.hpp"

Convection::Convection(){} 

Convection::Convection(double _z_step, double _z_size, Eigen::MatrixXd _S_hori, Eigen::MatrixXd& _Ez1){

	z_step = _z_step;
	z_size = _z_size;
	S_hori = _S_hori;
	Ez1 = _Ez1;
}

double Convection::calcFlux_superbee(int i , int j ,double mu, double De, double dt, int currlim, Eigen::MatrixXd n, int charge){
	double EPS0 = 8.85418781762e-12;
	double ECHARGE = 1.6e-19;

	double flux_left = 0.0;
	double flux_right = 0.0;

	if (j == 0) {

		double re = (n(i,j) - 0) / (n(i,j+1) - n(i,j) + 1e-8); 

        double ce = (-dt * mu * Ez1(i,j)* -charge) / z_step;

        double re_prime = (n(i,j+1) - n(i,j+2)) / (n(i,j) - n(i,j+1) + 1e-8);
        double ce_prime = -ce;

        double nem;
        if (mu * Ez1(i,j) * -charge  <= 0) {
            nem = n(i,j) + (phia(re) / 2) * (1 - ce) * (n(i,j+1) - n(i,j));
        } else {
            nem = n(i,j+1) + (phia(re_prime) / 2) * (1 - ce_prime) * (n(i,j) - n(i,j+1));
        }

        flux_left = 0.0;
        flux_right = -mu * Ez1(i,j)* -charge * nem - (De / z_step) * (n(i,j+1) - n(i,j));

        if (currlim == 1){
        	double EField_star_right = std::max(std::abs(Ez1(i,j)), (De * std::abs(n(i,j+1) - n(i,j))) / (mu*z_step * std::max({n(i,j), n(i,j+1), 1e-8})));

        	if (std::abs(flux_right) > (EPS0 * EField_star_right) / (ECHARGE * dt)) {
            	flux_right = (std::copysign(1.0, flux_right) * EPS0 * EField_star_right) / (ECHARGE * dt);
        	}
        }
        

    } else if (j == z_size - 1) {

        double rw = (n(i,j-1) - n(i,j-2)) / (n(i,j) - n(i,j-1) + 1e-8);
        double cw = (-dt * mu * Ez1(i,j-1)* -charge) / z_step;

        double rw_prime;
        if ( j == z_size - 2){
        	rw_prime = (n(i,j) - 0) / (n(i,j-1) - n(i,j) + 1e-8);
        } else {
        	rw_prime = (n(i,j) - 0) / (n(i,j-1) - n(i,j) + 1e-8);
        }

        double cw_prime = -cw;

        double nw;

        if (mu * Ez1(i,j-1)* -charge <= 0) {
            nw = n(i,j-1) + (phia(rw) / 2) * (1 - cw) * (n(i,j) - n(i,j-1));
        } else {
            nw = n(i,j) + (phia(rw_prime) / 2) * (1 - cw_prime) * (n(i,j-1) - n(i,j));
        }

        flux_left = -mu * Ez1(i,j-1)* -charge * nw - (De / (z_step) * (n(i,j) - n(i,j-1)));
        flux_right = 0;
        if (currlim == 1){
        	double EField_star_left = std::max( std::abs(Ez1(i,j-1)), (De * std::abs(n(i,j) - n(i,j-1))) / (mu*z_step * std::max({n(i,j-1), n(i,j), 1e-8})));

	        if (std::abs(flux_left) > (EPS0 * EField_star_left) / (ECHARGE * dt)) {
	            flux_left = (std::copysign(1.0, flux_left) * EPS0 * EField_star_left) / (ECHARGE * dt);
	        }
        }
        

    } else {
        // Initialize variables
    	double re, ce, rw, cw, re_prime, ce_prime, rw_prime, cw_prime;
    	double ne, nw;
    	double EField_star_right, EField_star_left;
	
    	// Calculate re, ce, rw, cw, and their derivatives
    	re = (n(i,j) - n(i,j-1)) / (n(i,j+1) - n(i,j) + 1e-8);
    	ce = (-dt * mu * Ez1(i,j)* -charge) / z_step;

		if (j == 1){
			rw = (n(i,j-1) - 0) / (n(i,j) - n(i,j-1) + 1e-8);
		} else {
			rw = (n(i,j-1) - n(i,j-2)) / (n(i,j) - n(i,j-1) + 1e-8);
		}
    	
    	cw = (-dt * mu * Ez1(i,j-1)* -charge) / z_step;

		if ( j == z_size - 2){
			re_prime = (n(i,j+1) - 0) / (n(i,j) - n(i,j+1) + 1e-8);
		} else {
			re_prime = (n(i,j+1) - n(i,j+2)) / (n(i,j) - n(i,j+1) + 1e-8);
		}
    	
    	ce_prime = -ce;
	
    	rw_prime = (n(i,j) - n(i,j+1)) / (n(i,j-1) - n(i,j) + 1e-8);
    	cw_prime = -cw;
	
    	// Compute ne
    	if (mu * Ez1(i,j)* -charge <= 0) { // Zero or positive velocity
    	    ne = n(i,j) + (phia(re) / 2.0) * (1.0 - ce) * (n(i,j+1) - n(i,j));
    	} else {
    	    ne = n(i,j+1) + (phia(re_prime) / 2.0) * (1.0 - ce_prime) * (n(i,j) - n(i,j+1));
    	}
	
    	// Compute nw
    	if (mu * Ez1(i,j-1)* -charge <= 0) {
    	    nw = n(i,j-1) + (phia(rw) / 2.0) * (1.0 - cw) * (n(i,j) - n(i,j-1));
    	} else {
    	    nw = n(i,j) + (phia(rw_prime) / 2.0) * (1.0 - cw_prime) * (n(i,j-1) - n(i,j));
    	}
	
    	// Compute fluxes
    	flux_left = -mu * Ez1(i,j-1)* -charge * nw - (De / z_step) * (n(i,j) - n(i,j-1));
	
    	flux_right = -mu * Ez1(i,j)* -charge * ne - (De / z_step) * (n(i,j+1) - n(i,j));
		
		if ( currlim == 1){
			// Compute EField_star_right and adjust flux_right if needed
    		EField_star_right = std::max(std::abs(Ez1(i,j)), (De * std::abs(n(i,j+1) - n(i,j))) / (mu*z_step * std::max({n(i,j),n(i,j+1), 1e-8})));
	
			//std::cout << flux_left<< " "<< flux_right <<std::endl;
	
    		if (std::abs(flux_right) > (EPS0 * EField_star_right) / (ECHARGE * dt)) {
    		    flux_right = std::copysign((EPS0 * EField_star_right) / (ECHARGE * dt), flux_right);
    		}
		
    		// Compute EField_star_left and adjust flux_left if needed
    		EField_star_left = std::max(std::abs(Ez1(i,j-1)),(De * std::abs(n(i,j) - n(i,j-1))) / (mu*z_step * std::max({n(i,j-1), n(i,j), 1e-8})));
		
    		if (std::abs(flux_left) > (EPS0 * EField_star_left) / (ECHARGE * dt)) {
    		    flux_left = std::copysign((EPS0 * EField_star_left) / (ECHARGE * dt), flux_left);
    		}
		}
    }
    //std::cout << flux_left<< " "<< flux_right <<std::endl;
    return (flux_left - flux_right)*S_hori(i,j);    
}

double Convection::calcFlux_UNO3(int i , int j ,double mu, double De, double dt, int currlim, Eigen::MatrixXd n, int charge, int grid_init, int grid_end){
	double EPS0 = 8.85418781762e-12;
	double ECHARGE = 1.6e-19;

	double flux_left = 0.0;
	double flux_right = 0.0;
    
	if (j == grid_init) { // Só precisamos de ne
		double g_c;
		double g_dc;
		double g_cu;
		double nem;
		double ve = -mu * Ez1(i,j)* -charge;

        if (mu * Ez1(i,j)* -charge <= 0) { // Velocidade nula ou positiva

        	g_dc = (n(i, j + 1) - n(i, j))/z_step;
			g_cu = (n(i, j) - 0)/z_step;

			g_c = calculate_g_c(g_dc, g_cu, ve, z_step, dt);
		    nem = n(i,j) + 0.5 * std::copysign(1.0, ve) * (z_step - std::abs(ve) * dt) * g_c;

        } else {

        	g_dc = -1.0 * (n(i, j) - n(i, j + 1))/z_step;
			g_cu = -1.0 *(n(i, j + 1) - n(i, j + 2))/z_step;

			g_c = calculate_g_c(g_dc, g_cu, ve, z_step, dt);
		    nem = n(i,j+1) +  0.5 * std::copysign(1.0, ve) * (z_step - std::abs(ve) * dt) * g_c;
        }

        flux_left = 0.0;
        //flux_left = -mu * Ez1(i,j) * nem - (De / z_step) * (n(i,j+1) - n(i,j)); // CHANGED
        flux_right = -mu * Ez1(i,j)* -charge * nem - (De / z_step) * (n(i,j+1) - n(i,j));

        if ( currlim == 1){
        	double EField_star_right = std::max(std::abs(Ez1(i,j)), (De * std::abs(n(i,j+1) - n(i,j))) / (mu*z_step * std::max({n(i,j), n(i,j+1), 1e-8})));

        	if (std::abs(flux_right) > (EPS0 * EField_star_right) / (ECHARGE * dt)) {
        	    flux_right = (std::copysign(1.0, flux_right) * EPS0 * EField_star_right) / (ECHARGE * dt);
        	    //flux_left = flux_right; // CHANGED
        	}
        }
        
    } else if (j == grid_end - 1) { // Só precisamos de nw
    	double g_c;
		double g_dc;
		double g_cu;
		double nw;
		double vw = -mu * Ez1(i,j-1)* -charge;
        if (mu * Ez1(i,j-1)* -charge <= 0) { // Velocidade nula ou positiva

			g_dc = (n(i, j) - n(i, j - 1))/z_step;
			g_cu = (n(i, j - 1) - n(i, j - 2))/z_step;

			g_c = calculate_g_c(g_dc, g_cu, vw, z_step, dt);
		    nw = n(i,j-1) + 0.5 * std::copysign(1.0, vw) * (z_step - std::abs(vw) * dt) * g_c;

        } else {

			g_cu = -1.0 * (n(i, j) - 0)/z_step;
			g_dc = -1.0 * (n(i, j - 1) - n(i, j))/z_step;

		    //G_A = calculate_g_c(G_WP, G_PE, vw, z_step, dt);
		    g_c = calculate_g_c(g_dc, g_cu, vw, z_step, dt);
		    nw = n(i,j) + 0.5 * std::copysign(1.0, vw) * (z_step - std::abs(vw) * dt) * g_c;
        }

        flux_left = -mu * Ez1(i,j-1)* -charge * nw - (De / z_step) * (n(i,j) - n(i,j-1));
    	flux_right = 0;
    	//flux_right = flux_left; // CHANGED
    	if ( currlim == 1){
    		double EField_star_left = std::max( std::abs(Ez1(i,j-1)), (De * std::abs(n(i,j) - n(i,j-1))) / (mu*z_step * std::max({n(i,j-1), n(i,j), 1e-8})));

        	if (std::abs(flux_left) > (EPS0 * EField_star_left) / (ECHARGE * dt)) {
        	    flux_left = (std::copysign(1.0, flux_left) * EPS0 * EField_star_left) / (ECHARGE * dt);
        	    //flux_right = flux_left; // CHANGED
        	}
    	}
	
    } else { // Precisamos de ne e de nw

		double ve = -mu * Ez1(i,j)* -charge;
		double vw = -mu * Ez1(i,j-1)* -charge;

		double g_c = 0;
		double g_cu = 0;
		double g_dc = 0;
		double nem = 0.0;
		double nw = 0.0;
		
		// Handle velocity for ne
		if (mu * Ez1(i,j)* -charge <= 0) { // Zero or positive velocity

			g_dc = (n(i, j + 1) - n(i, j))/z_step;
			g_cu = (n(i, j) - n(i, j - 1))/z_step;

			g_c = calculate_g_c(g_dc, g_cu, ve, z_step, dt);
		    nem = n(i,j) + 0.5 * std::copysign(1.0, ve) * (z_step - std::abs(ve) * dt) * g_c;



		} else {

			g_dc = -1.0 * (n(i, j) - n(i, j + 1))/z_step;

			if ( j == z_size - 2){
				g_cu = -1.0 *(n(i, j + 1) -0)/z_step;
			}else{
				g_cu = -1.0 *(n(i, j + 1) - n(i, j + 2))/z_step;
			}
	    	
		    g_c = calculate_g_c(g_dc, g_cu, ve, z_step, dt);

		    nem = n(i,j+1) +  0.5 * std::copysign(1.0, ve) * (z_step - std::abs(ve) * dt) * g_c;
		   	
		}
		
		// Handle velocity for nw
		if (mu * Ez1(i,j-1)* -charge <= 0) { // Zero or positive velocity
			
			g_dc = (n(i, j) - n(i, j - 1))/z_step;

			if ( j == 1){
				g_cu = (n(i, j - 1) - 0)/z_step;
			}else{
				g_cu = (n(i, j - 1) - n(i, j - 2))/z_step;
			}

		    g_c = calculate_g_c(g_dc, g_cu, vw, z_step, dt);
		    nw = n(i,j-1) + 0.5 * std::copysign(1.0, vw) * (z_step - std::abs(vw) * dt) * g_c;



		} else {
			g_cu = -1.0 * (n(i, j) - n(i, j + 1))/z_step;
			g_dc = -1.0 * (n(i, j - 1) - n(i, j))/z_step;

		    g_c = calculate_g_c(g_dc, g_cu, vw, z_step, dt);
		    nw = n(i,j) + 0.5 * std::copysign(1.0, vw) * (z_step - std::abs(vw) * dt) * g_c;
		}

		// Compute fluxes
    	flux_left = vw * nw - (De / z_step) * (n(i,j) - n(i,j-1));
	
    	flux_right = ve * nem - (De / z_step) * (n(i,j+1) - n(i,j));
		
		if ( currlim == 1){
			// Compute EField_star_right and adjust flux_right if needed
    		double EField_star_right = std::max(std::abs(Ez1(i,j)), (De * std::abs(n(i,j+1) - n(i,j))) / (mu*z_step * std::max({n(i,j),n(i,j+1), 1e-8})));
	
			//std::cout << flux_left<< " "<< flux_right <<std::endl;
	
    		if (std::abs(flux_right) > (EPS0 * EField_star_right) / (ECHARGE * dt)) {
    		    flux_right = std::copysign((EPS0 * EField_star_right) / (ECHARGE * dt), flux_right);
    		}
		
    		// Compute EField_star_left and adjust flux_left if needed
    		double EField_star_left = std::max(std::abs(Ez1(i,j-1)),(De * std::abs(n(i,j) - n(i,j-1))) / (mu*z_step * std::max({n(i,j-1), n(i,j), 1e-8})));
		
    		if (std::abs(flux_left) > (EPS0 * EField_star_left) / (ECHARGE * dt)) {
    		    flux_left = std::copysign((EPS0 * EField_star_left) / (ECHARGE * dt), flux_left);
    		}
		}
    }

    return (flux_left - flux_right)*S_hori(i,j);
}

double Convection::calcFlux_Koren(int i , int j ,double mu, double De, double dt, int currlim, Eigen::MatrixXd n, int charge){
	double EPS0 = 8.85418781762e-12;
	double ECHARGE = 1.6e-19;

	double flux_left = 0.0;
	double flux_right = 0.0;

	if (j == 0) { // Só precisamos de ne

		double ve = -mu * Ez1(i,j)* -charge;

        if (ve >= 0) { // Velocidade nula ou positiva

        	double r_i = (n(i,j) - 0.0 + 1e-308)/(n(i,j+1) - n(i,j) + 1e-308);
        	flux_right = ve * (n(i,j) + korenLimiter(r_i)*(n(i,j+1) - n(i,j)));

        } else {

        	double r_i = (n(i,j+1) - n(i, j) + 1e-308)/(n(i,j+2) - n(i,j+1) + 1e-308);
        	flux_right = ve * (n(i,j+1) - korenLimiter(1/r_i)*(n(i,j+1) - n(i,j)));

        }

        flux_left = 0.0;
        flux_right = flux_right - (De / z_step) * (n(i,j+1) - n(i,j));

        if ( currlim == 1){
        	double EField_star_right = std::max(std::abs(Ez1(i,j)), (De * std::abs(n(i,j+1) - n(i,j))) / (mu*z_step * std::max({n(i,j), n(i,j+1), 1e-8})));

        	if (std::abs(flux_right) > (EPS0 * EField_star_right) / (ECHARGE * dt)) {
        	    flux_right = (std::copysign(1.0, flux_right) * EPS0 * EField_star_right) / (ECHARGE * dt);
        	}
        }
        
    } else if (j == z_size - 1) { // Só precisamos de nw

		double vw = -mu * Ez1(i,j-1)* -charge;

        if (vw >= 0) { // Velocidade nula ou positiva

        	double r_i = (n(i,j - 1) - n(i, j - 2) + 1e-308)/(n(i,j) - n(i,j - 1) + 1e-308);
        	flux_left = vw * (n(i,j - 1) + korenLimiter(r_i)*(n(i,j) - n(i,j - 1)));

        } else {

        	double r_i = (n(i,j) - n(i, j - 1) + 1e-308)/(0.0 - n(i,j) + 1e-308);
        	flux_left = vw * (n(i,j) - korenLimiter(1/r_i)*(n(i,j) - n(i,j - 1)));

        }

        flux_left = flux_left - (De / z_step) * (n(i,j) - n(i,j-1));
    	flux_right = 0;

    	if ( currlim == 1){
    		double EField_star_left = std::max( std::abs(Ez1(i,j-1)), (De * std::abs(n(i,j) - n(i,j-1))) / (mu*z_step * std::max({n(i,j-1), n(i,j), 1e-8})));

        	if (std::abs(flux_left) > (EPS0 * EField_star_left) / (ECHARGE * dt)) {
        	    flux_left = (std::copysign(1.0, flux_left) * EPS0 * EField_star_left) / (ECHARGE * dt);
        	}
    	}
	
    } else { // Precisamos de ne e de nw

		double ve = -mu * Ez1(i,j)* -charge;
		double vw = -mu * Ez1(i,j-1)* -charge;
		
		// Handle velocity for ne
		if (ve >= 0) { // Zero or positive velocity

			double r_i = (n(i,j) - n(i, j -1) + 1e-308)/(n(i,j+1) - n(i,j) + 1e-308);
        	flux_right = ve * (n(i,j) + korenLimiter(r_i)*(n(i,j+1) - n(i,j)));

		} else {

			double r_i;
        	
			if ( j == z_size - 2){
				r_i = (n(i,j+1) - n(i, j)+ 1e-308)/(0.0 - n(i,j+1) + 1e-308);
			}else{
				r_i = (n(i,j+1) - n(i, j)+ 1e-308)/(n(i,j+2) - n(i,j+1) + 1e-308);
			}

			flux_right = ve * (n(i,j+1) - korenLimiter(1/r_i)*(n(i,j+1) - n(i,j)));
		}
		
		// Handle velocity for nw
		if (vw >= 0) { // Zero or positive velocity
			
			double r_i;

			if ( j == 1){
				r_i = (n(i,j - 1) - 0.0 + 1e-308)/(n(i,j) - n(i,j - 1) + 1e-308);
			}else{
				r_i = (n(i,j - 1) - n(i, j - 2) + 1e-308)/(n(i,j) - n(i,j - 1) + 1e-308);
			}

			
        	flux_left = vw * (n(i,j - 1) + korenLimiter(r_i)*(n(i,j) - n(i,j - 1)));

		} else {

			double r_i = (n(i,j) - n(i, j - 1) + 1e-308)/(n(i,j + 1) - n(i,j) + 1e-308);
        	flux_left = vw * (n(i,j) - korenLimiter(1/r_i)*(n(i,j) - n(i,j - 1)));

		}

		// Compute fluxes
    	flux_left = flux_left - (De / z_step) * (n(i,j) - n(i,j-1));
	
    	flux_right = flux_right - (De / z_step) * (n(i,j+1) - n(i,j));
		
		if ( currlim == 1){
			// Compute EField_star_right and adjust flux_right if needed
    		double EField_star_right = std::max(std::abs(Ez1(i,j)), (De * std::abs(n(i,j+1) - n(i,j))) / (mu*z_step * std::max({n(i,j),n(i,j+1), 1e-8})));
	
			//std::cout << flux_left<< " "<< flux_right <<std::endl;
	
    		if (std::abs(flux_right) > (EPS0 * EField_star_right) / (ECHARGE * dt)) {
    		    flux_right = std::copysign((EPS0 * EField_star_right) / (ECHARGE * dt), flux_right);
    		}
		
    		// Compute EField_star_left and adjust flux_left if needed
    		double EField_star_left = std::max(std::abs(Ez1(i,j-1)),(De * std::abs(n(i,j) - n(i,j-1))) / (mu*z_step * std::max({n(i,j-1), n(i,j), 1e-8})));
		
    		if (std::abs(flux_left) > (EPS0 * EField_star_left) / (ECHARGE * dt)) {
    		    flux_left = std::copysign((EPS0 * EField_star_left) / (ECHARGE * dt), flux_left);
    		}
		}
    }

    return (flux_left - flux_right)*S_hori(i,j);
}

double Convection::calcFlux_UNO2(int i , int j ,double mu, double De, double dt, int currlim, Eigen::MatrixXd n, int charge){
	double EPS0 = 8.85418781762e-12;
	double ECHARGE = 1.6e-19;

	double flux_left = 0.0;
	double flux_right = 0.0;
    
	if (j == 0) { // Só precisamos de ne
		double g_c;
		double g_dc;
		double g_cu;
		double nem;
		double ve = -mu * Ez1(i,j)* -charge;

        if (mu * Ez1(i,j)* -charge <= 0) { // Velocidade nula ou positiva

        	g_dc = (n(i, j + 1) - n(i, j))/z_step;
			g_cu = (n(i, j) - 0)/z_step;

			g_c = calculate_g_c_UNO2(g_dc, g_cu);
		    nem = n(i,j) + 0.5 * std::copysign(1.0, ve) * (z_step - std::abs(ve) * dt) * g_c;

        } else {

        	g_dc = -1.0 * (n(i, j) - n(i, j + 1))/z_step;
			g_cu = -1.0 *(n(i, j + 1) - n(i, j + 2))/z_step;

			g_c = calculate_g_c_UNO2(g_dc, g_cu);
		    nem = n(i,j+1) +  0.5 * std::copysign(1.0, ve) * (z_step - std::abs(ve) * dt) * g_c;
        }

        flux_left = 0.0;
        //flux_left = -mu * Ez1(i,j) * nem - (De / z_step) * (n(i,j+1) - n(i,j)); // CHANGED
        flux_right = -mu * Ez1(i,j)* -charge * nem - (De / z_step) * (n(i,j+1) - n(i,j));

        if ( currlim == 1){
        	double EField_star_right = std::max(std::abs(Ez1(i,j)), (De * std::abs(n(i,j+1) - n(i,j))) / (mu*z_step * std::max({n(i,j), n(i,j+1), 1e-8})));

        	if (std::abs(flux_right) > (EPS0 * EField_star_right) / (ECHARGE * dt)) {
        	    flux_right = (std::copysign(1.0, flux_right) * EPS0 * EField_star_right) / (ECHARGE * dt);
        	    //flux_left = flux_right; // CHANGED
        	}
        }
        
    } else if (j == z_size - 1) { // Só precisamos de nw
    	double g_c;
		double g_dc;
		double g_cu;
		double nw;
		double vw = -mu * Ez1(i,j-1)* -charge;
        if (mu * Ez1(i,j-1)* -charge <= 0) { // Velocidade nula ou positiva

			g_dc = (n(i, j) - n(i, j - 1))/z_step;
			g_cu = (n(i, j - 1) - n(i, j - 2))/z_step;

			g_c = calculate_g_c_UNO2(g_dc, g_cu);
		    nw = n(i,j-1) + 0.5 * std::copysign(1.0, vw) * (z_step - std::abs(vw) * dt) * g_c;

        } else {

			g_cu = -1.0 * (n(i, j) - 0)/z_step;
			g_dc = -1.0 * (n(i, j - 1) - n(i, j))/z_step;

		    //G_A = calculate_g_c(G_WP, G_PE, vw, z_step, dt);
		    g_c = calculate_g_c_UNO2(g_dc, g_cu);
		    nw = n(i,j) + 0.5 * std::copysign(1.0, vw) * (z_step - std::abs(vw) * dt) * g_c;
        }

        flux_left = -mu * Ez1(i,j-1)* -charge * nw - (De / z_step) * (n(i,j) - n(i,j-1));
    	flux_right = 0;
    	//flux_right = flux_left; // CHANGED
    	if ( currlim == 1){
    		double EField_star_left = std::max( std::abs(Ez1(i,j-1)), (De * std::abs(n(i,j) - n(i,j-1))) / (mu*z_step * std::max({n(i,j-1), n(i,j), 1e-8})));

        	if (std::abs(flux_left) > (EPS0 * EField_star_left) / (ECHARGE * dt)) {
        	    flux_left = (std::copysign(1.0, flux_left) * EPS0 * EField_star_left) / (ECHARGE * dt);
        	    //flux_right = flux_left; // CHANGED
        	}
    	}
	
    } else { // Precisamos de ne e de nw

		double ve = -mu * Ez1(i,j)* -charge;
		double vw = -mu * Ez1(i,j-1)* -charge;

		double g_c = 0;
		double g_cu = 0;
		double g_dc = 0;
		double nem = 0.0;
		double nw = 0.0;
		
		// Handle velocity for ne
		if (mu * Ez1(i,j)* -charge <= 0) { // Zero or positive velocity

			g_dc = (n(i, j + 1) - n(i, j))/z_step;
			g_cu = (n(i, j) - n(i, j - 1))/z_step;

			g_c = calculate_g_c_UNO2(g_dc, g_cu);
		    nem = n(i,j) + 0.5 * std::copysign(1.0, ve) * (z_step - std::abs(ve) * dt) * g_c;



		} else {

			g_dc = -1.0 * (n(i, j) - n(i, j + 1))/z_step;

			if ( j == z_size - 2){
				g_cu = -1.0 *(n(i, j + 1) -0)/z_step;
			}else{
				g_cu = -1.0 *(n(i, j + 1) - n(i, j + 2))/z_step;
			}
	    	
		    g_c = calculate_g_c_UNO2(g_dc, g_cu);

		    nem = n(i,j+1) +  0.5 * std::copysign(1.0, ve) * (z_step - std::abs(ve) * dt) * g_c;
		   	
		}
		
		// Handle velocity for nw
		if (mu * Ez1(i,j-1)* -charge <= 0) { // Zero or positive velocity
			
			g_dc = (n(i, j) - n(i, j - 1))/z_step;

			if ( j == 1){
				g_cu = (n(i, j - 1) - 0)/z_step;
			}else{
				g_cu = (n(i, j - 1) - n(i, j - 2))/z_step;
			}

		    g_c = calculate_g_c_UNO2(g_dc, g_cu);
		    nw = n(i,j-1) + 0.5 * std::copysign(1.0, vw) * (z_step - std::abs(vw) * dt) * g_c;



		} else {
			g_cu = -1.0 * (n(i, j) - n(i, j + 1))/z_step;
			g_dc = -1.0 * (n(i, j - 1) - n(i, j))/z_step;

		    g_c = calculate_g_c_UNO2(g_dc, g_cu);
		    nw = n(i,j) + 0.5 * std::copysign(1.0, vw) * (z_step - std::abs(vw) * dt) * g_c;
		}

		// Compute fluxes
    	flux_left = vw * nw - (De / z_step) * (n(i,j) - n(i,j-1));
	
    	flux_right = ve * nem - (De / z_step) * (n(i,j+1) - n(i,j));
		
		if ( currlim == 1){
			// Compute EField_star_right and adjust flux_right if needed
    		double EField_star_right = std::max(std::abs(Ez1(i,j)), (De * std::abs(n(i,j+1) - n(i,j))) / (mu*z_step * std::max({n(i,j),n(i,j+1), 1e-8})));
	
			//std::cout << flux_left<< " "<< flux_right <<std::endl;
	
    		if (std::abs(flux_right) > (EPS0 * EField_star_right) / (ECHARGE * dt)) {
    		    flux_right = std::copysign((EPS0 * EField_star_right) / (ECHARGE * dt), flux_right);
    		}
		
    		// Compute EField_star_left and adjust flux_left if needed
    		double EField_star_left = std::max(std::abs(Ez1(i,j-1)),(De * std::abs(n(i,j) - n(i,j-1))) / (mu*z_step * std::max({n(i,j-1), n(i,j), 1e-8})));
		
    		if (std::abs(flux_left) > (EPS0 * EField_star_left) / (ECHARGE * dt)) {
    		    flux_left = std::copysign((EPS0 * EField_star_left) / (ECHARGE * dt), flux_left);
    		}
		}
    }

    return (flux_left - flux_right)*S_hori(i,j);
}

double Convection::midWayFlux(double donor_cell, double u, double dx, double dt, double g_c){
	return donor_cell + 0.5 * std::copysign(1.0,u)*(dx - std::abs(u) * dt) * g_c;
}

double Convection::calculate_g_c(double g_dc, double g_cu, double u, double dx, double dt){
	
	if (std::abs(g_dc - g_cu) < 0.6*std::abs(g_dc + g_cu)){
		//std::cout <<"a"<<std::endl;
		//double xmf = std::copysign(1.0, u) * (dx - std::abs(u) * dt)/2.0;
		return g_dc - (dx*copysign(1.0,-u) + std::abs(u)*dt)/(1.5*std::copysign(1.0,u))*(g_dc - g_cu)/(2.0 *dx*copysign(1.0,u));
		//return g_dc - 4/3*(std::copysign(1.0, u)*dx - xmf)*(g_dc - g_cu)/(2.0 *dx* std::copysign(1.0, u));

	} else if (g_dc * g_cu > 0.0 ){
		//std::cout <<"b"<<std::endl;
		return std::copysign(1.0, g_dc)*2.0*std::min(std::abs(g_dc),std::abs(g_cu));
	
	}else{
		//std::cout <<"c"<<std::endl;
		return std::copysign(1.0, g_dc)*std::min(std::abs(g_dc),std::abs(g_cu));
	}
	
	//return std::copysign(1.0, g_dc)*2*(std::abs(g_dc * g_cu)/(std::abs(g_dc) + std::abs(g_cu)+ 1e-308));
}

double Convection::calculate_g_c_UNO2(double g_dc, double g_cu){

	return std::copysign(1.0, g_dc)*2*(std::abs(g_dc * g_cu)/(std::abs(g_dc) + std::abs(g_cu)+ 1e-308));
}

double Convection::phia(double x){
	//return std::max(0., std::min(std::min(1., (2+x)/6 ), x));
	return std::max({0., std::max(std::min(1., x*2), std::min(x, 2.))});
}

double Convection::korenLimiter(double x){
	return std::max(0., std::min({1.0, (2.0+x)/6.0 , x}));
}

