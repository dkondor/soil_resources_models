/*
 * pop_model1.cpp -- numerically solve analytic population models with
 * 	two components
 * 
 * Note: the GSL library is used for calculating numerical solutions;
 * see general info here:
 * https://www.gnu.org/software/gsl/
 * and the documentation on ODE solvers here:
 * https://www.gnu.org/software/gsl/doc/html/ode-initval.html
 * 
 * Typical compilation (on Linux, with GSL installed system-wide):
g++ -o pm1 pop_model1.cpp -lgsl -lgslcblas -O3 -march=native -lm
 * 
 * Copyright 2022 Daniel Kondor <kondor@csh.ac.at>
 * 
 */

#ifndef HAVE_INLINE
#define HAVE_INLINE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <random>
#include <time.h>
#include <stdexcept>

struct params {
	enum class type {
		NO_TYPE,
		LOGISTIC, /* logistic growth with no resources */
		LV, /* Lotka-Volterra model */
		TURCHIN, /* Turchin-Korotayev model of population and warfare (Model IV) */
		BT, /* Volterra model (Model I) */
		BT_LIN, /* Volterra model with linear resource renewal (Model II) */
		SOIL1, /* simple version of the soil model with resource depletion independent of resource levels (not used) */
		SOIL2, /* soil model (Model III in the SI) with resource depletion proportional to both resources and population */
		SOIL_M, /* soil model with manuring / other input proportional to population */
	};
	type t = type::NO_TYPE;
	double r = 0.01;
	double K = 1000.0;
	double a = 0.1;
	double b = 0.01;
	double c = 0.01;
	double d = 0.01;
	double e = 0.01;
	
	static bool try_parse_type(const char* x, type& t) {
		if(!x) return false;
		
		if(x[0] == 'l' || x[0] == 'L') {
			if(x[1] == 'o' || x[1] == 'O') {
				if(x[2] == 'g' || x[2] == 'G') {
					t = type::LOGISTIC;
					return true;
				}
				if(x[2] == 't' || x[2] == 'T') {
					t = type::LV;
					return true;
				}
			}
			if(x[1] == 'v' || x[1] == 'V') {
				t = type::LV;
				return true;
			}
			return false;
		}
		
		if(x[0] == 't' || x[1] == 'T') {
			t = type::TURCHIN;
			return true;
		}
		
		if(x[0] == 'b' || x[0] == 'B') {
			if(x[1] == 't' || x[1] == 'T') switch(x[2]) {
				case 'e':
				case 'E':
					t = type::BT_LIN;
					return true;
				case 0:
					t = type::BT;
					return true;
				default:
					return false;
			}
			return false;
		}
		
		if(x[0] == 's' || x[0] == 'S') switch(x[1]) {
			case '1':
				t = type::SOIL1;
				return true;
			case '2':
				t = type::SOIL2;
				return true;
			case 'm':
				t = type::SOIL_M;
				return true;
			default:
				return false;
		}
		
		return false;
	}
	
	// return if x should be limited to strictly below K
	static bool limit_x(type t) {
		switch(t) {
			case type::TURCHIN:
			case type::SOIL1:
			case type::SOIL2:
			case type::SOIL_M:
				return true;
			default:
				return false;
		}
	}
	bool limit_x() const { return limit_x(t); }
	
	// return if y should be limited to strictly below K
	static bool limit_y(type t) {
		switch(t) {
			case type::BT:
			case type::BT_LIN:
			case type::SOIL1:
			case type::SOIL2:
			case type::SOIL_M:
				return true;
			default:
				return false;
		}
	}
	bool limit_y() const { return limit_y(t); }
	
	// get the value of the stationary solution
	void get_sol(double x[2]) const {
		switch(t) {
			case type::LOGISTIC:
				x[0] = K;
				x[1] = 0.0;
				break;
			case type::BT:
				x[0] = c * K * ( 1.0 - d / b) / a;
				x[1] = d * K / b;
				break;
			case type::BT_LIN:
				x[0] = c * K * ( b / d - 1.0) / a;
				x[1] = d * K / b;
				break;
			case type::SOIL2:
				x[0] = K * c * (sqrt(1.0 + 4.0*a/c) - 1.0) / (2.0 * a);
				x[1] = x[0];
				break;
			case type::SOIL_M:
				x[0] = K * ( e - c + sqrt( (e + c)*(e + c) + 4.0 * a * c) ) / (2.0 * (a + e));
				x[1] = x[0];
				break;
			case type::TURCHIN:
				x[0] = 0.5 * b * r * (sqrt(1.0 + 4.0*a / (b*r)) - 1) / a;
				x[1] = r * (1.0 - x[0]);
				x[0] *= K;
				break;
			default:
				throw std::runtime_error("Stationary solution not implemented for the requested model type!\n");
		}
	}
	
	/*
	 * Calculate the eigenvalues of the Jacobian at the stationary solution,
	 * return the two parts of the expression and true / false based on
	 * whether there is an imaginary part.
	 * I.e. if the return value is true, than the eigenvalues are
	 * x[0] +/- i x[1]
	 * otherwise
	 * x[0] +/- x[1]
	 */
	bool get_jac_comp(double x[2]) const {
		double tmp1;
		double scale;
		switch(t) {
			case type::BT:
				tmp1 = c * d / b;
				x[0] = -0.5 * tmp1;
				tmp1 = tmp1 * tmp1 - 4 * tmp1 * (b - d);
				scale = 0.5;
				break;
			case type::SOIL2:
				{
					double xs = c * (sqrt(1.0 + 4.0*a/c) - 1.0) / (2.0 * a);
					tmp1 = -r - c - a * xs;
					x[0] = tmp1 / 2.0;
					tmp1 = tmp1 * tmp1 - 4.0 * c * r - 8.0 * a * r * xs;
					scale = 0.5;
				}
				break;
			case type::BT_LIN:
				tmp1 = b * c / d;
				x[0] = -1.0 * tmp1 / 2.0;
				tmp1 = tmp1 * tmp1 - 4.0 * c * (b - d);
				scale = 0.5;
				break;
			case type::TURCHIN:
				{
					double xs = 0.5 * b * r * (sqrt(1.0 + 4.0*a / (b*r)) - 1) / a;
					tmp1 = -b - r * xs;
					x[0] = 0.5 * tmp1;
					tmp1 = tmp1 * tmp1 - 8.0 * b * r + 4.0 * b * r * xs;
					scale = 0.5;
					break;
				}
			default:
				throw std::runtime_error("Jacobian eigenvalues not implemented for the requested model type!\n");
		}
		
		if(tmp1 >= 0.0) {
			x[1] = sqrt(tmp1) * scale;
			return false;
		}
		x[1] = sqrt(-1.0 * tmp1) * scale;
		return true;
	}
	
	/* print (the beginning) of a file header with the names of important parameters */
	void print_header(FILE* f, char sep, bool newline = false)  const {
		switch(t) {
			case type::BT:
			case type::BT_LIN:
				fprintf(f, "a%1$cb%1$cc%1$cd%1$cr%1$cK", sep);
				break;
			case type::SOIL2:
				fprintf(f, "a%1$cc%1$cr%1$cK", sep);
				break;
			case type::SOIL_M:
				fprintf(f, "a%1$cc%1$cr%1$ce%1$cK", sep);
				break;
			default:
				throw std::runtime_error("Parameter output not implemented for the requested model type!\n");
		}
		if(newline) putc('\n', f);
		else putc(sep, f);
	}
	
	/* print the current value of parameters that are used by the current model */
	void print_pars(FILE* f, char sep, bool newline = false)  const {
		switch(t) {
			case type::BT:
			case type::BT_LIN:
				fprintf(f,"%2$f%1$c%3$f%1$c%4$f%1$c%5$f%1$c%6$f%1$c%7$f", sep, a, b, c, d, r, K);
				break;
			case type::SOIL2:
				fprintf(f,"%2$f%1$c%3$f%1$c%4$f%1$c%5$f", sep, a, c, r, K);
				break;
			case type::SOIL_M:
				fprintf(f,"%2$f%1$c%3$f%1$c%4$f%1$c%5$f%1$c%6$f", sep, a, c, r, e, K);
				break;
			default:
				throw std::runtime_error("Parameter output not implemented for the requested model type!\n");
		}
		if(newline) putc('\n', f);
		else putc(sep, f);
	}
	
	
	/* set the a parameter automatically based on the values of b, c, d
	 * in the case of the Volterra model and its linear variant */
	void set_auto_a() {
		switch(t) {
			case type::BT:
				a = (1.0 - d / b) * c * b / d;
				r = b - d;
				break;
			case type::BT_LIN:
				a = (1.0 - d / b) * c * b * b / d / d;
				r = b - d;
				break;
			default:
				throw std::runtime_error("Auto-adjustment of the a parameter is not supported for this model variant!\n");
		}
	}
	
	/* set the b and d parameters automatically based on the values of a, c, r
	 * in the case of the Volterra model and its linear variant */
	void set_auto_bd() {
		switch(t) {
			case type::BT:
				d = r * c / a;
				b = d + r;
				break;
			case type::BT_LIN:
				d = 0.5 * c * r * (1.0 + sqrt(1 + 4.0*a/c)) / a;
				b = a*d*d/c/r;
				break;
			default:
				throw std::runtime_error("Auto-adjustment of the b and d parameters is not supported for this model variant!\n");
		}
	}
};


/* calculate the derivatives */
int dfunc2(double t, const double y1[], double dydt[], const params& pars) {
	double x = y1[0];
	double y = y1[1];
	
	switch(pars.t) {
		case params::type::LOGISTIC:
			dydt[0] = pars.r * x * (1 - x / pars.K);
			dydt[1] = 0.0;
			break;
		case params::type::LV:
			dydt[0] = pars.r * x - pars.d * x * y;
			dydt[1] = pars.a * x * y - pars.b * y;
			break;
		case params::type::TURCHIN:
			dydt[0] = pars.r * x * (1 - x / pars.K) - pars.c * x * y; // note: could be c * x * y
			dydt[1] = pars.a * x * x / (pars.K * pars.K) - pars.b * y;
			break;
		case params::type::BT:
			dydt[0] = -1.0 * pars.d * x + pars.b * x * y / pars.K;
			dydt[1] = pars.c * y * (1 - y / pars.K) - pars.a * x * y / pars.K;
			break;
		case params::type::BT_LIN:
			dydt[0] = -1.0 * pars.d * x + pars.b * x * y / pars.K;
			dydt[1] = pars.c * (pars.K - y) - pars.a * x * y / pars.K;
			break;
		case params::type::SOIL1:
			dydt[0] = pars.r * x * (1 - x / y);
			dydt[1] = pars.c * (pars.K - y) - pars.a * x;
			break;
		case params::type::SOIL2:
			dydt[0] = pars.r * x * (1 - x / y);
			dydt[1] = pars.c * (pars.K - y) - pars.a * x * y / pars.K;
			break;
		case params::type::SOIL_M:
			dydt[0] = pars.r * x * (1 - x / y);
			dydt[1] = (pars.c + pars.e * x / pars.K) * (pars.K - y) - pars.a * x * y / pars.K;
			break;
		default:
			return GSL_EBADFUNC;
	}
	
	return GSL_SUCCESS;
}

int dfunc(double t, const double y1[], double dydt[], void* params1) {
	const params& pars = *(const params*)params1;
	return dfunc2(t, y1, dydt, pars);
}


/* get which "segment" of the phase space we are based on the sign of
 * the derivatives */
int get_segment(double dydt[2]) {
	const double eps = 1e-15;
	if(fabs(dydt[0]) < eps || fabs(dydt[1]) < eps) return 0;
	if(dydt[0] > 0.0) {
		if(dydt[1] > 0.0) return 1;
		else return 2;
	}
	else if(dydt[1] > 0.0) return 4;
	return 3;
}


/* run one instance of the simulation from the given initial conditions using the given parameters */
int do_one_run(const params& pars, const double y0[2], double tmax, double dt, bool fixed_step, double noise_x,
		double noise_y, bool discrete_noise, FILE* fout, FILE* cout, bool cycles_summary, int target_segment,
		std::mt19937_64& rng, char sep, bool output_params) {
	
	std::normal_distribution<double> nx(0.0, noise_x > 0.0 ? noise_x : 1.0);
	std::normal_distribution<double> ny(0.0, noise_y > 0.0 ? noise_y : 1.0);
	
	gsl_odeiv2_system sys = {dfunc, nullptr, 2, const_cast<params*>(&pars)};
	gsl_odeiv2_driver* d = nullptr;
	if(fixed_step <= 0.0) d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
	
	double t = 0.0;
	double y1[] = {y0[0], y0[1]};
	double dydt[] = {0.0, 0.0};
	int segment = 0; // which segment of the phase space we are in
	double last_t = -1.0; // last time when segment == target_segment
	double ys[] = {0.0, 0.0}; // stationary solution
	double dst[] = {0.0, 0.0}; // distance from the stationary solution
	double step1 = -1.0; // decrease in distance in the first cycle
	double overshoot = -1.0; // maximum overshoot (only calculated if cycles_summary == true)
	const double dist_eps = 1e-6; // stop running in cycles_summary mode if we are this close to the stationary solution
	
	
	if(fout) fprintf(fout, "%f\t%f\t%f\n", t, y1[0], y1[1]);
	
	// keep track of the derivatives in dydt
	dfunc2(t, y1, dydt, pars);
	
	// calculate the eigenvalues of the Jacobian at the stationary solution if needed
	double j1[2];
	bool ji = false;
	if(cout || cycles_summary) {
		pars.get_sol(ys);
		segment = get_segment(dydt);
		if(segment == target_segment) last_t = t; // note: this is expected
		dst[0] = y1[0] - ys[0];
		dst[1] = y1[1] - ys[1];
		
		ji = pars.get_jac_comp(j1);
		if(!ji) {
			// real eigenvalues, we only care about the larger one
			// (which should still be negative)
			j1[0] += j1[1];
			j1[1] = 0.0;
		}
		if(!cycles_summary)
			fprintf(stderr, "%f\t%f\t%f\t%f\n", j1[0], j1[1], ys[0], ys[1]);
	}
	
	while(tmax > 0.0 && t < tmax) {
		// add noise sampled from a Gaussian to either or both variables
		if(noise_x > 0.0) y1[0] = std::max(0.0, y1[0] + nx(rng) * y1[0]);
		if(noise_y > 0.0) y1[1] = std::max(0.0, y1[1] + ny(rng) * y1[1]);
		
		if(fixed_step) {
			for(int i = 0; i < 2; i++) if(discrete_noise && dydt[i] != 0.0) {
				bool xn = false;
				if(dydt[i] < 0.0) {
					xn = true;
					dydt[i] *= -1.0;
				}
				std::exponential_distribution<double> ed1(1.0 / dydt[i]);
				dydt[i] = ed1(rng);
				if(xn) dydt[i] *= -1.0;
			}
			y1[0] += dydt[0];
			y1[1] += dydt[1];
			if(y1[0] < 0.0) y1[0] = 0.0; // all variables are constrained to >= 0
			if(y1[1] < 0.0) y1[1] = 0.0;
			if(pars.limit_x() && y1[0] > pars.K) y1[0] = pars.K;
			if(pars.limit_y() && y1[1] > pars.K) y1[1] = pars.K;
			t += dt;
		}
		else {
			if(noise_x > 0.0 || noise_y > 0.0) gsl_odeiv2_driver_reset(d);
			double t1 = t + dt;
			int r = gsl_odeiv2_driver_apply(d, &t, t1, y1);
			if(r != GSL_SUCCESS) {
				fprintf(stderr, "Error advancing the solution (%d)!\n", r);
				gsl_odeiv2_driver_free(d);
				return r;
			}
		}
		
		if(cycles_summary && y1[0] > ys[0] && y1[0] > overshoot) overshoot = y1[0];
		
		if(fout) fprintf(fout, "%f\t%f\t%f\n", t, y1[0], y1[1]);
		
		// calculate the new derivatives, test if they are in the same segment if needed
		if(cout || cycles_summary || fixed_step) dfunc2(t, y1, dydt, pars);
		if(cout || cycles_summary) {
			int seg2 = get_segment(dydt);
			
			double d2;
			double dst2[] = {0.0, 0.0};
			
			if(seg2 == target_segment || cycles_summary) {
				dst2[0] = y1[0] - ys[0];
				dst2[1] = y1[1] - ys[1];
				d2 = sqrt(dst2[0]*dst2[0] + dst2[1]*dst2[1]);
				if(cycles_summary && d2 < dist_eps) break;
			}
			
			if(seg2 != 0 && seg2 != segment) { // only process derivatives if they are not close to zero
				// check if the progression is "valid", i.e. we didn't skip a segment
				if(seg2 == segment + 1 || (seg2 == 1 && segment == 4)) {
					// valid progression, update segment
					segment = seg2;
				}
				else {
					// invalid progression (probably we are too close to the stationary solution),
					// stop evaluating this
					cout = nullptr;
					if(cycles_summary) break;
				}
				
				if(segment == target_segment) {
					if(last_t >= 0.0) {
						double tmp1 = dst2[0] / dst[0];
						double tmp2 = dst2[1] / dst[1];
						double d1 = sqrt(dst[0]*dst[0] + dst[1]*dst[1]);
						
						double tmp3 = d2 / d1;
						if(cycles_summary) {
							step1 = tmp3;
							break;
						}
						fprintf(cout, "%f\t%f\t%f\t%f\n", t - last_t, tmp1, tmp2, tmp3);
					}
					last_t = t;
					dst[0] = dst2[0];
					dst[1] = dst2[1];
				}
			}
		}
	}
	
	if(d) gsl_odeiv2_driver_free(d);
	
	if(cycles_summary) {
		if(overshoot > 0.0) overshoot = (overshoot - ys[0]) / ys[0];
		FILE* f1 = stdout;
		if(output_params) pars.print_pars(f1, sep, false);
		fprintf(f1, "%2$f%1$c%3$f%1$c%4$f%1$c%5$f\n", sep, j1[0], j1[1], step1 > 0.0 ? step1 : NAN, std::max(overshoot, 0.0));
	}
	
	return 0;
}


int main(int argc, char **argv) {
	params pars;
	double x = 100.0;
	double y = 10.0;
	double tmax = 500.0; // maximum time
	bool tmax_par = false; // whether the above was read as a command line parameter
	double dt = 1.0; // time step for output
	bool fixed_step = false; // if true, instead of running a solver,
	// use discrete-time steps (of e.g. 1 year), essentially evolving
	// the system as a difference equation
	double noise_x = 0.0; // add noise from a Gaussian with this relative scale
	double noise_y = 0.0; // to the x or y variables in each year
	bool discrete_noise = false; // add noise in the discrete model
	// directly to the yearly differences
	uint64_t seed = time(0);
	char* cycles_out = nullptr; // output info about convergence per cycles
	bool cycles_summary = false; // only output some summary info:
	// decrease of distance after the first cycle, eigenvalues and maximum overshoot
	int target_segment = 1; // which segment to consider as the target for convergence analysis
	
	bool auto_a = false; // adjust the a parameter automatically (only works for the BT and BT_LIN model versions)
	bool auto_bd = false; // adjust the b and d parameters automatically (only works for the BT and BT_LIN model versions)
	
	bool multi_run = false; // repeated runs with cycles_summary == true
	double amin, amax, astep = -1.0;
	double bmin, bmax, bstep = -1.0;
	double cmin, cmax, cstep = -1.0;
	
	bool flush = false;
	
	for(int i = 1; i < argc; i++) if(argv[i][0] == '-') switch(argv[i][1]) {
		case 'M':
			{
				multi_run = true;
				double min1 = atof(argv[i+1]);
				double step1 = atof(argv[i+2]);
				double max1 = atof(argv[i+3]);
				switch(argv[i][2]) {
					case 'a':
						amin = min1; amax = max1; astep = step1;
						break;
					case 'b':
						bmin = min1; bmax = max1; bstep = step1;
						break;
					case 'c':
						cmin = min1; cmax = max1; cstep = step1;
						break;
					default:
						fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
						return 1;
				}
				i += 3;
				break;
			}
		case 'a':
			if(argv[i+1][0] == '-' && argv[i+1][1] == 0) auto_a = true;
			else pars.a = atof(argv[i+1]);
			i++;
			break;
		case 'b':
			if(argv[i+1][0] == '-' && argv[i+1][1] == 0) auto_bd = true;
			pars.b = atof(argv[i+1]);
			i++;
			break;
		case 'c':
			pars.c = atof(argv[i+1]);
			i++;
			break;
		case 'd':
			if(argv[i+1][0] == '-' && argv[i+1][1] == 0) auto_bd = true;
			pars.d = atof(argv[i+1]);
			i++;
			break;
		case 'e':
			pars.e = atof(argv[i+1]);
			i++;
			break;
		case 'r':
			pars.r = atof(argv[i+1]);
			i++;
			break;
		case 'K':
			pars.K = atof(argv[i+1]);
			i++;
			break;
		case 't':
			if(!params::try_parse_type(argv[i+1], pars.t)) {
				fprintf(stderr, "Unknown model type: %s!\n", argv[i+1]);
				return 1;
			}
			i++;
			break;
		case 'x':
			x = atof(argv[i+1]);
			i++;
			break;
		case 'y':
			y = atof(argv[i+1]);
			i++;
			break;
		case 'T':
			tmax = atof(argv[i+1]);
			tmax_par = true;
			i++;
			break;
		case 'D':
			dt = atof(argv[i+1]);
			i++;
			break;
		case 'F':
			fixed_step = true;
			break;
		case 'N':
			switch(argv[i][2]) {
				case 'x':
					noise_x = atof(argv[i+1]);
					i++;
					break;
				case 'y':
					noise_y = atof(argv[i+1]);
					i++;
					break;
				case 0:
					noise_x = atof(argv[i+1]);
					noise_y = noise_y;
					i++;
					break;
				case 'd':
					discrete_noise = true;
					break;
				default:
					fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
					return 1;
			}
			break;
		case 's':
			seed = strtoul(argv[i+1], nullptr, 10);
			i++;
			break;
		case 'C':
			switch(argv[i][2]) {
				case 0:
					cycles_out = argv[i+1];
					i++;
					break;
				case 's':
					target_segment = atoi(argv[i+1]);
					i++;
					if(target_segment < 1 || target_segment > 4) {
						fprintf(stderr, "Invalid value for target segment: %s %s!\n", argv[i], argv[i+1]);
						return 1;
					}
					break;
				case 'S':
					cycles_summary = true;
					if(!tmax_par) tmax = -1.0;
					break;
				default:
					fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
					return 1;
			}
			break;
		case 'f':
			flush = true;
			break;
		default:
			fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
			return 1;
	}
	else {
		fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
		return 1;
	}
	
	if(pars.t == params::type::NO_TYPE) {
		fprintf(stderr, "Type of system not given!\n");
		return 1;
	}
	
	if(multi_run) {
		cycles_summary = true;
		switch(pars.t) {
			case params::type::BT:
			case params::type::BT_LIN:
				if((bstep <= 0.0 && astep <= 0.0) || cstep <= 0.0) throw std::runtime_error("Missing parameters!\n");
				if(bstep > 0.0 && astep > 0.0) throw std::runtime_error("Invalid parameter combination!\n");
				break;
			case params::type::SOIL2:
				if(astep <= 0.0 || cstep <= 0.0) throw std::runtime_error("Missing parameters!\n");
				break;
			default:
				throw std::runtime_error("Repeated runs not supported for the given model variant!\n");
		}
	}
	else {
		if(auto_a && auto_bd) throw std::runtime_error("Unsupported parameter combination (either a or b and d must be specified)!\n");
		if(auto_a) pars.set_auto_a();
		if(auto_bd) pars.set_auto_bd();
	}
	
	std::mt19937_64 rng(seed);
	
	FILE* fout = stdout;
	FILE* cout = nullptr;
	if(cycles_summary) fout = nullptr;
	else if(cycles_out) {
		if(cycles_out[0] == '-' && cycles_out[1] == 0) {
			cout = stdout;
			fout = nullptr;
		}
		else {
			cout = fopen(cycles_out, "w");
			if(!cout) {
				fprintf(stderr, "Error opening output file %s!\n", cycles_out);
				return 1;
			}
		}
	}
	
	double y0[] = {x, y};
	
	if(multi_run) {
		pars.print_header(stdout, ',', false);
		fprintf(stdout, "lambda_re,lambda_im,step1,overshoot\n");
		switch(pars.t) {
			case params::type::BT:
			case params::type::BT_LIN:
				if(bstep > 0.0) {
					for(pars.b = bmin; pars.b <= bmax; pars.b += bstep)
					for(pars.c = cmin; pars.c <= cmax; pars.c += cstep) {
						if(auto_a) pars.set_auto_a();
						do_one_run(pars, y0, tmax, dt, fixed_step, noise_x, noise_y, discrete_noise,
							fout, cout, cycles_summary, target_segment, rng, ',', true);
						if(flush) fflush(stdout);
					}
					break;
				}
				else auto_bd = true; // fallthrough in this case
			case params::type::SOIL2:
				{
					for(pars.a = amin; pars.a <= amax; pars.a += astep)
					for(pars.c = cmin; pars.c <= cmax; pars.c += cstep) {
						if(auto_bd) pars.set_auto_bd();
						do_one_run(pars, y0, tmax, dt, fixed_step, noise_x, noise_y, discrete_noise,
							fout, cout, cycles_summary, target_segment, rng, ',', true);
						if(flush) fflush(stdout);
					}
				}
				break;
		}
	}
	else do_one_run(pars, y0, tmax, dt, fixed_step, noise_x, noise_y, discrete_noise, fout, cout,
		cycles_summary, target_segment, rng, '\t', false);
	
	if(cout && cout != stdout) fclose(cout);
	
	return 0;
}

