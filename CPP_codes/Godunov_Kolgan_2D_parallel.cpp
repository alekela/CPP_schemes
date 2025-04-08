#include <iostream>
#include <fstream>
#include <string> 
#include <vector>
#include <map>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <sstream>


namespace fs = std::filesystem;
typedef std::vector<double> vec1d;
typedef std::vector<std::vector<double>> vec2d;


struct InitialState {
	int nx = 5;
	int ny = 5;
	int fict = 2;
	int Nx = nx + 2 * fict;
	int Ny = ny + 2 * fict;

	double x_start = 0;
	double x_end = 1;
	double y_start = 0;
	double y_end = 1;
	double hx = (x_end - x_start) / nx;
	double hy = (y_end - y_start) / ny;

	double t_end = 0.25;
	double time = 0;
	double tau;

	double CFL = 0.05;
	double gamma = 1.4;
	double Q = 2;

    double ux_left, ux_right, uy_left, uy_right, P_left, P_right, rho_left, rho_right;
	double ux_down, ux_up, uy_down, uy_up, P_down, P_up, rho_down, rho_up;
};


// from conservative
void cons_to_noncons(double gamma, double& p, double& vx, double& vy, double& rho, double& m, double& impx, double& impy, double& rhoe) {
	p = (gamma - 1.0) * (rhoe - 0.5 * (pow(impx, 2.0) + pow(impy, 2.0)) * m);
	vx = impx / m;
	vy = impy / m;
	rho = m;
}
// to conservative
void noncons_to_cons(double gamma, double& p, double& vx, double& vy, double& rho, double& m, double& impx, double& impy, double& rhoe) {
	m = rho;
	impx = rho * vx;
	rhoe = 0.5 * rho * (pow(vx, 2.0) + pow(vy, 2.0)) + p / (gamma - 1.0);
}


void grid(InitialState IS, vec1d& x, vec1d& y, vec1d& xc, vec1d& yc) {
	for (int i = 0; i < IS.Nx + 1; ++i) {
		x[i] = IS.x_start + (i - IS.fict) * IS.hx;
	}
	for (int i = 0; i < IS.Nx; ++i) {
		xc[i] = 0.5 * (x[i] + x[i + 1]);
	}
	for (int i = 0; i < IS.Ny + 1; ++i) {
		y[i] = IS.y_start + (i - IS.fict) * IS.hy;
	}

	for (int i = 0; i < IS.Ny; ++i) {
		yc[i] = 0.5 * (y[i] + y[i + 1]);
	}
}


void initial_state(int Nx, int Ny, double gamma, vec1d xc, vec1d yc, vec2d& p, vec2d& vx, vec2d& vy, vec2d& rho, vec2d& m, vec2d& impx, vec2d& impy, vec2d& rhoe) {
	double R = 0.1;
	double p1, vx1, vy1, rho1, p2, vx2, vy2, rho2, d;

	// toro test 1

	p1 = 1.0;
	vx1 = 0.;
	vy1 = 0.;
	rho1 = 1.0;

	p2 = 0.1;
	vx2 = 0.;
	vy2 = 0.;
	rho2 = 0.125;

	for (int i = 0; i < Ny; i++) {
		for (int j = 0; j < Nx; j++) {

			//d = (xc[i] - 0.5) * (xc[i] - 0.5) + (yc[j] - 0.5) * (yc[j] - 0.5);
			//if (d < R * R) {
			if (j < Nx / 2) {

				p[i][j] = p1;
				vx[i][j] = vx1;
				vy[i][j] = vy1;
				rho[i][j] = rho1;
			}
			else {
				p[i][j] = p2;
				vx[i][j] = vx2;
				vy[i][j] = vy2;
				rho[i][j] = rho2;
			}
			noncons_to_cons(gamma, p[i][j], vx[i][j], vy[i][j], rho[i][j], m[i][j], impx[i][j], impy[i][j], rhoe[i][j]);
		}
	}
}



// minmod functions: 0 - Kolgan,1972; 1 - Kolgan,1975; 2 (or another...) - Osher,1984
// func - type of minmod
double minmod(double a, double b, int func) {
	if (func == 0)
		return (abs(a) < abs(b) ? a : b);
	else if (func == 1)
	{
		double c = (a + b) / 2;
		return ((abs(a) < abs(b) && abs(a) < abs(c)) ? a : (abs(b) < abs(a) && abs(b) < abs(c) ? b : c));
	}
	else
		return (a * b < 0 ? 0 : (a * a < a * b ? a : b));
}


// is vacuum created?
bool is_vacuum(double gamma, double CL, double UL, double CR, double UR) {
	if (2.0 / (gamma - 1.0) * (CL + CR) < UR - UL) return true;
	return false;
}


// in F calculate U on contact, in FD calculate dU/dP on contact
void prefun(double gamma, double& F, double& FD, double P, double DK, double PK, double CK) {

	if (P <= PK) //rarefaction wave
	{
		F = 2. / (gamma - 1) * CK * (pow((P / PK), (gamma - 1.) / 2. / gamma) - 1.);
		FD = CK * gamma / PK * pow((P / PK), (gamma - 1) / gamma / 2. - 1);
	}
	else //shock wave
	{
		F = (P - PK) / sqrt(DK / 2 * ((gamma + 1) * P + (gamma - 1) * PK));
		FD = 1 / sqrt(2 * DK) * ((gamma + 1) * P + (3 * gamma - 1) * P) / pow(((gamma + 1) * P + (gamma - 1) * PK), 1.5);
	}
}


void Newton_find_P(double gamma, double Quser, double PL, double DL, double UL, double PR, double DR, double UR, double& P) {
	double CUP, GEL, GER, PMAX, PMIN, PPV, QMAX, EPS = 1.e-8;
	double CL = sqrt(gamma * PL / DL);
	double CR = sqrt(gamma * PR / DR);

	CUP = 0.25 * (DL + DR) * (CL + CR);

	PPV = 0.5 * (PL + PR) + 0.5 * (UL - UR) * CUP; // initial condition from linear task
	PPV = std::max(0.0, PPV);

	PMIN = std::min(PL, PR);
	PMAX = std::max(PL, PR);

	QMAX = PMAX / PMIN;

	double PM, UM;

	if (QMAX <= Quser && (PMIN <= PPV && PPV <= PMAX) || abs(PMIN - PPV) < EPS || abs(PMAX - PPV) < EPS) {
		PM = PPV;
	}
	else {
		if (PPV < PMIN) // two rarefaction waves
		{
			double gamma1 = 0.5 * (gamma - 1.0) / gamma;
			PM = pow(abs(CL + CR - 0.5 * (gamma - 1.0) * (UR - UL)) / abs(CL / pow(PL, gamma1) + CR / pow(PR, gamma1)), 1.0 / gamma1);
		}
		else //two shock waves
		{
			double gamma1 = 2.0 / (gamma + 1.0);
			double gamma2 = (gamma - 1.) / (gamma + 1.);

			GEL = sqrt(gamma1 / DL / (gamma2 * PL + PPV));
			GER = sqrt(gamma1 / DR / (gamma2 * PR + PPV));
			PM = abs(GEL * PL + GER * PR - (UR - UL)) / (GEL + GER);
		}
	}
	P = PM;
}


//  to compute the solution for pressure and velocity on contact
void Riemann_solver(double gamma, double Quser, double PL, double DL, double UL, double PR, double DR, double UR, double& UM, double& PM) {
	int NRITER = 30;
	double CHANGE = 1e6, FL, FLD, FR, FRD, POLD, TOLPRE = 1.0e-6, UDIFF;
	double CL = sqrt(gamma * PL / DL);
	double CR = sqrt(gamma * PR / DR);

	//std::cout << PL << " " << PR << std::endl;

	if (is_vacuum(gamma, CL, UL, CR, UR)) {
		std::cout << "Vacuum is generated" << std::endl << "Program stopped" << std::endl;
		exit(228);
	}
	Newton_find_P(gamma, Quser, PL, DL, UL, PR, DR, UR, POLD);		// initial condition for pressure

	if (POLD < 0.)
	{
		std::cout << "Pressure is negative" << std::endl;
		exit(105);
	}

	if (POLD != POLD)
	{
		std::cout << "Pressure is NaN after guess_p" << std::endl;
		exit(-1323415);
	}

	UDIFF = UR - UL;

	for (int i = 0; i < NRITER; ++i)
	{
		prefun(gamma, FL, FLD, POLD, DL, PL, CL);
		prefun(gamma, FR, FRD, POLD, DR, PR, CR);

		if (POLD < 0.)
		{
			std::cout << "Pressure is negative" << std::endl;
			exit(106);
		}
		if (POLD != POLD)
		{
			std::cout << "Pressure is NaN in iter" << std::endl;
			exit(1323415);
		}

		PM = POLD - (FL + FR + UDIFF) / (FLD + FRD);
		CHANGE = 2.0 * abs((PM - POLD) / (PM + POLD));
		if (CHANGE < TOLPRE) break;
	}
	UM = 0.5 * (UL + UR + FR - FL);

	return;
}


// to sample the solution throughout the wave pattern, PM, UM - contact's parameters
void find_in_which_part(double gamma, double PL, double DL, double UL, double PR, double DR, double UR,
	double UM, double PM, double S, double& D, double& U, double& P) {
	double C, CML, CMR, PML, PMR, SHL, SHR, SL, SR, STL, STR;
	//
	// C - local sound speed in rarefaction wave
	// CML, CMR - sound speed on the left / right side of contact
	// SHL, SHR - speeds of head of left/right rarefaction wave
	// SL, SR - left/right shock speed
	// STL, STR - speeds of tail of left/right rarefaction wave
	//
	double CL = sqrt(gamma * PL / DL);
	double CR = sqrt(gamma * PR / DR);

	if (S <= UM) { // point is left from contact
		if (PM <= PL) { //left rarefaction wave
			SHL = UL - CL;
			if (S <= SHL) { //point is left from rarefaction wave - left parameters
				D = DL;
				U = UL;
				P = PL;
			}
			else {
				CML = CL * pow((PM / PL), (gamma - 1.0) / (2.0 * gamma));
				STL = UM - CML;
				if (S > STL) { // point is a state with *
					D = DL * pow((PM / PL), (1. / gamma));
					U = UM;
					P = PM;
				}
				else { //point is in a left rarefaction wave
					U = 2.0 / (gamma + 1.0) * (CL + (gamma - 1.0) / 2.0 * UL + S);
					C = 2.0 / (gamma + 1.0) * (CL + (gamma - 1.0) / 2.0 * (UL - S));
					D = DL * pow((C / CL), 2.0 / (gamma - 1.0));
					P = PL * pow((C / CL), 2.0 * gamma / (gamma - 1.0));
				}
			}
		}
		else { // left shock wave
			PML = PM / PL;
			SL = UL - CL * sqrt((gamma + 1.0) / (2.0 * gamma) * PML + (gamma - 1.0) / (2.0 * gamma));
			if (S <= SL) { //point is left from shock wave - left parameters
				D = DL;
				U = UL;
				P = PL;
			}
			else { // point is a state with *
				D = DL * (PML + (gamma - 1.0) / (gamma + 1.0)) / (PML * (gamma - 1.0) / (gamma + 1.0) + 1.);
				U = UM;
				P = PM;
			}
		}
	}
	else { // point is right from contact
		if (PM > PR) { //right from shock wave
			PMR = PM / PR;
			SR = UR + CR * sqrt((gamma + 1.0) / (2.0 * gamma) * PMR + (gamma - 1.0) / (2.0 * gamma));
			if (S >= SR) { //point is right from shock wave - right parameters
				D = DR;
				U = UR;
				P = PR;
			}
			else { // point is a state with *
				D = DR * (PMR + (gamma - 1.0) / (gamma + 1.0)) / (PMR * (gamma - 1.0) / (gamma + 1.0) + 1.);
				U = UM;
				P = PM;
			}
		}
		else // right rarefaction wave
		{
			SHR = UR + CR;
			if (S >= SHR) //point is right from rarefaction wave - right parameters
			{
				D = DR;
				U = UR;
				P = PR;
			}
			else
			{
				CMR = CR * pow((PM / PR), (gamma - 1.0) / (2.0 * gamma));
				STR = UM + CMR;
				if (S <= STR) // point is a state with *
				{
					D = DR * pow((PM / PR), (1. / ((gamma - 1.0) / (2.0 * gamma))));
					U = UM;
					P = PM;
				}
				else // point is in a right rarefaction wave
				{
					U = 2.0 / (gamma + 1.0) * (-CR + (gamma - 1.0) / 2.0 * UR + S);
					C = 2.0 / (gamma + 1.0) * (CR - (gamma - 1.0) / 2.0 * (UR - S));
					D = DR * pow((C / CR), 2.0 / (gamma - 1.0));
					P = PR * pow((C / CR), 2.0 * gamma / (gamma - 1.0));
				}
			}
		}
	}
	return;
}


// Riemann solver
void Riemann_solver(double gamma, double Quser, double PL, double DL, double UL, double PR, double DR, double UR,
	double& D, double& U, double& P) {
	double CL = sqrt(gamma * PL / DL); // left sound speed
	double CR = sqrt(gamma * PR / DR); // right sound speed
	double PM, UM; // pressure and velocity on contact
	double S = 0.;

	/*if (DL <= 0.0 || DR <= 0.0 || PL <= 0.0 || PR <= 0.0) {
		std::cerr << "Invalid input to Riemann solver" << std::endl;
		exit(1);
	}*/

	// vacuum generation
	if (is_vacuum(gamma, CL, UL, CR, UR)) {
		std::cout << "Vacuum is generated" << std::endl << "Program stopped" << std::endl;
		exit(228);
	}
	Riemann_solver(gamma, Quser, PL, DL, UL, PR, DR, UR, UM, PM);

	// found results
	find_in_which_part(gamma, PL, DL, UL, PR, DR, UR, UM, PM, 0., D, U, P);

	if (P <= 0.0 || D <= 0.0 || !std::isfinite(P) || !std::isfinite(U)) {
		std::cerr << "Invalid output from Riemann solver" << std::endl;
		exit(1);
	}
}

void Godunov_flux_x(double gamma, double p, double vx, double vy, double rho, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
	double m, impx, impy, rhoe;
	noncons_to_cons(gamma, p, vx, vy, rho, m, impx, impy, rhoe);
	Fm = rho * vx;
	Fimpx = Fm * vx + p;
	Fimpy = Fm * vy;
	Fe = (p + rhoe) * vx;
}


void Godunov_flux_y(double gamma, double p, double vx, double vy, double rho, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
	double m, impx, impy, rhoe;
	noncons_to_cons(gamma, p, vx, vy, rho, m, impx, impy, rhoe);
	Fm = rho * vy;
	Fimpx = Fm * vx;
	Fimpy = Fm * vy + p;
	Fe = (p + rhoe) * vy;
}


void Godunov_method_x(double gamma, double Quser, double ml, double impxl, double impyl, double el,
	double mr, double impxr, double impyr, double er,
	double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
	double p, vx, vy, rho;
	double pl, vxl, vyl, rhol;
	double pr, vxr, vyr, rhor;

	cons_to_noncons(gamma, pl, vxl, vyl, rhol, ml, impxl, impyl, el);
	cons_to_noncons(gamma, pr, vxr, vyr, rhor, mr, impxr, impyr, er);
	Riemann_solver(gamma, Quser, pl, rhol, vxl, pr, rhor, vxr, rho, vx, p);

	if (vx >= 0)
		vy = impyl / ml;
	else
		vy = impyr / mr;

	Godunov_flux_x(gamma, p, vx, vy, rho, Fm, Fimpx, Fimpy, Fe);
}


void Godunov_method_y(double gamma, double Quser, double md, double impxd, double impyd, double ed,
	double mu, double impxu, double impyu, double eu,
	double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
	double p, vx, vy, rho;
	double pd, vxd, vyd, rhod;
	double pu, vxu, vyu, rhou;

	cons_to_noncons(gamma, pd, vxd, vyd, rhod, md, impxd, impyd, ed);
	cons_to_noncons(gamma, pu, vxu, vyu, rhou, mu, impxu, impyu, eu);
	Riemann_solver(gamma, Quser, pd, rhod, vyd, pu, rhou, vyu, rho, vy, p);

	if (vy >= 0)
		vx = impxd / md;
	else
		vx = impxu / mu;

	//Godunov_flux_x(double gamma, double p, double vx, double vy, double rho, double& Fm, double& Fimpx, double& Fimpy, double& Fe)
	Godunov_flux_y(gamma, p, vx, vy, rho, Fm, Fimpx, Fimpy, Fe);
}


void Boundary_x(double p, double vx, double vy, double rho, double& pb, double& vxb, double& vyb, double& rhob, int mode) {
	// wall
	if (mode == 0)
	{
		pb = p;
		vxb = -vx;
		vyb = vy;
		rhob = rho;
	}
	// free flux
	else
	{
		pb = p;
		vxb = vx;
		vyb = vy;
		rhob = rho;
	}
}

void Boundary_y(double p, double vx, double vy, double rho, double& pb, double& vxb, double& vyb, double& rhob, int mode) {
	// wall
	if (mode == 0)
	{
		pb = p;
		vxb = vx;
		vyb = -vy;
		rhob = rho;
	}
	// free flux
	else
	{
		pb = p;
		vxb = vx;
		vyb = vy;
		rhob = rho;
	}
}

void Boundary(InitialState IS, vec2d& P, vec2d& ux, vec2d& uy, vec2d& rho, int b_l, int b_u, int b_r, int b_d) {
	if (b_l > -1) {
		for (int i = IS.fict; i < IS.Ny - IS.fict; i++) {
			for (int j = 0; j < IS.fict; j++) {
				if (b_l == 0) {
					P[i][j] = IS.P_left;
					ux[i][j] = IS.ux_left;
					uy[i][j] = IS.uy_left;
					rho[i][j] = IS.rho_left;
				}
				else if (b_l == 1) {
					P[i][j] = P[i][IS.fict];
					ux[i][j] = ux[i][IS.fict];
					uy[i][j] = uy[i][IS.fict];
					rho[i][j] = rho[i][IS.fict];
				}
			}
		}
	}
	if (b_u > -1) {
		for (int i = IS.Ny - IS.fict; i < IS.Ny; i++) {
			for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
				if (b_u == 0) {
					P[i][j] = IS.P_up;
					ux[i][j] = IS.ux_up;
					uy[i][j] = IS.uy_up;
					rho[i][j] = IS.rho_up;
				}
				else if (b_u == 1) {
					P[i][j] = P[IS.Ny - IS.fict - 1][j];
					ux[i][j] = ux[IS.Ny - IS.fict - 1][j];
					uy[i][j] = uy[IS.Ny - IS.fict - 1][j];
					rho[i][j] = rho[IS.Ny - IS.fict - 1][j];
				}
			}
		}
	}
	if (b_r > -1) {
		for (int i = IS.fict; i < IS.Ny - IS.fict; i++) {
			for (int j = IS.Nx - IS.fict; j < IS.Nx; j++) {
				if (b_r == 0) {
					P[i][j] = IS.P_right;
					ux[i][j] = IS.ux_right;
					uy[i][j] = IS.uy_right;
					rho[i][j] = IS.rho_right;
				}
				else if (b_r == 1) {
					P[i][j] = P[i][IS.Nx - IS.fict - 1];
					ux[i][j] = ux[i][IS.Nx - IS.fict - 1];
					uy[i][j] = uy[i][IS.Nx - IS.fict - 1];
					rho[i][j] = rho[i][IS.Nx - IS.fict - 1];
				}
			}
		}
	}
	if (b_d > -1) {
		for (int i = 0; i < IS.fict; i++) {
			for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
				if (b_d == 0) {
					P[i][j] = IS.P_down;
					ux[i][j] = IS.ux_down;
					uy[i][j] = IS.uy_down;
					rho[i][j] = IS.rho_down;
				}
				else if (b_d == 1) {
					P[i][j] = P[IS.fict][j];
					ux[i][j] = ux[IS.fict][j];
					uy[i][j] = uy[IS.fict][j];
					rho[i][j] = rho[IS.fict][j];
				}
			}
		}
	}
}


void Godunov_Kolgan_solver_2D_one_step(InitialState IS, double* t_start, vec1d& x, vec1d& y, 
	vec2d& P, vec2d& ux, vec2d& uy, vec2d& rho, vec2d& mass, vec2d& Imp_x, vec2d& Imp_y, vec2d& rhoe, 
	int* neighbors, int bound_left, int bound_up, int bound_right, int bound_down, int minmod_type) {

	vec2d new_mass(IS.Ny, vec1d(IS.Nx)), new_Imp_x(IS.Ny, vec1d(IS.Nx)), new_Imp_y(IS.Ny, vec1d(IS.Nx)), new_rhoe(IS.Ny, vec1d(IS.Nx));
	double mb, impxb, impyb, eb, pb, vxb, vyb, rhob;
	double FmL, FimpxL, FimpyL, FeL, FmR, FimpxR, FimpyR, FeR;

	double pl, vl, rhol, pr, vr, rhor;
	double massl, impxl, impyl, rhoel, massr, impxr, impyr, rhoer;
	double dm, dimp, de;

	// local boundary conditions for each proc
	int b_l, b_r, b_u, b_d;
	b_l = neighbors[0] > 0 ? -1 : bound_left;
	b_u = neighbors[1] > 0 ? -1 : bound_up;
	b_r = neighbors[2] > 0 ? -1 : bound_right;
	b_d = neighbors[3] > 0 ? -1 : bound_down;


	double tau = 1.e6;
	double c, temp;
	for (int i = IS.fict; i < IS.Nx - IS.fict; ++i) {
		for (int j = IS.fict; j < IS.Ny - IS.fict; ++j) {
			c = sqrt(IS.gamma * P[i][j] / rho[i][j]);
			temp = std::min(abs(IS.CFL * (x[i + 1] - x[i]) / (abs(ux[i][j]) + c)), abs(IS.CFL * (y[j + 1] - y[j]) / (abs(uy[i][j]) + c)));
			if (temp < tau) {
				tau = temp;
			}
		}
	}
	(*t_start) += tau;

	Boundary(IS, P, ux, uy, rho, b_l, b_u, b_r, b_d);

	for (int i = IS.fict; i < IS.Ny - IS.fict + 1; i++) {
		for (int j = IS.fict; j < IS.Nx - IS.fict + 1; j++) {
			massl = mass[i - 1][j] + 0.5 * minmod(mass[i - 1][j] - mass[i - 2][j], mass[i][j] - mass[i - 1][j], minmod_type);
			impxl = Imp_x[i - 1][j] + 0.5 * minmod(Imp_x[i - 1][j] - Imp_x[i - 2][j], Imp_x[i][j] - Imp_x[i - 1][j], minmod_type);
			impyl = Imp_y[i - 1][j] + 0.5 * minmod(Imp_y[i - 1][j] - Imp_y[i - 2][j], Imp_y[i][j] - Imp_y[i - 1][j], minmod_type);
			rhoel = rhoe[i - 1][j] + 0.5 * minmod(rhoe[i - 1][j] - rhoe[i - 2][j], rhoe[i][j] - rhoe[i - 1][j], minmod_type);
			
			massr = mass[i][j] - 0.5 * minmod(mass[i][j] - mass[i - 1][j], mass[i + 1][j] - mass[i][j], minmod_type);
			impxr = Imp_x[i][j] - 0.5 * minmod(Imp_x[i][j] - Imp_x[i - 1][j], Imp_x[i + 1][j] - Imp_x[i][j], minmod_type);
			impyr = Imp_y[i][j] - 0.5 * minmod(Imp_y[i][j] - Imp_y[i - 1][j], Imp_y[i + 1][j] - Imp_y[i][j], minmod_type);
			rhoer = rhoe[i][j] - 0.5 * minmod(rhoe[i][j] - rhoe[i - 1][j], rhoe[i + 1][j] - rhoe[i][j], minmod_type);

			Godunov_method_y(IS.gamma, IS.Q, massl, impxl, impyl, rhoel, massr, impxr, impyr, rhoer, FmL, FimpxL, FimpyL, FeL);

			/*
			massl = mass[i][j] + 0.5 * minmod(mass[i][j] - mass[i - 1][j], mass[i + 1][j] - mass[i][j], minmod_type);
			impxl = Imp_x[i][j] + 0.5 * minmod(Imp_x[i][j] - Imp_x[i - 1][j], Imp_x[i + 1][j] - Imp_x[i][j], minmod_type);
			impyl = Imp_y[i][j] + 0.5 * minmod(Imp_y[i][j] - Imp_y[i - 1][j], Imp_y[i + 1][j] - Imp_y[i][j], minmod_type);
			rhoel = rhoe[i][j] + 0.5 * minmod(rhoe[i][j] - rhoe[i - 1][j], rhoe[i + 1][j] - rhoe[i][j], minmod_type);

			massr = mass[i + 1][j] - 0.5 * minmod(mass[i + 1][j] - mass[i][j], mass[i + 2][j] - mass[i + 1][j], minmod_type);
			impxr = Imp_x[i + 1][j] - 0.5 * minmod(Imp_x[i + 1][j] - Imp_x[i][j], Imp_x[i + 2][j] - Imp_x[i + 1][j], minmod_type);
			impyr = Imp_y[i + 1][j] - 0.5 * minmod(Imp_y[i + 1][j] - Imp_y[i][j], Imp_y[i + 2][j] - Imp_y[i + 1][j], minmod_type);
			rhoer = rhoe[i + 1][j] - 0.5 * minmod(rhoe[i + 1][j] - rhoe[i][j], rhoe[i + 2][j] - rhoe[i + 1][j], minmod_type);
			Godunov_method_y(IS.gamma, IS.Q, massl, impxl, impyl, rhoel, massr, impxr, impyr, rhoer, FmR, FimpxR, FimpyR, FeR);
			*/

			//std::cout << FmL << " " << FmR << std::endl;

			new_mass[i][j] = mass[i][j] - tau * (FmR - FmL) / IS.hx;
			new_Imp_x[i][j] = Imp_x[i][j] - tau * (FimpxR - FimpxL) / IS.hx;
			new_Imp_y[i][j] = Imp_y[i][j] - tau * (FimpyR - FimpyL) / IS.hx;
			new_rhoe[i][j] = rhoe[i][j] - tau * (FeR - FeL) / IS.hx;

			
			massl = mass[i][j - 1] + 0.5 * minmod(mass[i][j - 1] - mass[i][j - 2], mass[i][j] - mass[i][j - 1], minmod_type);
			impxl = Imp_x[i][j - 1] + 0.5 * minmod(Imp_x[i][j - 1] - Imp_x[i][j - 2], Imp_x[i][j] - Imp_x[i][j - 1], minmod_type);
			impyl = Imp_y[i][j - 1] + 0.5 * minmod(Imp_y[i][j - 1] - Imp_y[i][j - 2], Imp_y[i][j] - Imp_y[i][j - 1], minmod_type);
			rhoel = rhoe[i][j - 1] + 0.5 * minmod(rhoe[i][j - 1] - rhoe[i][j - 2], rhoe[i][j] - rhoe[i][j - 1], minmod_type);

			massr = mass[i][j] - 0.5 * minmod(mass[i][j] - mass[i][j - 1], mass[i][j + 1] - mass[i][j], minmod_type);
			impxr = Imp_x[i][j] - 0.5 * minmod(Imp_x[i][j] - Imp_x[i][j - 1], Imp_x[i][j + 1] - Imp_x[i][j], minmod_type);
			impyr = Imp_y[i][j] - 0.5 * minmod(Imp_y[i][j] - Imp_y[i][j - 1], Imp_y[i][j + 1] - Imp_y[i][j], minmod_type);
			rhoer = rhoe[i][j] - 0.5 * minmod(rhoe[i][j] - rhoe[i][j - 1], rhoe[i][j + 1] - rhoe[i][j], minmod_type);
			Godunov_method_x(IS.gamma, IS.Q, massl, impxl, impyl, rhoel, massr, impxr, impyr, rhoer, FmL, FimpxL, FimpyL, FeL);

			/*
			massl = mass[i][j] + 0.5 * minmod(mass[i][j] - mass[i][j - 1], mass[i][j + 1] - mass[i][j], minmod_type);
			impxl = Imp_x[i][j] + 0.5 * minmod(Imp_x[i][j] - Imp_x[i][j - 1], Imp_x[i][j + 1] - Imp_x[i][j], minmod_type);
			impyl = Imp_y[i][j] + 0.5 * minmod(Imp_y[i][j] - Imp_y[i][j - 1], Imp_y[i][j + 1] - Imp_y[i][j], minmod_type);
			rhoel = rhoe[i][j] + 0.5 * minmod(rhoe[i][j] - rhoe[i][j - 1], rhoe[i][j + 1] - rhoe[i][j], minmod_type);

			massr = mass[i][j + 1] - 0.5 * minmod(mass[i][j + 1] - mass[i][j], mass[i][j + 2] - mass[i][j + 1], minmod_type);
			impxr = Imp_x[i][j + 1] - 0.5 * minmod(Imp_x[i][j + 1] - Imp_x[i][j], Imp_x[i][j + 2] - Imp_x[i][j + 1], minmod_type);
			impyr = Imp_y[i][j + 1] - 0.5 * minmod(Imp_y[i][j + 1] - Imp_y[i][j], Imp_y[i][j + 2] - Imp_y[i][j + 1], minmod_type);
			rhoer = rhoe[i][j + 1] - 0.5 * minmod(rhoe[i][j + 1] - rhoe[i][j], rhoe[i][j + 2] - rhoe[i][j + 1], minmod_type);
			Godunov_method_x(IS.gamma, IS.Q, massl, impxl, impyl, rhoel, massr, impxr, impyr, rhoer, FmR, FimpxR, FimpyR, FeR);
			*/

			new_mass[i][j] = new_mass[i][j] - tau * (FmR - FmL) / IS.hy;
			new_Imp_x[i][j] = new_Imp_x[i][j] - tau * (FimpxR - FimpxL) / IS.hy;
			new_Imp_y[i][j] = new_Imp_y[i][j] - tau * (FimpyR - FimpyL) / IS.hy;
			new_rhoe[i][j] = new_rhoe[i][j] - tau * (FeR - FeL) / IS.hy;

		}
	}

	// update parameters
	for (int i = 0; i < IS.Nx; ++i) {
		for (int j = 0; j < IS.Ny; ++j) {
			mass[i][j] = new_mass[i][j];
			Imp_x[i][j] = new_Imp_x[i][j];
			Imp_y[i][j] = new_Imp_y[i][j];
			rhoe[i][j] = new_rhoe[i][j];

			cons_to_noncons(IS.gamma, P[i][j], ux[i][j], uy[i][j], rho[i][j], mass[i][j], Imp_x[i][j], Imp_y[i][j], rhoe[i][j]);
		}
	}
	return;
}


void mpi_get_px_py(int p, int Nx, int Ny, int* px, int* py) {
	if (p == 1) {
		(*px) = 1;
		(*py) = 1;
	}
	else {
		int best_px, best_py;
		int best_p;
		int mera = -1;
		for (int i = 2; i < p + 1; i++) {
			int first, second;
			for (int j = 1; j < sqrt(i) + 1; j++) {
				if (i % j == 0) {
					first = j;
					second = i / j;
				}
			}
			if (i - abs(second - first) > mera) {
				mera = i - abs(second - first);
				best_px = first;
				best_py = second;
				best_p = i;
			}
		}
		if ((Nx % best_px) * (Ny % best_py) < (Nx % best_py) * (Ny % best_px)) {
			(*px) = best_px;
			(*py) = best_py;
		}
		else {
			(*px) = best_py;
			(*py) = best_px;
		}
	}
}


void sending(InitialState IS, vec2d& P, vec2d& ux, vec2d& uy, vec2d& rho, int* neighbors, int rank, MPI_Status* Status) {
	int length;
	int recv_len;
	vec1d send_arr;
	for (int d = 0; d < 4; d++) {
		if (neighbors[d] > 0) {
			// packing data to send
			if (d == 0) {
				length = 4 * (IS.Ny - 2 * IS.fict) * IS.fict;
				send_arr.reserve(length);
				for (int i = IS.fict; i < IS.Ny - IS.fict; i++) {
					for (int j = IS.fict; j < IS.fict + IS.fict; j++) {
						for (int k = 0; k < 4; k++) {
							if (k == 0) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = P[i][j];
							if (k == 1) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = ux[i][j];
							if (k == 2) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = uy[i][j];
							if (k == 3) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = rho[i][j];
						}
					}
				}
			}
			else if (d == 1) {
				length = 4 * (IS.Nx - 2 * IS.fict) * IS.fict;
				send_arr.reserve(length);
				for (int i = IS.Ny - 2 * IS.fict; i < IS.Ny - IS.fict; i++) {
					for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
						for (int k = 0; k < 4; k++) {
							if (k == 0) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = P[i][j];
							if (k == 1) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = ux[i][j];
							if (k == 2) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = uy[i][j];
							if (k == 3) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = rho[i][j];
						}
					}
				}
			}
			else if (d == 2) {
				length = 4 * (IS.Ny - 2 * IS.fict) * IS.fict;
				send_arr.reserve(length);
				for (int i = IS.fict; i < IS.Ny - IS.fict; i++) {
					for (int j = IS.Nx - 2 * IS.fict; j < IS.Nx - IS.fict; j++) {
						for (int k = 0; k < 4; k++) {
							if (k == 0) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = P[i][j];
							if (k == 1) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = ux[i][j];
							if (k == 2) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = uy[i][j];
							if (k == 3) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = rho[i][j];
						}
					}
				}
			}
			else if (d == 3) {
				length = 4 * (IS.Nx - 2 * IS.fict) * IS.fict;
				send_arr.reserve(length);
				for (int i = IS.fict; i < IS.fict + IS.fict; i++) {
					for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
						for (int k = 0; k < 4; k++) {
							if (k == 0) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = P[i][j];
							if (k == 1) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = ux[i][j];
							if (k == 2) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = uy[i][j];
							if (k == 3) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = rho[i][j];
						}
					}
				}
			}
			MPI_Send(&(send_arr[0]), length, MPI_DOUBLE, neighbors[d], neighbors[d], MPI_COMM_WORLD);
			MPI_Recv(&(recv_len), 1, MPI_INT, neighbors[d], rank, MPI_COMM_WORLD, Status);
			std::vector<double> recv_arr;
			MPI_Recv(&(recv_arr[0]), recv_len, MPI_DOUBLE, neighbors[d], rank, MPI_COMM_WORLD, Status);
			// unpacking received data
			if (d == 0) {
				if (recv_len != 4 * (IS.Ny - 2 * IS.fict) * IS.fict) {
					std::cerr << "Error in sending" << std::endl;
				}
				for (int i = IS.fict; i < IS.Ny - IS.fict; i++) {
					for (int j = IS.Nx - IS.fict; j < IS.Nx; j++) {
						for (int k = 0; k < 4; k++) {
							if (k == 0) P[i][j] = recv_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
							if (k == 1) ux[i][j] = recv_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
							if (k == 2) uy[i][j] = recv_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
							if (k == 3) rho[i][j] = recv_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
						}
					}
				}
			}
			else if (d == 1) {
				if (recv_len != 4 * (IS.Nx - 2 * IS.fict) * IS.fict) {
					std::cerr << "Error in sending" << std::endl;
				}
				for (int i = 0; i < IS.fict; i++) {
					for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
						for (int k = 0; k < 4; k++) {
							if (k == 0) P[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
							if (k == 1) ux[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
							if (k == 2) uy[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
							if (k == 3) rho[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
						}
					}
				}
			}
			else if (d == 2) {
				if (recv_len != 4 * (IS.Ny - 2 * IS.fict) * IS.fict) {
					std::cerr << "Error in sending" << std::endl;
				}
				for (int i = IS.fict; i < IS.Ny - IS.fict; i++) {
					for (int j = 0; j < IS.fict; j++) {
						for (int k = 0; k < 4; k++) {
							if (k == 0) P[i][j] = send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
							if (k == 1) ux[i][j] = send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
							if (k == 2) uy[i][j] = send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
							if (k == 3) rho[i][j] = send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
						}
					}
				}
			}
			else if (d == 3) {
				if (recv_len != 4 * (IS.Nx - 2 * IS.fict) * IS.fict) {
					std::cerr << "Error in sending" << std::endl;
				}
				for (int i = IS.Ny - IS.fict; i < IS.Ny; i++) {
					for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
						for (int k = 0; k < 4; k++) {
							if (k == 0) P[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
							if (k == 1) ux[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
							if (k == 2) uy[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
							if (k == 3) rho[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
						}
					}
				}
			}
		}
	}
}


void writeCSV(std::string filename, int iter, vec1d x, vec1d y, vec2d ux, vec2d uy, vec2d P, vec2d rho, double t, int Nx, int Ny, int fict) {
	std::string name = "CSVs\\";
	name += filename;
	name += "\\Iter=";
	name += std::to_string(iter);
	name += ".csv";

	std::ofstream outfile(name);
	std::string tmp;
	tmp = "Time,X,Y,Rho,P,Ux,Uy\n";
	outfile << tmp;
	for (int i = fict; i < Ny - fict; i++) {
		for (int j = fict; j < Nx - fict; j++) {
			tmp = "";
			tmp += std::to_string(t) + ',';
			tmp += std::to_string(x[j]) + ',';
			tmp += std::to_string(y[i]) + ',';
			tmp += std::to_string(rho[i][j]) + ',';
			tmp += std::to_string(P[i][j]) + ',';
			tmp += std::to_string(ux[i][j]) + ',';
			tmp += std::to_string(uy[i][j]) + '\n';
			outfile << tmp;
		}
	}
	outfile.close();
}


void writeCSV_p(std::string filename, InitialState& IS, int iter, vec1d& xc, vec1d& yc, vec2d& ux, vec2d& uy, vec2d& P, vec2d& rho,
				double _time, int myrank, MPI_Comm comm) {
	if (myrank == 0) {
		if (!fs::exists(out_dir)) {
			fs::create_directory(out_dir);
		}
	}
	MPI_Barrier(comm);

	std::string name = "CSVs\\";
	name += filename;
	name += "\\Iter=";
	name += std::to_string(iter);
	name += ".csv";

	std::ostringstream buffer;

	if (myrank == 0) {
		buffer << (std::to_string(_time) + "\n") << "X,Y,Rho,P,Ux,Uy\n";
	}
	for (int i = IS.fict; i < IS.Ny - IS.fict; ++i) {
		for (int j = IS.fict; j < IS.Nx - IS.fict; ++j) {
			buffer << xc[i] << "," << yc[j] << "," << rho[i][j] << ',' << P[i][j] << ',' << ux[i][j] << ',' << uy[i][j] << "\n";
		}
	}

	std::string local_data = buffer.str();
	int local_size = local_data.size();

	int offset = 0;
	MPI_Exscan(&local_size, &offset, 1, MPI_INT, MPI_SUM, comm);
	MPI_File fh;
	MPI_File_open(comm, name.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	MPI_File_write_at_all(fh, offset, local_data.c_str(), local_size, MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);
}


int main(int argc, char* argv[]) {
	// Godunov_Kolgan_solver_2D();

	MPI_Status Status;
	MPI_Request Request;
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// global boundary conditions
	int bound_left = 1, bound_right = 1, bound_down = 0, bound_up = 0;
	int minmod_type = 1;

	InitialState IS;
	int px, py;
	mpi_get_px_py(size, IS.nx, IS.ny, &px, &py);
	std::cout << "px " << px << "; py " << py << std::endl;

	int neighbors[4]; // 0 - left, 1 - up, 2 - right, 3 - down
	// -1 - left, -2 - up, -3 - right, -4 - down
	if (rank == 0) {
		neighbors[0] = -1;
		neighbors[1] = py == 1 ? -2 : rank + px;
		neighbors[2] = px == 1 ? -3 : rank + 1;
		neighbors[3] = -4;
	}
	else if (rank == px - 1) {
		neighbors[0] = px == 1 ? -1 : rank - 1;
		neighbors[1] = py == 1 ? -2 : rank + px;
		neighbors[2] = -3;
		neighbors[3] = -4;
	}
	else if (rank == px * (py - 1)) {
		neighbors[0] = -1;
		neighbors[1] = -2;
		neighbors[2] = px == 1 ? -3 : rank + 1;
		neighbors[3] = py == 1 ? -4 : rank - px;
	}
	else if (rank == px * py - 1) {
		neighbors[0] = px == 1 ? -1 : rank - 1;
		neighbors[1] = -2;
		neighbors[2] = -3;
		neighbors[3] = py == 1 ? -4 : rank - px;
	}
	else if (rank < px) {
		neighbors[0] = rank - 1;
		neighbors[1] = rank + px;
		neighbors[2] = rank + 1;
		neighbors[3] = -4;
	}
	else if (rank % px == 0) {
		neighbors[0] = -1;
		neighbors[1] = rank + px;
		neighbors[2] = rank + 1;
		neighbors[3] = rank - px;
	}
	else if (rank % px == px - 1) {
		neighbors[0] = rank - 1;
		neighbors[1] = rank + px;
		neighbors[2] = -3;
		neighbors[3] = rank - px;
	}
	else if (rank + px > px * py - 1) {
		neighbors[0] = rank - 1;
		neighbors[1] = -2;
		neighbors[2] = rank + 1;
		neighbors[3] = rank - px;
	}
	else {
		neighbors[0] = rank - 1;
		neighbors[1] = rank + px;
		neighbors[2] = rank + 1;
		neighbors[3] = rank - px;
	}


	int start_x, stop_x;
	int start_y, stop_y;
	int rank_x, rank_y;
	rank_x = rank % px;
	rank_y = rank / py;
	start_x = IS.nx / px * rank_x;
	stop_x = start_x + IS.nx / px;
	start_y = IS.ny / py * rank_y;
	stop_y = start_y + IS.ny / py;
	if (IS.nx % px) {
		if (rank_x >= px - IS.nx % px) {
			stop_x += IS.nx % px - (px - 1 - rank_x);
			start_x += IS.nx % px - 1 - (px - 1 - rank_x);
		}
	}
	if (IS.ny % py) {
		if (rank_y >= py - IS.ny % py) {
			stop_y += IS.ny % py - (py - 1 - rank_y);
			start_y += IS.ny % py - 1 - (py - 1 - rank_y);
		}
	}
	IS.nx = stop_x - start_x;
	IS.ny = stop_y - start_y;
	IS.x_start = start_x * IS.hx;
	IS.x_end = stop_x * IS.hx;
	IS.y_start = start_y * IS.hy;
	IS.y_end = stop_y * IS.hy;
	IS.Nx = IS.nx + 2 * IS.fict;
	IS.Ny = IS.ny + 2 * IS.fict;


	vec1d xc(IS.Nx);
	vec1d x(IS.Nx + 1);
	vec1d yc(IS.Ny);
	vec1d y(IS.Ny + 1);

	vec2d P(IS.Ny, vec1d(IS.Nx)), ux(IS.Ny, vec1d(IS.Nx)), uy(IS.Ny, vec1d(IS.Nx)), rho(IS.Ny, vec1d(IS.Nx));
	vec2d mass(IS.Ny, vec1d(IS.Nx)), Imp_x(IS.Ny, vec1d(IS.Nx)), Imp_y(IS.Ny, vec1d(IS.Nx)), rhoe(IS.Ny, vec1d(IS.Nx));

	int iter = 0, max_Nt = 101;
	double t_start = 0.0;
	double tau;
	std::string filename = "Godunov_2D";

	int iterwrite = 100;

	grid(IS, x, y, xc, yc);
	initial_state(IS.Nx, IS.Ny, IS.gamma, xc, yc, P, ux, uy, rho, mass, Imp_x, Imp_y, rhoe);

	// writeCSV(filename, iter, xc, yc, ux, uy, P, rho, t_start, IS.Nx, IS.Ny, IS.fict);
	writeCSV_p(filename, IS, iter, xc, yc, ux, uy, P, rho, t_start, rank, MPI_COMM_WORLD);

	while (t_start < IS.t_end && iter < max_Nt) {
		Godunov_Kolgan_solver_2D_one_step(IS, &t_start, x, y, P, ux, uy, rho, mass, Imp_x, Imp_y, rhoe, 
			neighbors, bound_left, bound_up, bound_right, bound_down, minmod_type);
		if (iter % iterwrite == 0) {
			// writeCSV(filename, iter, xc, yc, ux, uy, P, rho, t_start, IS.Nx, IS.Ny, IS.fict);
			writeCSV_p(filename, IS, iter, xc, yc, ux, uy, P, rho, t_start, rank, MPI_COMM_WORLD);
		}
		iter += 1;
	}
	//writeCSV(filename, iter, xc, yc, ux, uy, P, rho, t_start, IS.Nx, IS.Ny, IS.fict);
	writeCSV_p(filename, IS, iter, xc, yc, ux, uy, P, rho, t_start, rank, MPI_COMM_WORLD);

	if (t_start >= IS.t_end) {
		std::cout << "Solve stopped by time\n";
		std::cout << "On iter = " << iter << std::endl;
	}
	else {
		std::cout << "Solve stopped by iterations number";
	}
	return 0;
}