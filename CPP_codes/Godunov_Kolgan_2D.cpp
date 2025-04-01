#include <iostream>
#include <fstream>
#include <string> 
#include <vector>
#include <map>
#include <cmath>
#include <stdlib.h>

typedef std::vector<double> vec1d;
typedef std::vector<std::vector<double>> vec2d;


struct InitialState {
	int nx = 5;
	int ny = 5;
	int fict = 1;
	int Nx = nx + 2 * fict;
	int Ny = ny + 2 * fict;

	double x_start = 0;
	double x_end = 1;
	double y_start = 0;
	double y_end = 1;
	double Lx = x_end - x_start;
	double Ly = y_end - y_start;
	double hx = Lx / nx;
	double hy = Ly / ny;

	double t_end = 0.25;
	double time = 0;
	double tau;

	double CFL = 0.05;
	double gamma = 1.4;
	double Q = 2;

	double u_left, u_right, P_left, P_right, rho_left, rho_right;
};


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
	vx1 = 0.75;
	vy1 = 0.;
	rho1 = 1.0;

	p2 = 1.0;
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
			/*

			double PQ = pow(PL / PR, (gamma - 1.0) / (2.0 * gamma));
			UM = (PQ * UL / CL + UR / CR + 2.0 / (gamma - 1.0) * (PQ - 1.0f) / (PQ / CL + 1.0f / CR));
			double PTL = 1.0 + (gamma - 1.0) / 2.0 * (UL - UM) / CL;
			double PTR = 1.0 + (gamma - 1.0) / 2.0 * (UM - UR) / CR;
			PM = 0.5 * (PL * pow(PTL, 2.0 * gamma / (gamma - 1.0)) + PR * pow(PTR, 2.0 * gamma / (gamma - 1.0)));
			*/
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
void Rieman_solver(double gamma, double Quser, double PL, double DL, double UL, double PR, double DR, double UR, double& UM, double& PM) {
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
	Rieman_solver(gamma, Quser, PL, DL, UL, PR, DR, UR, UM, PM);

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


void Godunov_Kolgan_solver_2D() {
	InitialState IS;

	vec1d xc(IS.Nx);
	vec1d x(IS.Nx + 1);
	vec1d yc(IS.Ny);
	vec1d y(IS.Ny + 1);

	vec2d P(IS.Ny, vec1d(IS.Nx)), ux(IS.Ny, vec1d(IS.Nx)), uy(IS.Ny, vec1d(IS.Nx)), rho(IS.Ny, vec1d(IS.Nx));
	vec2d mass(IS.Ny, vec1d(IS.Nx)), Imp_x(IS.Ny, vec1d(IS.Nx)), Imp_y(IS.Ny, vec1d(IS.Nx)), rhoe(IS.Ny, vec1d(IS.Nx));
	vec2d new_mass(IS.Ny, vec1d(IS.Nx)), new_Imp_x(IS.Ny, vec1d(IS.Nx)), new_Imp_y(IS.Ny, vec1d(IS.Nx)), new_rhoe(IS.Ny, vec1d(IS.Nx));


	double mb, impxb, impyb, eb, pb, vxb, vyb, rhob;
	double FmL, FimpxL, FimpyL, FeL, FmR, FimpxR, FimpyR, FeR;

	double pl, vl, rhol, pr, vr, rhor;
	double massl, impxl, impyl, rhoel, massr, impxr, impyr, rhoer;
	double dm, dimp, de;

	int iter = 0, max_Nt = 101;
	double t_start = 0.0;
	double tau;
	std::string filename = "Godunov_2D";

	int iterwrite = 100;
	int minmod_type = 1, bound_left = 1, bound_right = 1, bound_down = 0, bound_up = 0;

	grid(IS, x, y, xc, yc);
	initial_state(IS.Nx, IS.Ny, IS.gamma, xc, yc, P, ux, uy, rho, mass, Imp_x, Imp_y, rhoe);

	writeCSV(filename, iter, xc, yc, ux, uy, P, rho, t_start, IS.Nx, IS.Ny, IS.fict);

	while (t_start < IS.t_end && iter < max_Nt) {
		tau = 1.e6;
		double c, temp;
		for (int i = IS.fict; i < IS.Ny - IS.fict; ++i) {
			for (int j = IS.fict; j < IS.Nx - IS.fict; ++j) {
				c = sqrt(IS.gamma * P[i][j] / rho[i][j]);
				temp = std::min(abs(IS.CFL * (x[i + 1] - x[i]) / (abs(ux[i][j]) + c)), abs(IS.CFL * (y[j + 1] - y[j]) / (abs(uy[i][j]) + c)));
				if (temp < tau) {
					tau = temp;
				}
			}
		}

		for (int i = 0; i < IS.Ny; i++) {
			for (int j = 0; j < IS.Nx; j++) {
				if (i == 0) {
					Boundary_x(P[0][j], ux[0][j], uy[0][j], rho[0][j], pb, vxb, vyb, rhob, bound_left);
					noncons_to_cons(IS.gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);
					massl = mb;
					impxl = impxb;
					impyl = impyb;
					rhoel = eb;

					massr = mass[0][j] - 0.5 * minmod(mass[0][j] - mb, mass[1][j] - mass[0][j], minmod_type);
					impxr = Imp_x[0][j] - 0.5 * minmod(Imp_x[0][j] - impxb, Imp_x[1][j] - Imp_x[0][j], minmod_type);
					impyr = Imp_y[0][j] - 0.5 * minmod(Imp_y[0][j] - impyb, Imp_y[1][j] - Imp_y[0][j], minmod_type);
					rhoer = rhoe[0][j] - 0.5 * minmod(rhoe[0][j] - eb, rhoe[1][j] - rhoe[0][j], minmod_type);
				}

				else if (i == 1) {
					Boundary_x(P[0][j], ux[0][j], uy[0][j], rho[0][j], pb, vxb, vyb, rhob, bound_left);
					noncons_to_cons(IS.gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					massl = mass[0][j] + 0.5 * minmod(mass[0][j] - mb, mass[1][j] - mass[0][j], minmod_type);
					impxl = Imp_x[0][j] + 0.5 * minmod(Imp_x[0][j] - impxb, Imp_x[1][j] - Imp_x[0][j], minmod_type);
					impyl = Imp_y[0][j] + 0.5 * minmod(Imp_y[0][j] - impyb, Imp_y[1][j] - Imp_y[0][j], minmod_type);
					rhoel = rhoe[0][j] + 0.5 * minmod(rhoe[0][j] - eb, rhoe[1][j] - rhoe[0][j], minmod_type);

					massr = mass[1][j] - 0.5 * minmod(mass[i][j] - mass[i - 1][j], mass[i + 1][j] - mass[i][j], minmod_type);
					impxr = Imp_x[1][j] - 0.5 * minmod(Imp_x[i][j] - Imp_x[i - 1][j], Imp_x[i + 1][j] - Imp_x[i][j], minmod_type);
					impyr = Imp_y[1][j] - 0.5 * minmod(Imp_y[i][j] - Imp_y[i - 1][j], Imp_y[i + 1][j] - Imp_y[i][j], minmod_type);
					rhoer = rhoe[1][j] - 0.5 * minmod(rhoe[i][j] - rhoe[i - 1][j], rhoe[i + 1][j] - rhoe[i][j], minmod_type);
				}

				else if (i == IS.Nx - 1) {
					Boundary_x(P[IS.Nx - 1][j], ux[IS.Nx - 1][j], uy[IS.Nx - 1][j], rho[IS.Nx - 1][j], pb, vxb, vyb, rhob, bound_right);
					noncons_to_cons(IS.gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					massl = mass[i - 1][j] + 0.5 * minmod(mass[i - 1][j] - mass[i - 2][j], mass[i][j] - mass[i - 1][j], minmod_type);
					impxl = Imp_x[i - 1][j] + 0.5 * minmod(Imp_x[i - 1][j] - Imp_x[i - 2][j], Imp_x[i][j] - Imp_x[i - 1][j], minmod_type);
					impyl = Imp_y[i - 1][j] + 0.5 * minmod(Imp_y[i - 1][j] - Imp_y[i - 2][j], Imp_y[i][j] - Imp_y[i - 1][j], minmod_type);
					rhoel = rhoe[i - 1][j] + 0.5 * minmod(rhoe[i - 1][j] - rhoe[i - 2][j], rhoe[i][j] - rhoe[i - 1][j], minmod_type);

					massr = mass[i][j] - 0.5 * minmod(mass[i][j] - mass[i - 1][j], mb - mass[i][j], minmod_type);
					impxr = Imp_x[i][j] - 0.5 * minmod(Imp_x[i][j] - Imp_x[i - 1][j], impxb - Imp_x[i][j], minmod_type);
					impyr = Imp_y[i][j] - 0.5 * minmod(Imp_y[i][j] - Imp_y[i - 1][j], impyb - Imp_y[i][j], minmod_type);
					rhoer = rhoe[i][j] - 0.5 * minmod(rhoe[i][j] - rhoe[i - 1][j], eb - rhoe[i][j], minmod_type);
				}

				else {
					massl = mass[i - 1][j] + 0.5 * minmod(mass[i - 1][j] - mass[i - 2][j], mass[i][j] - mass[i - 1][j], minmod_type);
					impxl = Imp_x[i - 1][j] + 0.5 * minmod(Imp_x[i - 1][j] - Imp_x[i - 2][j], Imp_x[i][j] - Imp_x[i - 1][j], minmod_type);
					impyl = Imp_y[i - 1][j] + 0.5 * minmod(Imp_y[i - 1][j] - Imp_y[i - 2][j], Imp_y[i][j] - Imp_y[i - 1][j], minmod_type);
					rhoel = rhoe[i - 1][j] + 0.5 * minmod(rhoe[i - 1][j] - rhoe[i - 2][j], rhoe[i][j] - rhoe[i - 1][j], minmod_type);

					massr = mass[i][j] - 0.5 * minmod(mass[i][j] - mass[i - 1][j], mass[i + 1][j] - mass[i][j], minmod_type);
					impxr = Imp_x[i][j] - 0.5 * minmod(Imp_x[i][j] - Imp_x[i - 1][j], Imp_x[i + 1][j] - Imp_x[i][j], minmod_type);
					impyr = Imp_y[i][j] - 0.5 * minmod(Imp_y[i][j] - Imp_y[i - 1][j], Imp_y[i + 1][j] - Imp_y[i][j], minmod_type);
					rhoer = rhoe[i][j] - 0.5 * minmod(rhoe[i][j] - rhoe[i - 1][j], rhoe[i + 1][j] - rhoe[i][j], minmod_type);
				}
				Godunov_method_x(IS.gamma, IS.Q, massl, impxl, impyl, rhoel, massr, impxr, impyr, rhoer, FmL, FimpxL, FimpyL, FeL);

				//		����� ����� ������ ����� ������
				if (i == IS.Nx - 1) {
					Boundary_x(P[IS.Nx - 1][j], ux[IS.Nx - 1][j], uy[IS.Nx - 1][j], rho[IS.Nx - 1][j], pb, vxb, vyb, rhob, bound_right);
					noncons_to_cons(IS.gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					massl = mass[i][j] + 0.5 * minmod(mass[i][j] - mass[i - 1][j], mb - mass[i][j], minmod_type);
					impxl = Imp_x[i][j] + 0.5 * minmod(Imp_x[i][j] - Imp_x[i - 1][j], impxb - Imp_x[i][j], minmod_type);
					impyl = Imp_y[i][j] + 0.5 * minmod(Imp_y[i][j] - Imp_y[i - 1][j], impyb - Imp_y[i][j], minmod_type);
					rhoel = rhoe[i][j] + 0.5 * minmod(rhoe[i][j] - rhoe[i - 1][j], eb - rhoe[i][j], minmod_type);

					massr = mb; //m[IS.Nx - 1][j];
					impxr = impxb; //impx[IS.Nx - 1][j];
					impyr = impyb; //impy[IS.Nx - 1][j];
					rhoer = eb; //rhoe[IS.Nx - 1][j];

				}
				else if (i == IS.Nx - 2) {
					Boundary_x(P[IS.Nx - 1][j], ux[IS.Nx - 1][j], uy[IS.Nx - 1][j], rho[IS.Nx - 1][j], pb, vxb, vyb, rhob, bound_right);
					noncons_to_cons(IS.gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					massl = mass[i][j] + 0.5 * minmod(mass[i][j] - mass[i - 1][j], mass[i + 1][j] - mass[i][j], minmod_type);
					impxl = Imp_x[i][j] + 0.5 * minmod(Imp_x[i][j] - Imp_x[i - 1][j], Imp_x[i + 1][j] - Imp_x[i][j], minmod_type);
					impyl = Imp_y[i][j] + 0.5 * minmod(Imp_y[i][j] - Imp_y[i - 1][j], Imp_y[i + 1][j] - Imp_y[i][j], minmod_type);
					rhoel = rhoe[i][j] + 0.5 * minmod(rhoe[i][j] - rhoe[i - 1][j], rhoe[i + 1][j] - rhoe[i][j], minmod_type);

					massr = mass[i + 1][j] - 0.5 * minmod(mass[i + 1][j] - mass[i][j], mb - mass[i + 1][j], minmod_type);
					impxr = Imp_x[i + 1][j] - 0.5 * minmod(Imp_x[i + 1][j] - Imp_x[i][j], impxb - Imp_x[i + 1][j], minmod_type);
					impyr = Imp_y[i + 1][j] - 0.5 * minmod(Imp_y[i + 1][j] - Imp_y[i][j], impyb - Imp_y[i + 1][j], minmod_type);
					rhoer = rhoe[i + 1][j] - 0.5 * minmod(rhoe[i + 1][j] - rhoe[i][j], eb - rhoe[i + 1][j], minmod_type);
				}
				else if (i == 0) {
					Boundary_x(P[0][j], ux[0][j], uy[0][j], rho[0][j], pb, vxb, vyb, rhob, bound_left);
					noncons_to_cons(IS.gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					massl = mass[i][j] + 0.5 * minmod(mass[i][j] - mb, mass[i + 1][j] - mass[i][j], minmod_type);
					impxl = Imp_x[i][j] + 0.5 * minmod(Imp_x[i][j] - impxb, Imp_x[i + 1][j] - Imp_x[i][j], minmod_type);
					impyl = Imp_y[i][j] + 0.5 * minmod(Imp_y[i][j] - impyb, Imp_y[i + 1][j] - Imp_y[i][j], minmod_type);
					rhoel = rhoe[i][j] + 0.5 * minmod(rhoe[i][j] - eb, rhoe[i + 1][j] - rhoe[i][j], minmod_type);

					massr = mass[i + 1][j] - 0.5 * minmod(mass[i + 1][j] - mass[i][j], mass[i + 2][j] - mass[i + 1][j], minmod_type);
					impxr = Imp_x[i + 1][j] - 0.5 * minmod(Imp_x[i + 1][j] - Imp_x[i][j], Imp_x[i + 2][j] - Imp_x[i + 1][j], minmod_type);
					impyr = Imp_y[i + 1][j] - 0.5 * minmod(Imp_y[i + 1][j] - Imp_y[i][j], Imp_y[i + 2][j] - Imp_y[i + 1][j], minmod_type);
					rhoer = rhoe[i + 1][j] - 0.5 * minmod(rhoe[i + 1][j] - rhoe[i][j], rhoe[i + 2][j] - rhoe[i + 1][j], minmod_type);
				}
				else {
					massl = mass[i][j] + 0.5 * minmod(mass[i][j] - mass[i - 1][j], mass[i + 1][j] - mass[i][j], minmod_type);
					impxl = Imp_x[i][j] + 0.5 * minmod(Imp_x[i][j] - Imp_x[i - 1][j], Imp_x[i + 1][j] - Imp_x[i][j], minmod_type);
					impyl = Imp_y[i][j] + 0.5 * minmod(Imp_y[i][j] - Imp_y[i - 1][j], Imp_y[i + 1][j] - Imp_y[i][j], minmod_type);
					rhoel = rhoe[i][j] + 0.5 * minmod(rhoe[i][j] - rhoe[i - 1][j], rhoe[i + 1][j] - rhoe[i][j], minmod_type);

					massr = mass[i + 1][j] - 0.5 * minmod(mass[i + 1][j] - mass[i][j], mass[i + 2][j] - mass[i + 1][j], minmod_type);
					impxr = Imp_x[i + 1][j] - 0.5 * minmod(Imp_x[i + 1][j] - Imp_x[i][j], Imp_x[i + 2][j] - Imp_x[i + 1][j], minmod_type);
					impyr = Imp_y[i + 1][j] - 0.5 * minmod(Imp_y[i + 1][j] - Imp_y[i][j], Imp_y[i + 2][j] - Imp_y[i + 1][j], minmod_type);
					rhoer = rhoe[i + 1][j] - 0.5 * minmod(rhoe[i + 1][j] - rhoe[i][j], rhoe[i + 2][j] - rhoe[i + 1][j], minmod_type);
				}
				Godunov_method_x(IS.gamma, IS.Q, massl, impxl, impyl, rhoel, massr, impxr, impyr, rhoer, FmR, FimpxR, FimpyR, FeR);

				//std::cout << FmL << " " << FmR << std::endl;

				new_mass[i][j] = mass[i][j] - tau * (FmR - FmL) / IS.hx;
				new_Imp_x[i][j] = Imp_x[i][j] - tau * (FimpxR - FimpxL) / IS.hx;
				new_Imp_y[i][j] = Imp_y[i][j] - tau * (FimpyR - FimpyL) / IS.hx;
				new_rhoe[i][j] = rhoe[i][j] - tau * (FeR - FeL) / IS.hx;

				if (j == 0) {
					Boundary_x(P[i][0], ux[i][0], uy[i][0], rho[i][0], pb, vxb, vyb, rhob, bound_down);
					noncons_to_cons(IS.gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);
					massl = mb;
					impxl = impxb;
					impyl = impyb;
					rhoel = eb;

					massr = mass[i][0] - 0.5 * minmod(mass[i][0] - mb, mass[i][1] - mass[i][0], minmod_type);
					impxr = Imp_x[i][0] - 0.5 * minmod(Imp_x[i][0] - impxb, Imp_x[i][1] - Imp_x[i][0], minmod_type);
					impyr = Imp_y[i][0] - 0.5 * minmod(Imp_y[i][0] - impyb, Imp_y[i][1] - Imp_y[i][0], minmod_type);
					rhoer = rhoe[i][0] - 0.5 * minmod(rhoe[i][0] - eb, rhoe[i][1] - rhoe[i][0], minmod_type);
				}

				else if (j == 1) {
					Boundary_y(P[i][0], ux[i][0], uy[i][0], rho[i][0], pb, vxb, vyb, rhob, bound_down);
					noncons_to_cons(IS.gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					massr = mass[i][0] + 0.5 * minmod(mass[i][0] - mb, mass[i][1] - mass[i][0], minmod_type);
					impxr = Imp_x[i][0] + 0.5 * minmod(Imp_x[i][0] - impxb, Imp_x[i][1] - Imp_x[i][0], minmod_type);
					impyr = Imp_y[i][0] + 0.5 * minmod(Imp_y[i][0] - impyb, Imp_y[i][1] - Imp_y[i][0], minmod_type);
					rhoer = rhoe[i][0] + 0.5 * minmod(rhoe[i][0] - eb, rhoe[i][1] - rhoe[i][0], minmod_type);

					massr = mass[i][1] - 0.5 * minmod(mass[i][j] - mass[i][j - 1], mass[i][j + 1] - mass[i][j], minmod_type);
					impxr = Imp_x[i][1] - 0.5 * minmod(Imp_x[i][j] - Imp_x[i][j - 1], Imp_x[i][j + 1] - Imp_x[i][j], minmod_type);
					impyr = Imp_y[i][1] - 0.5 * minmod(Imp_y[i][j] - Imp_y[i][j - 1], Imp_y[i][j + 1] - Imp_y[i][j], minmod_type);
					rhoer = rhoe[i][1] - 0.5 * minmod(rhoe[i][j] - rhoe[i][j - 1], rhoe[i][j + 1] - rhoe[i][j], minmod_type);
				}

				else if (j == IS.Ny - 1) {
					Boundary_y(P[i][IS.Ny - 1], ux[i][IS.Ny - 1], uy[i][IS.Ny - 1], rho[i][IS.Ny - 1], pb, vxb, vyb, rhob, bound_up);
					noncons_to_cons(IS.gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					massl = mass[i][j - 1] + 0.5 * minmod(mass[i][j - 1] - mass[i][j - 2], mass[i][j] - mass[i][j - 1], minmod_type);
					impxl = Imp_x[i][j - 1] + 0.5 * minmod(Imp_x[i][j - 1] - Imp_x[i][j - 2], Imp_x[i][j] - Imp_x[i][j - 1], minmod_type);
					impyl = Imp_y[i][j - 1] + 0.5 * minmod(Imp_y[i][j - 1] - Imp_y[i][j - 2], Imp_y[i][j] - Imp_y[i][j - 1], minmod_type);
					rhoel = rhoe[i][j - 1] + 0.5 * minmod(rhoe[i][j - 1] - rhoe[i][j - 2], rhoe[i][j] - rhoe[i][j - 1], minmod_type);

					massr = mass[i][j] - 0.5 * minmod(mass[i][j] - mass[i][j - 1], mb - mass[i][j], minmod_type);
					impxr = Imp_x[i][j] - 0.5 * minmod(Imp_x[i][j] - Imp_x[i][j - 1], impxb - Imp_x[i][j], minmod_type);
					impyr = Imp_y[i][j] - 0.5 * minmod(Imp_y[i][j] - Imp_y[i][j - 1], impyb - Imp_y[i][j], minmod_type);
					rhoer = rhoe[i][j] - 0.5 * minmod(rhoe[i][j] - rhoe[i][j - 1], eb - rhoe[i][j], minmod_type);
				}

				else {
					massl = mass[i][j - 1] + 0.5 * minmod(mass[i][j - 1] - mass[i][j - 2], mass[i][j] - mass[i][j - 1], minmod_type);
					impxl = Imp_x[i][j - 1] + 0.5 * minmod(Imp_x[i][j - 1] - Imp_x[i][j - 2], Imp_x[i][j] - Imp_x[i][j - 1], minmod_type);
					impyl = Imp_y[i][j - 1] + 0.5 * minmod(Imp_y[i][j - 1] - Imp_y[i][j - 2], Imp_y[i][j] - Imp_y[i][j - 1], minmod_type);
					rhoel = rhoe[i][j - 1] + 0.5 * minmod(rhoe[i][j - 1] - rhoe[i][j - 2], rhoe[i][j] - rhoe[i][j - 1], minmod_type);

					massr = mass[i][j] - 0.5 * minmod(mass[i][j] - mass[i][j - 1], mass[i][j + 1] - mass[i][j], minmod_type);
					impxr = Imp_x[i][j] - 0.5 * minmod(Imp_x[i][j] - Imp_x[i][j - 1], Imp_x[i][j + 1] - Imp_x[i][j], minmod_type);
					impyr = Imp_y[i][j] - 0.5 * minmod(Imp_y[i][j] - Imp_y[i][j - 1], Imp_y[i][j + 1] - Imp_y[i][j], minmod_type);
					rhoer = rhoe[i][j] - 0.5 * minmod(rhoe[i][j] - rhoe[i][j - 1], rhoe[i][j + 1] - rhoe[i][j], minmod_type);
				}
				Godunov_method_y(IS.gamma, IS.Q, massl, impxl, impyl, rhoel, massr, impxr, impyr, rhoer, FmL, FimpxL, FimpyL, FeL);

				if (j == IS.Ny - 1) {
					Boundary_y(P[i][IS.Ny - 1], ux[i][IS.Ny - 1], uy[i][IS.Ny - 1], rho[i][IS.Ny - 1], pb, vxb, vyb, rhob, bound_up);
					noncons_to_cons(IS.gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					massl = mass[i][j] + 0.5 * minmod(mass[i][j] - mass[i][j - 1], mb - mass[i][j], minmod_type);
					impxl = Imp_x[i][j] + 0.5 * minmod(Imp_x[i][j] - Imp_x[i][j - 1], impxb - Imp_x[i][j], minmod_type);
					impyl = Imp_y[i][j] + 0.5 * minmod(Imp_y[i][j] - Imp_y[i][j - 1], impyb - Imp_y[i][j], minmod_type);
					rhoel = rhoe[i][j] + 0.5 * minmod(rhoe[i][j] - rhoe[i][j - 1], eb - rhoe[i][j], minmod_type);

					massr = mb;
					impxr = impxb;
					impyr = impyb;
					rhoer = eb;
				}
				else if (j == IS.Ny - 2) {
					Boundary_y(P[i][IS.Ny - 1], ux[i][IS.Ny - 1], uy[i][IS.Ny - 1], rho[i][IS.Ny - 1], pb, vxb, vyb, rhob, bound_up);
					noncons_to_cons(IS.gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					massl = mass[i][j] + 0.5 * minmod(mass[i][j] - mass[i][j - 1], mass[i][j + 1] - mass[i][j], minmod_type);
					impxl = Imp_x[i][j] + 0.5 * minmod(Imp_x[i][j] - Imp_x[i][j - 1], Imp_x[i][j + 1] - Imp_x[i][j], minmod_type);
					impyl = Imp_y[i][j] + 0.5 * minmod(Imp_y[i][j] - Imp_y[i][j - 1], Imp_y[i][j + 1] - Imp_y[i][j], minmod_type);
					rhoel = rhoe[i][j] + 0.5 * minmod(rhoe[i][j] - rhoe[i][j - 1], rhoe[i][j + 1] - rhoe[i][j], minmod_type);

					massr = mass[i][j + 1] - 0.5 * minmod(mass[i][j + 1] - mass[i][j], mb - mass[i][j + 1], minmod_type);
					impxr = Imp_x[i][j + 1] - 0.5 * minmod(Imp_x[i][j + 1] - Imp_x[i][j], impxb - Imp_x[i][j + 1], minmod_type);
					impyr = Imp_y[i][j + 1] - 0.5 * minmod(Imp_y[i][j + 1] - Imp_y[i][j], impyb - Imp_y[i][j + 1], minmod_type);
					rhoer = rhoe[i][j + 1] - 0.5 * minmod(rhoe[i][j + 1] - rhoe[i][j], eb - rhoe[i][j + 1], minmod_type);
				}
				else if (j == 0) {
					Boundary_y(P[i][0], ux[i][0], uy[i][0], rho[i][0], pb, vxb, vyb, rhob, bound_down);
					noncons_to_cons(IS.gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					massl = mass[i][j] + 0.5 * minmod(mass[i][j] - mb, mass[i][j + 1] - mass[i][j], minmod_type);
					impxl = Imp_x[i][j] + 0.5 * minmod(Imp_x[i][j] - impxb, Imp_x[i][j + 1] - Imp_x[i][j], minmod_type);
					impyl = Imp_y[i][j] + 0.5 * minmod(Imp_y[i][j] - impyb, Imp_y[i][j + 1] - Imp_y[i][j], minmod_type);
					rhoel = rhoe[i][j] + 0.5 * minmod(rhoe[i][j] - eb, rhoe[i][j + 1] - rhoe[i][j], minmod_type);

					massr = mass[i][j + 1] - 0.5 * minmod(mass[i][j + 1] - mass[i][j], mass[i][j + 2] - mass[i][j + 1], minmod_type);
					impxr = Imp_x[i][j + 1] - 0.5 * minmod(Imp_x[i][j + 1] - Imp_x[i][j], Imp_x[i][j + 2] - Imp_x[i][j + 1], minmod_type);
					impyr = Imp_y[i][j + 1] - 0.5 * minmod(Imp_y[i][j + 1] - Imp_y[i][j], Imp_y[i][j + 2] - Imp_y[i][j + 1], minmod_type);
					rhoer = rhoe[i][j + 1] - 0.5 * minmod(rhoe[i][j + 1] - rhoe[i][j], rhoe[i][j + 2] - rhoe[i][j + 1], minmod_type);
				}
				else {
					massl = mass[i][j] + 0.5 * minmod(mass[i][j] - mass[i][j - 1], mass[i][j + 1] - mass[i][j], minmod_type);
					impxl = Imp_x[i][j] + 0.5 * minmod(Imp_x[i][j] - Imp_x[i][j - 1], Imp_x[i][j + 1] - Imp_x[i][j], minmod_type);
					impyl = Imp_y[i][j] + 0.5 * minmod(Imp_y[i][j] - Imp_y[i][j - 1], Imp_y[i][j + 1] - Imp_y[i][j], minmod_type);
					rhoel = rhoe[i][j] + 0.5 * minmod(rhoe[i][j] - rhoe[i][j - 1], rhoe[i][j + 1] - rhoe[i][j], minmod_type);

					massr = mass[i][j + 1] - 0.5 * minmod(mass[i][j + 1] - mass[i][j], mass[i][j + 2] - mass[i][j + 1], minmod_type);
					impxr = Imp_x[i][j + 1] - 0.5 * minmod(Imp_x[i][j + 1] - Imp_x[i][j], Imp_x[i][j + 2] - Imp_x[i][j + 1], minmod_type);
					impyr = Imp_y[i][j + 1] - 0.5 * minmod(Imp_y[i][j + 1] - Imp_y[i][j], Imp_y[i][j + 2] - Imp_y[i][j + 1], minmod_type);
					rhoer = rhoe[i][j + 1] - 0.5 * minmod(rhoe[i][j + 1] - rhoe[i][j], rhoe[i][j + 2] - rhoe[i][j + 1], minmod_type);
				}
				Godunov_method_y(IS.gamma, IS.Q, massl, impxl, impyl, rhoel, massr, impxr, impyr, rhoer, FmR, FimpxR, FimpyR, FeR);
				//std::cout << FmL << " " << FmR << std::endl;

				new_mass[i][j] = new_mass[i][j] - tau * (FmR - FmL) / IS.hy;
				new_Imp_x[i][j] = new_Imp_x[i][j] - tau * (FimpxR - FimpxL) / IS.hy;
				new_Imp_y[i][j] = new_Imp_y[i][j] - tau * (FimpyR - FimpyL) / IS.hy;
				new_rhoe[i][j] = new_rhoe[i][j] - tau * (FeR - FeL) / IS.hy;

			}
		}

		// update parameters
		for (int i = 0; i < IS.Ny; ++i) {
			for (int j = 0; j < IS.Nx; ++j) {
				mass[i][j] = new_mass[i][j];
				Imp_x[i][j] = new_Imp_x[i][j];
				Imp_y[i][j] = new_Imp_y[i][j];
				rhoe[i][j] = new_rhoe[i][j];

				cons_to_noncons(IS.gamma, P[i][j], ux[i][j], uy[i][j], rho[i][j], mass[i][j], Imp_x[i][j], Imp_y[i][j], rhoe[i][j]);
			}
		}
		//data_to_file(int Nx, int Ny, double time, vec x, vec y, mtrx p, mtrx vx, mtrx vy, mtrx rho)
		if (iter % iterwrite == 0) {
			writeCSV(filename, iter, xc, yc, ux, uy, P, rho, t_start, IS.Nx, IS.Ny, IS.fict);
		}
		iter += 1;
		t_start += tau;

	}
	writeCSV(filename, iter, xc, yc, ux, uy, P, rho, t_start, IS.Nx, IS.Ny, IS.fict);

	if (t_start >= IS.t_end) {
		std::cout << "Solve stopped by time\n";
		std::cout << "On iter = " << iter << std::endl;
	}
	else {
		std::cout << "Solve stopped by iterations number";
	}
	return;
}


int main() {
	Godunov_Kolgan_solver_2D();
	return 0;
}