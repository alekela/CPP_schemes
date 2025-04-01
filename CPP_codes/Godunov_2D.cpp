#include <iostream>
#include <fstream>
#include <string> 
#include <vector>
#include <map>
#include <cmath>
#include <stdlib.h>

typedef std::vector<double> vec;
typedef std::vector<std::vector<double>> mtrx;

struct Params
{
	double CFL;
	double gamma;
	double dx, dy;
	int scheme_type;
	double Lx, Ly;
	double T;
	int Nx, Ny;
	double Quser;	
	int fict;			

	Params()
	{
		this->CFL = 0.05;
		this->gamma = 1.4;
		this->Nx = 20;
		this->Lx = 1.0;
		this->Ly = 1.0;
		this->Ny = 20;
		this->T = 1.;
		this->Quser = 2.;
		this->dx = Lx / Nx;
		this->dy = Ly / Ny;
		this->scheme_type = 1;
		if (scheme_type == 0) this->fict = 1;
		else this->fict = 2;
	}

	Params(int Nx_, int Ny_, int type, double cfl_, double Lx_, double Ly_, double T_, double gamma_, double Quser_)
	{
		this->CFL = cfl_;
		this->gamma = gamma_;
		this->scheme_type = type;

		if (type == 0) this->fict = 1;
		else this->fict = 2;

		this->Nx = Nx_;// +2 * this->fict;
		this->Lx = Lx_;
		this->Ly = Ly_;
		this->Ny = Ny_;// +2 * this->fict;
		this->T = T_;
		this->Quser = Quser_;
		this->dx = Lx / Nx_;
		this->dy = Ly / Ny_;
	}
};

void data_to_file(int Nx, int Ny, int fict, double gamma, double time, vec x, vec y, mtrx p, mtrx vx, mtrx vy, mtrx rho) {
	char name_file[100];
	//sprintf_s(name_file, "out/csv/out_%f_.csv", time);
	sprintf(name_file, "out/csv/out_%f_.csv", time);
	std::ofstream outfile(name_file, std::ios::app);
	outfile << "x;y;p;vx;vy;r;e;\n";
	for (int i = fict; i < Nx - fict; ++i) {
		for (int j = fict; j < Ny - fict; ++j) {
			outfile << x[i] << ";" << y[j] << ";" << p[i][j] << ";" << vx[i][j] << ";" << vy[i][j] << ";" << rho[i][j] << ";" << p[i][j] / ((gamma - 1.0) * rho[i][j]) << "\n";
		}
	}
}

// from conservative
void convert_from_conservative(double gamma, double& p, double& vx, double& vy, double& rho, double& m, double& impx, double& impy, double& e)
{
	p = (gamma - 1.0) * (e - 0.5 * (pow(impx, 2.0) + pow(impy, 2.0)) * m);
	vx = impx / m;
	vy = impy / m;
	rho = m;
}
// to conservative
void convert_to_conservative(double gamma, double& p, double& vx, double& vy, double& rho, double& m, double& impx, double& impy, double& e)
{
	m = rho;
	impx = rho * vx;
	e = 0.5 * rho * (pow(vx, 2.0) + pow(vy, 2.0)) + p / (gamma - 1.0);
}


void boundary_cond_x(double p, double vx, double vy, double rho, double& pb, double& vxb, double& vyb, double& rhob, int mode)
{
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

void boundary_cond_y(double p, double vx, double vy, double rho, double& pb, double& vxb, double& vyb, double& rhob, int mode)
{
	// wall
	if (mode == 0)
	{
		pb = p;
		vxb = vx;
		vyb = -vy;
		rhob = rho;
	}
	// free flux, ���� ������� �� ���������
	else
	{
		pb = p;
		vxb = vx;
		vyb = vy;
		rhob = rho;
	}
}

void grid(Params* params, vec& x, vec& y, vec& xc, vec& yc)
{
	for (int i = 0; i < params->Nx + 1; ++i) {
		x[i] = 0.0 + (i - params->fict) * params->dx;
	}
	for (int i = 0; i < params->Nx; ++i) {
		xc[i] = 0.5 * (x[i] + x[i + 1]);
	}
	for (int i = 0; i < params->Ny + 1; ++i) {
		y[i] = 0.0 + (i - params->fict) * params->dy;
	}

	for (int i = 0; i < params->Ny; ++i) {
		yc[i] = 0.5 * (y[i] + y[i + 1]);
	}
}

void init_krujok(int Nx, int Ny, double gamma, vec xc, vec yc, mtrx& p, mtrx& vx, mtrx& vy, mtrx& rho, mtrx& m, mtrx& impx, mtrx& impy, mtrx& e) {
	double R = 0.1;
	double p1, vx1, vy1, rho1, p2, vx2, vy2, rho2, d;

	// toro test 1

	p1 = 1.0; // ����������� ��������
	vx1 = 0.75;
	vy1 = 0.;
	rho1 = 1.0; // ��������� �������

	p2 = 1.0; // ����������� ��������
	vx2 = 0.;
	vy2 = 0.;
	rho2 = 0.125; // ��������� �������

	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			d = (xc[i] - 0.5) * (xc[i] - 0.5) + (yc[j] - 0.5) * (yc[j] - 0.5);
			if (d < R * R) {

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
			convert_to_conservative(gamma, p[i][j], vx[i][j], vy[i][j], rho[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
		}
	}
}

double get_dt(Params* params, mtrx p, mtrx rho, mtrx vx, mtrx vy, vec x, vec y)
{
	double dt = 10.0e6;
	double c, dt_step;
	for (int i = params->fict; i < params->Nx - params->fict; ++i) {
		for (int j = params->fict; j < params->Ny - params->fict; ++j) {
			c = sqrt(params->gamma * p[i][j] / rho[i][j]);
			dt_step = std::min(abs(params->CFL * (x[i + 1] - x[i]) / (abs(vx[i][j]) + c)), abs(params->CFL * (y[j + 1] - y[j]) / (abs(vy[i][j]) + c)));
			if (dt_step < dt) {
				dt = dt_step;
			}
		}
	}
	return dt;
}

/*



*/

// minmod functions: 0 - Kolgan,1972; 1 - Kolgan,1975; 2 (or another...) - Osher,1984
double minmod(double a, double b, int func) // func - type of minmod we use
{
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
bool is_vacuum(double gamma, double CL, double UL, double CR, double UR)// C - sound's speed
{
	if (2.0 / (gamma - 1.0) * (CL + CR) < UR - UL) return true;
	return false;
}

// in F calculate U on contact, in FD calculate dU/dP on contact
void prefun(double gamma, double& F, double& FD, double P, double DK, double PK, double CK) // PK  =  PR or PL
{

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

void guess_p(double gamma, double Quser, double PL, double DL, double UL, double PR, double DR, double UR, double& P)
{
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
void starpu(double gamma, double Quser, double PL, double DL, double UL, double PR, double DR, double UR, double& UM, double& PM)
{
	int NRITER = 30;
	double CHANGE = 1e6, FL, FLD, FR, FRD, POLD, TOLPRE = 1.0e-6, UDIFF;
	double CL = sqrt(gamma * PL / DL);
	double CR = sqrt(gamma * PR / DR);

	//std::cout << PL << " " << PR << std::endl;

	if (is_vacuum(gamma, CL, UL, CR, UR)) {
		std::cout << "Vacuum is generated" << std::endl << "Program stopped" << std::endl;
		exit(228);
	}
	guess_p(gamma, Quser, PL, DL, UL, PR, DR, UR, POLD);		// initial condition for pressure

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
void sample(double gamma, double PL, double DL, double UL, double PR, double DR, double UR,
	double UM, double PM, double S, double& D, double& U, double& P)
{
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

	if (S <= UM) // point is left from contact
	{
		if (PM <= PL) //left rarefaction wave
		{
			SHL = UL - CL;
			if (S <= SHL) //point is left from rarefaction wave - left parameters
			{
				D = DL;
				U = UL;
				P = PL;
			}
			else
			{
				CML = CL * pow((PM / PL), (gamma - 1.0) / (2.0 * gamma));
				STL = UM - CML;
				if (S > STL) // point is a state with *
				{
					D = DL * pow((PM / PL), (1. / gamma));
					U = UM;
					P = PM;
				}
				else //point is in a left rarefaction wave
				{
					U = 2.0 / (gamma + 1.0) * (CL + (gamma - 1.0) / 2.0 * UL + S);
					C = 2.0 / (gamma + 1.0) * (CL + (gamma - 1.0) / 2.0 * (UL - S));
					D = DL * pow((C / CL), 2.0 / (gamma - 1.0));
					P = PL * pow((C / CL), 2.0 * gamma / (gamma - 1.0));
				}
			}
		}
		else // left shock wave
		{
			PML = PM / PL;
			SL = UL - CL * sqrt((gamma + 1.0) / (2.0 * gamma) * PML + (gamma - 1.0) / (2.0 * gamma));
			if (S <= SL) //point is left from shock wave - left parameters
			{
				D = DL;
				U = UL;
				P = PL;
			}
			else // point is a state with *
			{
				D = DL * (PML + (gamma - 1.0) / (gamma + 1.0)) / (PML * (gamma - 1.0) / (gamma + 1.0) + 1.);
				U = UM;
				P = PM;
			}
		}
	}
	else // point is right from contact
	{
		if (PM > PR) //right from shock wave
		{
			PMR = PM / PR;
			SR = UR + CR * sqrt((gamma + 1.0) / (2.0 * gamma) * PMR + (gamma - 1.0) / (2.0 * gamma));
			if (S >= SR) //point is right from shock wave - right parameters
			{
				D = DR;
				U = UR;
				P = PR;
			}
			else // point is a state with *
			{
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
	double& D, double& U, double& P)
{
	double CL = sqrt(gamma * PL / DL); // left sound speed
	double CR = sqrt(gamma * PR / DR); // right sound speed
	double PM, UM; // pressure and velocity on contact
	double S = 0.;

	// �������� ������� ������
	/*if (DL <= 0.0 || DR <= 0.0 || PL <= 0.0 || PR <= 0.0) {
		std::cerr << "Invalid input to Riemann solver" << std::endl;
		exit(1);
	}*/

	// vacuum generation
	if (is_vacuum(gamma, CL, UL, CR, UR)) {
		std::cout << "Vacuum is generated" << std::endl << "Program stopped" << std::endl;
		exit(228);
	}
	// iteration pressure and velocity
	//void starpu(double gamma, double Quser, double PL, double DL, double UL, double PR, double DR, double UR, double& UM, double& PM)

	starpu(gamma, Quser, PL, DL, UL, PR, DR, UR, UM, PM);

	// found results
	sample(gamma, PL, DL, UL, PR, DR, UR, UM, PM, 0., D, U, P);

	if (P <= 0.0 || D <= 0.0 || !std::isfinite(P) || !std::isfinite(U)) {
		std::cerr << "Invalid output from Riemann solver" << std::endl;
		exit(1);
	}
}

void Godunov_flux_x(double gamma, double p, double vx, double vy, double rho, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
	double m, impx, impy, e;
	//convert_to_conservative(double gamma, double& p, double& vx, double& vy, double& rho, double& m, double& impx, double& impy, double& e)
	convert_to_conservative(gamma, p, vx, vy, rho, m, impx, impy, e);
	Fm = rho * vx;
	Fimpx = Fm * vx + p;
	Fimpy = Fm * vy;
	Fe = (p + e) * vx;
}

void Godunov_flux_y(double gamma, double p, double vx, double vy, double rho, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
	double m, impx, impy, e;
	//convert_to_conservative(double gamma, double& p, double& vx, double& vy, double& rho, double& m, double& impx, double& impy, double& e)
	convert_to_conservative(gamma, p, vx, vy, rho, m, impx, impy, e);
	Fm = rho * vy;
	Fimpx = Fm * vx;
	Fimpy = Fm * vy + p;
	Fe = (p + e) * vy;
}

void Godunov_method_x(double gamma, double Quser, double ml, double impxl, double impyl, double el,
	double mr, double impxr, double impyr, double er,
	double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
	double p, vx, vy, rho;
	double pl, vxl, vyl, rhol;
	double pr, vxr, vyr, rhor;

	//convert_to_conservative(gamma, p, vx, vy, rho, m, impx, impy, e);

	convert_from_conservative(gamma, pl, vxl, vyl, rhol, ml, impxl, impyl, el);
	convert_from_conservative(gamma, pr, vxr, vyr, rhor, mr, impxr, impyr, er);
	/*
	void Riemann_solver(double gamma, double Quser, double PL, double DL, double UL, double PR, double DR, double UR,
	double& D, double& U, double& P)
	*/
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

	//convert_to_conservative(gamma, p, vx, vy, rho, m, impx, impy, e);

	convert_from_conservative(gamma, pd, vxd, vyd, rhod, md, impxd, impyd, ed);
	convert_from_conservative(gamma, pu, vxu, vyu, rhou, mu, impxu, impyu, eu);
	/*
	void Riemann_solver(double gamma, double Quser, double PL, double DL, double UL, double PR, double DR, double UR,
	double& D, double& U, double& P)
	*/

	Riemann_solver(gamma, Quser, pd, rhod, vyd, pu, rhou, vyu, rho, vy, p);

	if (vy >= 0)
		vx = impxd / md;
	else
		vx = impxu / mu;

	//Godunov_flux_x(double gamma, double p, double vx, double vy, double rho, double& Fm, double& Fimpx, double& Fimpy, double& Fe)
	Godunov_flux_y(gamma, p, vx, vy, rho, Fm, Fimpx, Fimpy, Fe);
}

void GK_2d()
{
	Params* params = new Params(100, 40, 1, 0.15, 1.0, 1.0, 1., 1.4, 2.0);

	// ������������ ��� ����� ������� (+- fict � �������) ��������

	vec xc(params->Nx);
	vec x(params->Nx + 1);
	vec yc(params->Ny);
	vec y(params->Ny + 1);

	mtrx p(params->Nx, vec(params->Ny)), vx(params->Nx, vec(params->Ny)), vy(params->Nx, vec(params->Ny)), rho(params->Nx, vec(params->Ny));
	mtrx m(params->Nx, vec(params->Ny)), impx(params->Nx, vec(params->Ny)), impy(params->Nx, vec(params->Ny)), e(params->Nx, vec(params->Ny));
	mtrx m_next(params->Nx, vec(params->Ny)), impx_next(params->Nx, vec(params->Ny)), impy_next(params->Nx, vec(params->Ny)), e_next(params->Nx, vec(params->Ny));


	double mb, impxb, impyb, eb, pb, vxb, vyb, rhob;
	double FmL, FimpxL, FimpyL, FeL, FmR, FimpxR, FimpyR, FeR;

	double pl, vl, rhol, pr, vr, rhor;
	double ml, impxl, impyl, el, mr, impxr, impyr, er;
	double dm, dimp, de;

	int step = 0, max_step = 1001;// 201;
	double time = 0.0;
	double dt;

	int out = 100;
	int minmod_type = 1, bound_left = 1, bound_right = 1, bound_down = 0, bound_up = 0;

	grid(params, x, y, xc, yc);
	init_krujok(params->Nx, params->Ny, params->gamma, xc, yc, p, vx, vy, rho, m, impx, impy, e);

	data_to_file(params->Nx, params->Ny, params->fict, params->gamma, 0.0, x, y, p, vx, vy, rho);

	while (time < params->T && step < max_step)
	{

		dt = get_dt(params, p, rho, vx, vy, x, y);
		time += dt;

		//boundary_cond_x(double p, double vxb, double vyb, double rhob, double& p, double& vx, double& vy, double& rho, int mode);

		/*   START CALC   */
		for (int i = 0; i < params->Nx; i++)
		{
			for (int j = 0; j < params->Ny; j++)
			{
				//		������ ����� �
				// 
				//		����� ����� ����� ����� ������	
				// 

				if (i == 0)
				{
					boundary_cond_x(p[0][j], vx[0][j], vy[0][j], rho[0][j], pb, vxb, vyb, rhob, bound_left);
					convert_to_conservative(params->gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);
					ml = mb;
					impxl = impxb;
					impyl = impyb;
					el = eb;

					mr = m[0][j] - 0.5 * minmod(m[0][j] - mb, m[1][j] - m[0][j], minmod_type);
					impxr = impx[0][j] - 0.5 * minmod(impx[0][j] - impxb, impx[1][j] - impx[0][j], minmod_type);
					impyr = impy[0][j] - 0.5 * minmod(impy[0][j] - impyb, impy[1][j] - impy[0][j], minmod_type);
					er = e[0][j] - 0.5 * minmod(e[0][j] - eb, e[1][j] - e[0][j], minmod_type);
				}

				else if (i == 1)
				{
					boundary_cond_x(p[0][j], vx[0][j], vy[0][j], rho[0][j], pb, vxb, vyb, rhob, bound_left);
					convert_to_conservative(params->gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					ml = m[0][j] + 0.5 * minmod(m[0][j] - mb, m[1][j] - m[0][j], minmod_type);
					impxl = impx[0][j] + 0.5 * minmod(impx[0][j] - impxb, impx[1][j] - impx[0][j], minmod_type);
					impyl = impy[0][j] + 0.5 * minmod(impy[0][j] - impyb, impy[1][j] - impy[0][j], minmod_type);
					el = e[0][j] + 0.5 * minmod(e[0][j] - eb, e[1][j] - e[0][j], minmod_type);

					mr = m[1][j] - 0.5 * minmod(m[i][j] - m[i - 1][j], m[i + 1][j] - m[i][j], minmod_type);
					impxr = impx[1][j] - 0.5 * minmod(impx[i][j] - impx[i - 1][j], impx[i + 1][j] - impx[i][j], minmod_type);
					impyr = impy[1][j] - 0.5 * minmod(impy[i][j] - impy[i - 1][j], impy[i + 1][j] - impy[i][j], minmod_type);
					er = e[1][j] - 0.5 * minmod(e[i][j] - e[i - 1][j], e[i + 1][j] - e[i][j], minmod_type);
				}

				else if (i == params->Nx - 1)
				{
					boundary_cond_x(p[params->Nx - 1][j], vx[params->Nx - 1][j], vy[params->Nx - 1][j], rho[params->Nx - 1][j], pb, vxb, vyb, rhob, bound_right);
					convert_to_conservative(params->gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					ml = m[i - 1][j] + 0.5 * minmod(m[i - 1][j] - m[i - 2][j], m[i][j] - m[i - 1][j], minmod_type);
					impxl = impx[i - 1][j] + 0.5 * minmod(impx[i - 1][j] - impx[i - 2][j], impx[i][j] - impx[i - 1][j], minmod_type);
					impyl = impy[i - 1][j] + 0.5 * minmod(impy[i - 1][j] - impy[i - 2][j], impy[i][j] - impy[i - 1][j], minmod_type);
					el = e[i - 1][j] + 0.5 * minmod(e[i - 1][j] - e[i - 2][j], e[i][j] - e[i - 1][j], minmod_type);

					mr = m[i][j] - 0.5 * minmod(m[i][j] - m[i - 1][j], mb - m[i][j], minmod_type);
					impxr = impx[i][j] - 0.5 * minmod(impx[i][j] - impx[i - 1][j], impxb - impx[i][j], minmod_type);
					impyr = impy[i][j] - 0.5 * minmod(impy[i][j] - impy[i - 1][j], impyb - impy[i][j], minmod_type);
					er = e[i][j] - 0.5 * minmod(e[i][j] - e[i - 1][j], eb - e[i][j], minmod_type);
				}

				else
				{
					ml = m[i - 1][j] + 0.5 * minmod(m[i - 1][j] - m[i - 2][j], m[i][j] - m[i - 1][j], minmod_type);
					impxl = impx[i - 1][j] + 0.5 * minmod(impx[i - 1][j] - impx[i - 2][j], impx[i][j] - impx[i - 1][j], minmod_type);
					impyl = impy[i - 1][j] + 0.5 * minmod(impy[i - 1][j] - impy[i - 2][j], impy[i][j] - impy[i - 1][j], minmod_type);
					el = e[i - 1][j] + 0.5 * minmod(e[i - 1][j] - e[i - 2][j], e[i][j] - e[i - 1][j], minmod_type);

					mr = m[i][j] - 0.5 * minmod(m[i][j] - m[i - 1][j], m[i + 1][j] - m[i][j], minmod_type);
					impxr = impx[i][j] - 0.5 * minmod(impx[i][j] - impx[i - 1][j], impx[i + 1][j] - impx[i][j], minmod_type);
					impyr = impy[i][j] - 0.5 * minmod(impy[i][j] - impy[i - 1][j], impy[i + 1][j] - impy[i][j], minmod_type);
					er = e[i][j] - 0.5 * minmod(e[i][j] - e[i - 1][j], e[i + 1][j] - e[i][j], minmod_type);
				}
				Godunov_method_x(params->gamma, params->Quser, ml, impxl, impyl, el, mr, impxr, impyr, er, FmL, FimpxL, FimpyL, FeL);

				//		����� ����� ������ ����� ������
				if (i == params->Nx - 1)
				{
					boundary_cond_x(p[params->Nx - 1][j], vx[params->Nx - 1][j], vy[params->Nx - 1][j], rho[params->Nx - 1][j], pb, vxb, vyb, rhob, bound_right);
					convert_to_conservative(params->gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					ml = m[i][j] + 0.5 * minmod(m[i][j] - m[i - 1][j], mb - m[i][j], minmod_type);
					impxl = impx[i][j] + 0.5 * minmod(impx[i][j] - impx[i - 1][j], impxb - impx[i][j], minmod_type);
					impyl = impy[i][j] + 0.5 * minmod(impy[i][j] - impy[i - 1][j], impyb - impy[i][j], minmod_type);
					el = e[i][j] + 0.5 * minmod(e[i][j] - e[i - 1][j], eb - e[i][j], minmod_type);

					mr = mb; //m[params->Nx - 1][j];
					impxr = impxb; //impx[params->Nx - 1][j];
					impyr = impyb; //impy[params->Nx - 1][j];
					er = eb; //e[params->Nx - 1][j];

				}
				else if (i == params->Nx - 2)
				{
					boundary_cond_x(p[params->Nx - 1][j], vx[params->Nx - 1][j], vy[params->Nx - 1][j], rho[params->Nx - 1][j], pb, vxb, vyb, rhob, bound_right);
					convert_to_conservative(params->gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					ml = m[i][j] + 0.5 * minmod(m[i][j] - m[i - 1][j], m[i + 1][j] - m[i][j], minmod_type);
					impxl = impx[i][j] + 0.5 * minmod(impx[i][j] - impx[i - 1][j], impx[i + 1][j] - impx[i][j], minmod_type);
					impyl = impy[i][j] + 0.5 * minmod(impy[i][j] - impy[i - 1][j], impy[i + 1][j] - impy[i][j], minmod_type);
					el = e[i][j] + 0.5 * minmod(e[i][j] - e[i - 1][j], e[i + 1][j] - e[i][j], minmod_type);

					mr = m[i + 1][j] - 0.5 * minmod(m[i + 1][j] - m[i][j], mb - m[i + 1][j], minmod_type);
					impxr = impx[i + 1][j] - 0.5 * minmod(impx[i + 1][j] - impx[i][j], impxb - impx[i + 1][j], minmod_type);
					impyr = impy[i + 1][j] - 0.5 * minmod(impy[i + 1][j] - impy[i][j], impyb - impy[i + 1][j], minmod_type);
					er = e[i + 1][j] - 0.5 * minmod(e[i + 1][j] - e[i][j], eb - e[i + 1][j], minmod_type);
				}
				else if (i == 0)
				{
					boundary_cond_x(p[0][j], vx[0][j], vy[0][j], rho[0][j], pb, vxb, vyb, rhob, bound_left);
					convert_to_conservative(params->gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					ml = m[i][j] + 0.5 * minmod(m[i][j] - mb, m[i + 1][j] - m[i][j], minmod_type);
					impxl = impx[i][j] + 0.5 * minmod(impx[i][j] - impxb, impx[i + 1][j] - impx[i][j], minmod_type);
					impyl = impy[i][j] + 0.5 * minmod(impy[i][j] - impyb, impy[i + 1][j] - impy[i][j], minmod_type);
					el = e[i][j] + 0.5 * minmod(e[i][j] - eb, e[i + 1][j] - e[i][j], minmod_type);

					mr = m[i + 1][j] - 0.5 * minmod(m[i + 1][j] - m[i][j], m[i + 2][j] - m[i + 1][j], minmod_type);
					impxr = impx[i + 1][j] - 0.5 * minmod(impx[i + 1][j] - impx[i][j], impx[i + 2][j] - impx[i + 1][j], minmod_type);
					impyr = impy[i + 1][j] - 0.5 * minmod(impy[i + 1][j] - impy[i][j], impy[i + 2][j] - impy[i + 1][j], minmod_type);
					er = e[i + 1][j] - 0.5 * minmod(e[i + 1][j] - e[i][j], e[i + 2][j] - e[i + 1][j], minmod_type);
				}
				else
				{
					ml = m[i][j] + 0.5 * minmod(m[i][j] - m[i - 1][j], m[i + 1][j] - m[i][j], minmod_type);
					impxl = impx[i][j] + 0.5 * minmod(impx[i][j] - impx[i - 1][j], impx[i + 1][j] - impx[i][j], minmod_type);
					impyl = impy[i][j] + 0.5 * minmod(impy[i][j] - impy[i - 1][j], impy[i + 1][j] - impy[i][j], minmod_type);
					el = e[i][j] + 0.5 * minmod(e[i][j] - e[i - 1][j], e[i + 1][j] - e[i][j], minmod_type);

					mr = m[i + 1][j] - 0.5 * minmod(m[i + 1][j] - m[i][j], m[i + 2][j] - m[i + 1][j], minmod_type);
					impxr = impx[i + 1][j] - 0.5 * minmod(impx[i + 1][j] - impx[i][j], impx[i + 2][j] - impx[i + 1][j], minmod_type);
					impyr = impy[i + 1][j] - 0.5 * minmod(impy[i + 1][j] - impy[i][j], impy[i + 2][j] - impy[i + 1][j], minmod_type);
					er = e[i + 1][j] - 0.5 * minmod(e[i + 1][j] - e[i][j], e[i + 2][j] - e[i + 1][j], minmod_type);
				}
				Godunov_method_x(params->gamma, params->Quser, ml, impxl, impyl, el, mr, impxr, impyr, er, FmR, FimpxR, FimpyR, FeR);

				//std::cout << FmL << " " << FmR << std::endl;

				m_next[i][j] = m[i][j] - dt * (FmR - FmL) / params->dx;
				impx_next[i][j] = impx[i][j] - dt * (FimpxR - FimpxL) / params->dx;
				impy_next[i][j] = impy[i][j] - dt * (FimpyR - FimpyL) / params->dx;
				e_next[i][j] = e[i][j] - dt * (FeR - FeL) / params->dx;

				//		������ ����� y
				// 
				//		����� ����� ������ �����
				//
				if (j == 0)
				{
					boundary_cond_x(p[i][0], vx[i][0], vy[i][0], rho[i][0], pb, vxb, vyb, rhob, bound_down);
					convert_to_conservative(params->gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);
					ml = mb;
					impxl = impxb;
					impyl = impyb;
					el = eb;

					mr = m[i][0] - 0.5 * minmod(m[i][0] - mb, m[i][1] - m[i][0], minmod_type);
					impxr = impx[i][0] - 0.5 * minmod(impx[i][0] - impxb, impx[i][1] - impx[i][0], minmod_type);
					impyr = impy[i][0] - 0.5 * minmod(impy[i][0] - impyb, impy[i][1] - impy[i][0], minmod_type);
					er = e[i][0] - 0.5 * minmod(e[i][0] - eb, e[i][1] - e[i][0], minmod_type);
				}

				else if (j == 1)
				{
					boundary_cond_y(p[i][0], vx[i][0], vy[i][0], rho[i][0], pb, vxb, vyb, rhob, bound_down);
					convert_to_conservative(params->gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					mr = m[i][0] + 0.5 * minmod(m[i][0] - mb, m[i][1] - m[i][0], minmod_type);
					impxr = impx[i][0] + 0.5 * minmod(impx[i][0] - impxb, impx[i][1] - impx[i][0], minmod_type);
					impyr = impy[i][0] + 0.5 * minmod(impy[i][0] - impyb, impy[i][1] - impy[i][0], minmod_type);
					er = e[i][0] + 0.5 * minmod(e[i][0] - eb, e[i][1] - e[i][0], minmod_type);

					mr = m[i][1] - 0.5 * minmod(m[i][j] - m[i][j - 1], m[i][j + 1] - m[i][j], minmod_type);
					impxr = impx[i][1] - 0.5 * minmod(impx[i][j] - impx[i][j - 1], impx[i][j + 1] - impx[i][j], minmod_type);
					impyr = impy[i][1] - 0.5 * minmod(impy[i][j] - impy[i][j - 1], impy[i][j + 1] - impy[i][j], minmod_type);
					er = e[i][1] - 0.5 * minmod(e[i][j] - e[i][j - 1], e[i][j + 1] - e[i][j], minmod_type);
				}

				else if (j == params->Ny - 1)
				{
					boundary_cond_y(p[i][params->Ny - 1], vx[i][params->Ny - 1], vy[i][params->Ny - 1], rho[i][params->Ny - 1], pb, vxb, vyb, rhob, bound_up);
					convert_to_conservative(params->gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					ml = m[i][j - 1] + 0.5 * minmod(m[i][j - 1] - m[i][j - 2], m[i][j] - m[i][j - 1], minmod_type);
					impxl = impx[i][j - 1] + 0.5 * minmod(impx[i][j - 1] - impx[i][j - 2], impx[i][j] - impx[i][j - 1], minmod_type);
					impyl = impy[i][j - 1] + 0.5 * minmod(impy[i][j - 1] - impy[i][j - 2], impy[i][j] - impy[i][j - 1], minmod_type);
					el = e[i][j - 1] + 0.5 * minmod(e[i][j - 1] - e[i][j - 2], e[i][j] - e[i][j - 1], minmod_type);

					mr = m[i][j] - 0.5 * minmod(m[i][j] - m[i][j - 1], mb - m[i][j], minmod_type);
					impxr = impx[i][j] - 0.5 * minmod(impx[i][j] - impx[i][j - 1], impxb - impx[i][j], minmod_type);
					impyr = impy[i][j] - 0.5 * minmod(impy[i][j] - impy[i][j - 1], impyb - impy[i][j], minmod_type);
					er = e[i][j] - 0.5 * minmod(e[i][j] - e[i][j - 1], eb - e[i][j], minmod_type);
				}

				else
				{
					ml = m[i][j - 1] + 0.5 * minmod(m[i][j - 1] - m[i][j - 2], m[i][j] - m[i][j - 1], minmod_type);
					impxl = impx[i][j - 1] + 0.5 * minmod(impx[i][j - 1] - impx[i][j - 2], impx[i][j] - impx[i][j - 1], minmod_type);
					impyl = impy[i][j - 1] + 0.5 * minmod(impy[i][j - 1] - impy[i][j - 2], impy[i][j] - impy[i][j - 1], minmod_type);
					el = e[i][j - 1] + 0.5 * minmod(e[i][j - 1] - e[i][j - 2], e[i][j] - e[i][j - 1], minmod_type);

					mr = m[i][j] - 0.5 * minmod(m[i][j] - m[i][j - 1], m[i][j + 1] - m[i][j], minmod_type);
					impxr = impx[i][j] - 0.5 * minmod(impx[i][j] - impx[i][j - 1], impx[i][j + 1] - impx[i][j], minmod_type);
					impyr = impy[i][j] - 0.5 * minmod(impy[i][j] - impy[i][j - 1], impy[i][j + 1] - impy[i][j], minmod_type);
					er = e[i][j] - 0.5 * minmod(e[i][j] - e[i][j - 1], e[i][j + 1] - e[i][j], minmod_type);
				}
				Godunov_method_y(params->gamma, params->Quser, ml, impxl, impyl, el, mr, impxr, impyr, er, FmL, FimpxL, FimpyL, FeL);

				//		����� ����� ������� ����� ������
				if (j == params->Ny - 1)
				{
					boundary_cond_y(p[i][params->Ny - 1], vx[i][params->Ny - 1], vy[i][params->Ny - 1], rho[i][params->Ny - 1], pb, vxb, vyb, rhob, bound_up);
					convert_to_conservative(params->gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					ml = m[i][j] + 0.5 * minmod(m[i][j] - m[i][j - 1], mb - m[i][j], minmod_type);
					impxl = impx[i][j] + 0.5 * minmod(impx[i][j] - impx[i][j - 1], impxb - impx[i][j], minmod_type);
					impyl = impy[i][j] + 0.5 * minmod(impy[i][j] - impy[i][j - 1], impyb - impy[i][j], minmod_type);
					el = e[i][j] + 0.5 * minmod(e[i][j] - e[i][j - 1], eb - e[i][j], minmod_type);

					mr = mb;
					impxr = impxb;
					impyr = impyb;
					er = eb;
				}
				else if (j == params->Ny - 2)
				{
					boundary_cond_y(p[i][params->Ny - 1], vx[i][params->Ny - 1], vy[i][params->Ny - 1], rho[i][params->Ny - 1], pb, vxb, vyb, rhob, bound_up);
					convert_to_conservative(params->gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					ml = m[i][j] + 0.5 * minmod(m[i][j] - m[i][j - 1], m[i][j + 1] - m[i][j], minmod_type);
					impxl = impx[i][j] + 0.5 * minmod(impx[i][j] - impx[i][j - 1], impx[i][j + 1] - impx[i][j], minmod_type);
					impyl = impy[i][j] + 0.5 * minmod(impy[i][j] - impy[i][j - 1], impy[i][j + 1] - impy[i][j], minmod_type);
					el = e[i][j] + 0.5 * minmod(e[i][j] - e[i][j - 1], e[i][j + 1] - e[i][j], minmod_type);

					mr = m[i][j + 1] - 0.5 * minmod(m[i][j + 1] - m[i][j], mb - m[i][j + 1], minmod_type);
					impxr = impx[i][j + 1] - 0.5 * minmod(impx[i][j + 1] - impx[i][j], impxb - impx[i][j + 1], minmod_type);
					impyr = impy[i][j + 1] - 0.5 * minmod(impy[i][j + 1] - impy[i][j], impyb - impy[i][j + 1], minmod_type);
					er = e[i][j + 1] - 0.5 * minmod(e[i][j + 1] - e[i][j], eb - e[i][j + 1], minmod_type);
				}
				else if (j == 0)
				{
					boundary_cond_y(p[i][0], vx[i][0], vy[i][0], rho[i][0], pb, vxb, vyb, rhob, bound_down);
					convert_to_conservative(params->gamma, pb, vxb, vyb, rhob, mb, impxb, impyb, eb);

					ml = m[i][j] + 0.5 * minmod(m[i][j] - mb, m[i][j + 1] - m[i][j], minmod_type);
					impxl = impx[i][j] + 0.5 * minmod(impx[i][j] - impxb, impx[i][j + 1] - impx[i][j], minmod_type);
					impyl = impy[i][j] + 0.5 * minmod(impy[i][j] - impyb, impy[i][j + 1] - impy[i][j], minmod_type);
					el = e[i][j] + 0.5 * minmod(e[i][j] - eb, e[i][j + 1] - e[i][j], minmod_type);

					mr = m[i][j + 1] - 0.5 * minmod(m[i][j + 1] - m[i][j], m[i][j + 2] - m[i][j + 1], minmod_type);
					impxr = impx[i][j + 1] - 0.5 * minmod(impx[i][j + 1] - impx[i][j], impx[i][j + 2] - impx[i][j + 1], minmod_type);
					impyr = impy[i][j + 1] - 0.5 * minmod(impy[i][j + 1] - impy[i][j], impy[i][j + 2] - impy[i][j + 1], minmod_type);
					er = e[i][j + 1] - 0.5 * minmod(e[i][j + 1] - e[i][j], e[i][j + 2] - e[i][j + 1], minmod_type);
				}
				else
				{
					ml = m[i][j] + 0.5 * minmod(m[i][j] - m[i][j - 1], m[i][j + 1] - m[i][j], minmod_type);
					impxl = impx[i][j] + 0.5 * minmod(impx[i][j] - impx[i][j - 1], impx[i][j + 1] - impx[i][j], minmod_type);
					impyl = impy[i][j] + 0.5 * minmod(impy[i][j] - impy[i][j - 1], impy[i][j + 1] - impy[i][j], minmod_type);
					el = e[i][j] + 0.5 * minmod(e[i][j] - e[i][j - 1], e[i][j + 1] - e[i][j], minmod_type);

					mr = m[i][j + 1] - 0.5 * minmod(m[i][j + 1] - m[i][j], m[i][j + 2] - m[i][j + 1], minmod_type);
					impxr = impx[i][j + 1] - 0.5 * minmod(impx[i][j + 1] - impx[i][j], impx[i][j + 2] - impx[i][j + 1], minmod_type);
					impyr = impy[i][j + 1] - 0.5 * minmod(impy[i][j + 1] - impy[i][j], impy[i][j + 2] - impy[i][j + 1], minmod_type);
					er = e[i][j + 1] - 0.5 * minmod(e[i][j + 1] - e[i][j], e[i][j + 2] - e[i][j + 1], minmod_type);
				}
				Godunov_method_y(params->gamma, params->Quser, ml, impxl, impyl, el, mr, impxr, impyr, er, FmR, FimpxR, FimpyR, FeR);
				//std::cout << FmL << " " << FmR << std::endl;

				m_next[i][j] = m_next[i][j] - dt * (FmR - FmL) / params->dy;
				impx_next[i][j] = impx_next[i][j] - dt * (FimpxR - FimpxL) / params->dy;
				impy_next[i][j] = impy_next[i][j] - dt * (FimpyR - FimpyL) / params->dy;
				e_next[i][j] = e_next[i][j] - dt * (FeR - FeL) / params->dy;

			}
		}
		/*    END CALC    */

		// update parameters
		for (int i = 0; i < params->Nx; ++i)
		{
			for (int j = 0; j < params->Ny; ++j)
			{
				m[i][j] = m_next[i][j];
				impx[i][j] = impx_next[i][j];
				impy[i][j] = impy_next[i][j];
				e[i][j] = e_next[i][j];

				convert_from_conservative(params->gamma, p[i][j], vx[i][j], vy[i][j], rho[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
			}
		}
		//data_to_file(int Nx, int Ny, double time, vec x, vec y, mtrx p, mtrx vx, mtrx vy, mtrx rho)
		if (step % out == 0)
			data_to_file(params->Nx, params->Ny, params->fict, params->gamma, time, x, y, p, vx, vy, rho);
		step += 1;
		//std::cout << step << std::endl;
	}

	return;
}


int main()
{
	GK_2d();

	return 0;
}