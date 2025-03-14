#include <iostream>
#include <fstream>
#include <string> 
#include <vector>
#include <math.h>


void writeCSV(std::string filename, int iter, std::vector<double> grid, std::vector<double> u, std::vector<double> P, std::vector<double> rho, double t, int Nx, int fict) {
	std::string name = "CSVs\\";
	name += filename;
	name += "\\Iter=";
	name += std::to_string(iter);
	name += ".csv";

	std::ofstream outfile(name);
	std::string tmp;
	tmp = "Time,X,Rho,P,U\n";
	outfile << tmp;
	for (int i = fict + 1; i < Nx - fict; i++) {
		tmp = "";
		tmp += std::to_string(t) + ',';
		tmp += std::to_string((grid[i] + grid[i + 1]) / 2.) + ',';
		tmp += std::to_string(rho[i]) + ',';
		tmp += std::to_string(P[i]) + ',';
		tmp += std::to_string((u[i] + u[i + 1]) / 2.) + '\n';
		outfile << tmp;
	}
	outfile.close();
}


void Boundary(std::vector<double>* u, double v0) {
	(*u)[0] = v0;
	(*u)[1] = v0;
	(*u)[u->size() - 1] = (*u)[u->size() - 2];
}


void Artificial_viscosity(int Nx, std::vector<double> rho, std::vector<double> vb, std::vector<double>* dP) {
	const double mu_0 = 8.;
	for (int i = 0; i < Nx; ++i) {
		if ((vb[i + 1] - vb[i]) < 0) {
			(*dP)[i] = -mu_0 * 0.5 * rho[i] * std::abs(vb[i + 1] + vb[i]) * (vb[i + 1] - vb[i]);
		}
		else {
			(*dP)[i] = 0.0;
		}
	}
}


int main() {
	int N = 200, fict = 1; 
	int Nx = N + 2 * fict;
	int iter = 0, Nt = 500, iterwrite = 5;
	double L = 1., rho0 = 2., P0 = 101325.;
	double Q = 2 * 1.e6, gamma = 1.4;

	double t = 0., t_end = 1.e-2, tau;
	double CFL = 0.3;
	double h = L / N, rho_cj, u_cj;
	rho_cj = (gamma + 1) / gamma * rho0;
	u_cj = sqrt(2 * (gamma - 1) / (gamma + 1) * Q);

	std::string filename = "KrestDetonate";

	std::vector<double> grid(Nx + 1);
	std::vector<double> u(Nx + 1), P(Nx), rho(Nx);
	std::vector<double> mass(Nx), I(Nx), dP_vis(Nx), W(Nx), Temp(Nx);

	std::vector<double> new_grid(Nx + 1);
	std::vector<double> new_u(Nx + 1), new_P(Nx), new_rho(Nx);
	std::vector<double> new_I(Nx), new_W(Nx), new_Temp(Nx);

	// init non-zeros massives
	for (int j = 0; j < Nx + 1; j++) {
		grid[j] = h * (j - fict);
	}

	Boundary(&u, u_cj);
	
	for (int j = 0; j < Nx; j++) {
		rho[j] = rho0;
		mass[j] = rho[j] * (grid[j + 1] - grid[j]);
		P[j] = P0;
		W[j] = 1.;
		I[j] = 0.5 * rho[j] * (pow(u[j], 2.0)) + P[j] / (gamma - 1.0);
	}

	writeCSV(filename, iter, grid, u, P, rho, t, Nx, fict);

	while (t < t_end && iter < Nt) {
		// get dt
		tau = 1.0e3;
		double tmp;
		for (int j = 1; j < N; ++j) {
			double sz = sqrt(gamma * abs(P[j] / rho[j]));
			tmp = CFL * (grid[j + 1] - grid[j]) / (abs(u[j]) + sz);
			if (tmp < tau) {
				tau = tmp;
			}
		}

		Artificial_viscosity(Nx, rho, u, &dP_vis);

		// velocity
		for (int j = 1; j < Nx; j++) {
			new_u[j] = u[j] + tau / (0.5 * (mass[j - 1] + mass[j])) * (P[j - 1] - P[j] + dP_vis[j - 1] - dP_vis[j]);
		}
		Boundary(&new_u, u_cj);

		// Euler coordinates
		for (int j = 0; j < Nx + 1; j++) {
			new_grid[j] = grid[j] + new_u[j] * tau;
		}

		// density
		for (int j = fict; j < Nx - fict; j++) {
			new_rho[j] = mass[j] / abs(new_grid[j + 1] - new_grid[j]);
		}
		new_rho[0] = new_rho[1];
		new_rho[Nx - fict] = new_rho[Nx - fict - 1];

		// energy
		for (int j = fict; j < Nx - fict; j++) {
			new_I[j] = I[j] + tau / mass[j] * (new_u[j] * ((mass[j] * P[j - 1] + mass[j - 1] * P[j]) /
 			 (mass[j] + mass[j - 1]) + 0.5 * (dP_vis[j] + dP_vis[j - 1])) - new_u[j + 1] * ((mass[j + 1] * P[j] + mass[j] * P[j + 1]) /
  			  (mass[j] + mass[j + 1]) + 0.5 * (dP_vis[j] + dP_vis[j + 1]))) + 1.0 / 8.0 * (pow(u[j + 1] + u[j], 2) - pow(new_u[j + 1] + new_u[j], 2));	
		}

		// massovaya dolya
		for (int j = 0; j < Nx; j++) {
			new_W[j] = 1. - (1. / rho0 - 1. / new_rho[j]) / (1. / rho0 - 1. / rho_cj);
			if (new_W[j] > W[j] && new_W[j] < 0.9) {
				new_W[j] = 0.;
			}
			if (new_W[j] < 0.0) {
				new_W[j] = 0.0;
			}
			if (W[j] <= 1.e-8) {
				new_W[j] = 0.0;
			}
		}
		new_W[0] = new_W[1];
		new_W[N + fict] = new_W[N + fict - 1];

		for (int j = fict; j < Nx; j++) {
			if (new_W[j] < 0.99) {
				new_P[j] = (1 - new_W[j]) * (gamma - 1) * new_rho[j] * (new_I[j] + Q);
			}
			else {
				new_P[j] = P0;
			}
		}

		// temperature
		/*
		for (int i = 1; i < Nx; i++)
		{
			Temp_new[i] = Temp[i] + dt / (0.5 * (m[i - 1] + m[i])) * (v[i + 1] - v[i]);
		}
		*/
		
		// updating values
		for (int j = 0; j < Nx + 1; j++) {
			u[j] = new_u[j];
			grid[j] = new_grid[j];
		}
		for (int j = 0; j < Nx; j++) {
			P[j] = new_P[j];
			rho[j] = new_rho[j];
			W[j] = new_W[j];
			I[j] = new_I[j];
			Temp[j] = new_Temp[j];
		}

		t += tau;
		iter++;
		if (iter % iterwrite == 0)
		{
			writeCSV(filename, iter, grid, u, P, rho, t, Nx, fict);
		}
	}
	writeCSV(filename, iter, grid, u, P, rho, t, Nx, fict);

	if (t >= t_end) {
		std::cout << "Solve stopped by time\n";
		std::cout << "On iter = " << iter << std::endl;
	}
	else {
		std::cout << "Solve stopped by iterations number";
	}
	return 0;
}