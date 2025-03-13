#include <iostream>
#include <fstream>
#include <string> 
#include <vector>


void writeCSV(std::string name, std::vector<double> grid, std::vector<double> u, std::vector<double> P, std::vector<double> rho, double t, int Nx, int fict) {
    std::ofstream outfile(name);
    std::string tmp;
    tmp = "Time,X,Rho,P,U\n";
    outfile << tmp;
    for (int i = fict; i < Nx - fict; i++) {
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


void Artificial_viscosity(int Nx, std::vector<double> rho, std::vector<double> vb, std::vector<double>* dP, double u0) {
	const double mu_0 = 600. / u0;
	for (int i = 0; i < Nx; ++i) {
		if ((vb[i + 1] - vb[i]) < 0) {
			(*dP)[i] = -mu_0 * 0.5 * rho[i] * std::abs(vb[i + 1] + vb[i]) * (vb[i + 1] - vb[i]);
		}
		else {
			(*dP)[i] = 0.0;
		}
	}
}


int krest_tension() {
	int N = 300;
	int fict = 1;
	int Nx = N + 2 * fict;
	int Nt = 1000;
	double Y0 = 0.3 * 1e9, K = 156. * 1e9, mu = 79.3 * 1e9, rho0 = 7850.;
	double u0 = 150.;

	double t_start = 0., t_end = 0.0001, L = 1., CFL = 0.3;
	double h_start = L / N;
	double tau;

	int iter = 0;
	int iterwrite = 20;
	std::string filename = "KrestSteel";

	std::vector<double> grid(Nx + 1), center_grid(Nx);
	std::vector<double> u(Nx + 1), P(Nx), rho(Nx, 7850.);
	std::vector<double> mass(Nx), I(Nx), dP_vis(Nx), S_x(Nx), Eps(Nx);

	// new massives
	std::vector<double> new_grid(Nx + 1);
	std::vector<double> new_u(Nx + 1), new_P(Nx), new_rho(Nx), new_I(Nx), new_S_x(Nx), new_Eps(Nx);

	// init non-zeros massives
	for (int i = 0; i < Nx + 1; i++) {
		grid[i] = h_start * (i - fict);
	}
	for (int i = 0; i < Nx; i++) {
		rho[i] = rho0;
		mass[i] = rho[i] * (grid[i + 1] - grid[i]);
	}

	std::string title = "CSVs\\";
	title += filename;
	title += "\\Iter=0.csv";
	writeCSV(title, grid, u, P, rho, t_start, Nx, 1);

	while (t_start < t_end && iter < Nt) {
		// get dt
		tau = 1.e3;
		for (int j = fict; j < Nx - fict; j++) {
			tau = std::min(CFL * abs(grid[j + 1] - grid[j]) / (abs(u[j]) + 5000.), tau);
		}

		// art viscosity
		Artificial_viscosity(Nx, rho, u, &dP_vis, u0);

		// velocity
		for (int j = fict; j < Nx; j++){
			new_u[j] = u[j] + tau / (0.5 * (mass[j - 1] + mass[j])) * (P[j - 1] - P[j] + dP_vis[j - 1] - dP_vis[j]);
		}
		Boundary(&new_u, u0);

		// Euler coordinates
		for (int j = 0; j < Nx ; j++) {
			new_grid[j] = grid[j] + new_u[j] * tau;
		}

		// density
		for (int j = fict; j < Nx; j++) {
			new_rho[j] = mass[j] / (new_grid[j + 1] - new_grid[j]);
		}
		new_rho[0] = new_rho[1];
		new_rho[N + fict] = new_rho[N + fict - 1];

		// energy
		for (int j = fict; j < Nx - 1; j++) {
			double first, second, third;
			first = (mass[j] * P[j - 1] + mass[j - 1] * P[j]) / (mass[j] + mass[j - 1]) + 0.5 * (dP_vis[j] + dP_vis[j - 1]);
			second = (mass[j + 1] * P[j + 1] + mass[j] * P[j + 1]) / (mass[j + 1] + mass[j]) + 0.5 * (dP_vis[j] + dP_vis[j + 1]);
			third = (u[j + 1] + u[j]) * (u[j + 1] + u[j]) - (new_u[j + 1] + new_u[j]) * (new_u[j + 1] + new_u[j]);
			new_I[j] = I[j] + tau / mass[j] * (first * new_u[j] - new_u[j + 1] * second) + 1. / 8. * third;
		}

		// Sx
		for (int j = fict; j < Nx - 1; j++) {
			new_S_x[j] = S_x[j] + 2 * mu * (-tau * (new_u[j + 1] - new_u[j]) / (grid[j + 1] - grid[j]) + 2.0 / 3.0 * (1. / new_rho[j] - 1. / rho[j]) / (1. / new_rho[j] + 1. / rho[j]));
		}
		new_S_x[0] = new_S_x[1];
		new_S_x[N + fict] = new_S_x[N + fict - 1];

		// p
		for (int j = fict; j < Nx; j++) {
			if (abs(new_S_x[j]) < abs(2. / 3. * Y0)) {
				new_P[j] = K * (1 - rho0 / new_rho[j]) + new_S_x[j];
			}
			else {
				if (new_S_x[j] > 0) {
					new_P[j] = K * (1. - rho0 / new_rho[j]) + 2. / 3. * Y0;
				}
				else {
					new_P[j] = K * (1. - rho0 / new_rho[j]) - 2. / 3. * Y0;
				}
			}
		}
		new_P[0] = new_P[1];


		// new data
		for (int i = 0; i < Nx; i++)
		{
			u[i] = new_u[i];
			P[i] = new_P[i];
			if (new_u[i] != new_u[i]) { std::cerr << "nan ......." << std::endl; exit(99); }
			grid[i] = new_grid[i];
			rho[i] = new_rho[i];
			S_x[i] = new_S_x[i];
			I[i] = new_I[i];
		}

		t_start += tau;
		iter++;
		if (iter % iterwrite == 0)
		{
			title = "CSVs\\";
            title += filename;
            title += "\\Iter=";
            title += std::to_string(iter);
            title += ".csv";
			writeCSV(title, grid, u, P, rho, t_start, Nx, 1);
		}
	}
	title = "CSVs\\";
	title += filename;
	title += "\\Iter=";
	title += std::to_string(iter);
	title += ".csv";
	writeCSV(title, grid, u, P, rho, t_start, Nx, 1);

	if (t_start >= t_end) {
		std::cout << "Solve stopped by time\n";
		std::cout << "On iter = " << iter << std::endl;
	}
	else {
		std::cout << "Solve stopped by iterations number";
	}
	return 0;
}


int main() {
	krest_tension();
}