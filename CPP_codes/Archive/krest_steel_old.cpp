#include <iostream>
#include<fstream>
#include <string>
#include <vector>
#include <functional>


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
    const double mu_0 = 2. / u0;
    for (int i = 0; i < Nx; ++i) {
        if ((vb[i + 1] - vb[i]) < 0) {
            (*dP)[i] = -mu_0 * 0.5 * rho[i] * std::abs(vb[i + 1] + vb[i]) * (vb[i + 1] - vb[i]);
        }
        else {
            (*dP)[i] = 0.0;
        }
    }
}


void krest_steel() {
    int N = 300;
    int fict = 1;

    int Nx = N + 2 * fict;
    int Nt = 1000;
    double L = 1.;
    double t_start = 0.;
    double t_end = 0.0001;
    double tau;
    double h_start = L / N;

    double mu = 79.3 * 1e9; //Pa
    double Y0 = 0.3 * 1e9; //Pa
    double K = 156. * 1e9;
    double rho0 = 7850.;
    double u0 = 15.;
    std::string filename = "KrestSteel";

    double CFL = 0.3;

    std::vector<double> grid(Nx + 1);
    for (int j = 0; j < Nx + 1; j++) {
        grid[j] = h_start * (j - fict);
    }

    std::vector<double> u(Nx + 1), P(Nx), rho(Nx), mass(Nx), I(Nx), dP_vis(Nx), S_x(Nx);

    //initial conditions
    for (int j = 0; j < Nx + 1; j++) {
        u[j] = 0;
        if (j != Nx) {
            P[j] = 0;
            rho[j] = rho0;
        }
    }
    Boundary(&u, u0);

    for (int j = 0; j < Nx; j++) {
        I[j] = 0;
        mass[j] = (grid[j + 1] - grid[j]) * rho[j];
        S_x[j] = 0;
    }
    std::string title = "CSVs\\";
    title += filename;
    title += "\\Iter=0.csv";
    writeCSV(title, grid, u, P, rho, t_start, Nx, 1);

    double t = 0;
    int iter = 0;
    int iterwrite = 50;
    while (t < t_end && iter < Nt) {
        tau = 1e3;
        for (int j = 0; j < Nx; j++) {
            double sz = std::sqrt(K * rho0) / rho[j];
            sz = 5000.;
            tau = std::min(tau, CFL * (grid[j + 1] - grid[j]) / (std::abs(u[j]) + sz));
        }

        std::vector<double> new_grid(Nx + 1), new_u(Nx + 1), new_P(Nx), new_rho(Nx), new_mass(Nx), new_I(Nx), new_S_x(Nx);

        Artificial_viscosity(Nx, rho, u, &dP_vis, u0);

        for (int j = fict; j < Nx; j++) {
            new_u[j] = u[j] + tau / (0.5 * (mass[j - 1] + mass[j])) * (P[j - 1] - P[j] + dP_vis[j - 1] - dP_vis[j]);
        }
        Boundary(&new_u, u0);

        for (int j = 0; j < Nx + 1; j++) {
            new_grid[j] = new_u[j] * tau + grid[j];
        }

        for (int j = fict; j < Nx; j++) {
            new_rho[j] = mass[j] / (new_grid[j + 1] - new_grid[j]);
        }
        new_rho[0] = new_rho[1];
        new_rho[Nx - 2] = new_rho[Nx - 3];

        for (int j = fict; j < Nx - fict; j++) {
            double first, second, third;
            first = (mass[j] * P[j - 1] + mass[j - 1] * P[j]) / (mass[j] + mass[j - 1]) + 0.5 * (dP_vis[j] + dP_vis[j - 1]);
            second = (mass[j + 1] * P[j + 1] + mass[j] * P[j + 1]) / (mass[j + 1] + mass[j]) + 0.5 * (dP_vis[j] + dP_vis[j + 1]);
            third = (u[j + 1] + u[j]) * (u[j + 1] + u[j]) - (new_u[j + 1] + new_u[j]) * (new_u[j + 1] + new_u[j]);
            new_I[j] = I[j] + tau / mass[j] * (first * new_u[j] - new_u[j + 1] * second) + 1. / 8. * third;
        }

        for (int j = fict; j < Nx - fict; j++) {
            new_S_x[j] = S_x[j] + 2 * mu * (-tau * (new_u[j + 1] - new_u[j]) / (grid[j + 1] - grid[j]) + 2.0 / 3.0 * (1. / new_rho[j] - 1. / rho[j]) / (1. / new_rho[j] + 1. / rho[j]));

        }
        new_S_x[0] = new_S_x[1];
        new_S_x[Nx - 2] = new_S_x[Nx - 3];

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

        for (int j = 0; j < Nx + 1; j++) {
            u[j] = new_u[j];
        }
        for (int j = 0; j < Nx - 1; j++) {
            rho[j] = new_rho[j];
            P[j] = new_P[j];
            I[j] = new_I[j];
        }

        for (int j = 0; j < Nx + 1; j++) {
            grid[j] = new_grid[j];
        }


        t += tau;
        iter++;
        if (iter % iterwrite == 0) {
            title = "CSVs\\";
            title += filename;
            title += "\\Iter=";
            title += std::to_string(iter);
            title += ".csv";
            writeCSV(title, grid, u, P, rho, t, Nx, 1);
        }
    }
    if (t >= t_end) {
        std::cout << "Solve stopped by time\n";
        std::cout << "On iter = " << iter << std::endl;
    }
    else {
        std::cout << "Solve stopped by iterations number";
    }
    title = "CSVs\\";
    title += filename;
    title += "\\Iter=";
    title += std::to_string(iter);
    title += ".csv";
    writeCSV(title, grid, u, P, rho, t, Nx, 1);
}


int main() {
    krest_steel();
}