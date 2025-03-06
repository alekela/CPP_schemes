#include <iostream>
#include<fstream>
#include <string>
#include <vector>
#include <functional>
#include <cmath>


void writeCSV(std::string name, std::vector<long double> grid_P, std::vector<long double> u, std::vector<long double> P, std::vector<long double> rho, std::vector<long double> W, long double t, int fict) {
    std::ofstream outfile(name);
    std::string tmp;
    tmp = "Time,X,Rho,P,U,W\n";
    outfile << tmp;
    for (int i = fict; i < P.size() - fict; i++) {
        tmp = "";
        tmp += std::to_string(t) + ',';
        tmp += std::to_string(grid_P[i]) + ',';
        tmp += std::to_string(rho[i]) + ',';
        tmp += std::to_string(P[i]) + ',';
        tmp += std::to_string(u[i]) + ',';
        tmp += std::to_string(W[i]) + '\n';
        outfile << tmp;
    }
    outfile.close();
}


void Boundary(int state, std::vector<long double>* u, long double v0) {
    (*u)[0] = v0;
    (*u)[1] = v0;
    if (state == 0) { // wall
        (*u)[u->size() - 1] = 0;
        (*u)[u->size() - 2] = 0;
    }
    else if (state == 1) { // free end
        (*u)[u->size() - 1] = (*u)[u->size() - 2];
    }
}


void Artificial_viscosity(int Nx, long double* rho, long double* vb, long double* dP, int state) { //s - mass
    const long double nu_0 = 1.0;
    const long double mu_0 = 2e-5;
    long double fict = 1;
    if (state == 1) { // linear
        for (int i = 0; i < Nx; ++i) {
            if ((vb[i + 1] - vb[i]) < 0) {
                dP[i] = -nu_0 * rho[i] * (vb[i + 1] - vb[i]);
            }
            else {
                dP[i] = 0.0;
            }
        }
    }
    else if (state == 2) { // Latter
        for (int i = 0; i < Nx; ++i) {
            if ((vb[i + 1] - vb[i]) < 0) {
                dP[i] = -mu_0 * 0.5 * rho[i] * (vb[i + 1] + vb[i]) * (vb[i + 1] - vb[i]);
            }
            else {
                dP[i] = 0.0;
            }
        }
    }
}


void krest_detonate(int fict) {
    int Nx = 1002;
    int Nt = 250;
    long double x_start = 0;
    long double x_end = 100;
    long double t_start = 0;
    long double t_end = 0.2;
    long double tau;
    long double h_start = (x_end - x_start) / (Nx - 2);
    long double gamma = 1.4;
    long double P0 = 1e5; // Pa
    long double rho0 = 44;
    long double Q = 2e6; // J/kg
    long double Vcj = gamma / (gamma + 1) / rho0;
    long double u0 = sqrt(2 * (gamma - 1) / (gamma + 1) * Q);
    std::string filename = "KrestDetonate";

    long double CFL = 0.5;

    // Boundary conditions
    long double u_left, u_right, P_left, P_right, rho_left, rho_right;
    // Test1(&rho_left, &u_left, &P_left, &rho_right, &u_right, &P_right);


    std::vector<long double> grid(Nx + 1);
    std::vector<long double> center_grid(Nx);
    for (int j = 0; j < Nx + 1; j++) {
        grid[j] = (x_start + h_start * (j - fict));
        if (j != Nx) {
            center_grid[j] = (x_start + h_start * (j - fict) + h_start / 2);
        }
    }

    /*
    // theoretic graphs
    std::vector<long double> theoretic_P(Nx), theoretic_rho(Nx), theoretic_u(Nx);
    for (int j = 0; j < Nx; j++) {
        AnalitSolveTest1(center_grid[j], t_end, &theoretic_P[j], &theoretic_u[j], &theoretic_rho[j]);
    }
    writeCSV("CSVs\\Theory_graph.csv", center_grid, theoretic_u, theoretic_P, theoretic_rho, t_end);
    */

    std::vector<long double> u(Nx + 1), P(Nx), rho(Nx), mass(Nx), I(Nx), dP_vis(Nx), W(Nx), T(Nx);

    //initial conditions
    for (int j = 0; j < Nx + 1; j++) {
        u[j] = 0;
        if (j != Nx) {
            P[j] = P0;
            rho[j] = rho0;
        }
    }
    Boundary(0, &u, u0);

    for (int j = fict; j < Nx - fict; j++) {
        I[j] = P[j] / (gamma - 1) / rho[j];
        W[j] = 1;
        mass[j] = (h_start * rho[j]);
    }
    std::string title = "CSVs\\";
    title += filename;
    title += "\\Iter=0.csv";
    writeCSV(title, center_grid, u, P, rho, W, t_start, 1);

    long double t = 0;
    int iter = 0;
    int iterwrite = 50;
    while (t < t_end && iter < Nt) {
        tau = 1e6;
        for (int j = fict; j < Nx - fict; j++) {
            long double sz = sqrt(gamma * P[j] / rho[j]);
            double temp = CFL * (grid[j + 1] - grid[j]) / (std::abs(u[j] + u[j + 1]) / 2. + sz);
            if (temp < tau) {
                tau = temp;
            }
        }

        std::vector<long double> new_u(Nx + 1), new_P(Nx), new_rho(Nx), new_mass(Nx), new_I(Nx), new_W(Nx), new_T(Nx);

        Artificial_viscosity(Nx, rho.data(), u.data(), dP_vis.data(), 2); // 1 - linear, 2 - Latter

        for (int j = fict + 1; j < Nx + 1 - fict; j++) {
            new_u[j] = u[j] - tau * 2. / (mass[j] + mass[j - 1]) * (P[j] + dP_vis[j] - P[j - 1] - dP_vis[j - 1]);
        }
        Boundary(0, &new_u, u0);

        for (int j = 0; j < Nx + 1; j++) {
            grid[j] = new_u[j] * tau + grid[j];
        }

        for (int j = 0; j < Nx; j++) {
            center_grid[j] = (grid[j] + grid[j + 1]) / 2.;
        }

        for (int j = fict; j < Nx - fict; j++) {
            new_rho[j] = mass[j] / (grid[j + 1] - grid[j]);
        }

        for (int j = fict; j < Nx - fict; j++) {
            long double first, second, third;
            first = (mass[j] * P[j - 1] + mass[j - 1] * P[j]) / (mass[j] + mass[j - 1]) + 0.5 * (dP_vis[j] + dP_vis[j - 1]);
            second = (mass[j + 1] * P[j] + mass[j] * P[j + 1]) / (mass[j + 1] + mass[j]) + 0.5 * (dP_vis[j] + dP_vis[j + 1]);
            third = (u[j + 1] + u[j]) * (u[j + 1] + u[j]) - (new_u[j + 1] + new_u[j]) * (new_u[j + 1] + new_u[j]);
            new_I[j] = I[j] + tau / mass[j] * (first * new_u[j] - new_u[j + 1] * second) + 1. / 8. * third;
        }

        long double tmp;
        for (int j = fict; j < Nx - fict; j++) {
            tmp = 1. - (1. / rho0 - 1. / new_rho[j]) / (1. / rho0 - Vcj);
            if ((tmp > W[j]) && (tmp < 0.9)) {
                new_W[j] = 0;
            }
            else if (tmp < 0) {
                new_W[j] = 0;
            }
            else if (tmp <= W[j]) {
                new_W[j] = tmp;
            }
        }

        for (int j = fict; j < Nx - fict; j++) {
            if (new_W[j] < 0.99) {
                new_P[j] = (1 - new_W[j]) * (gamma - 1) * new_rho[j] * (new_I[j] + Q);
            }
            else {
                new_P[j] = P0;
            }
        }

        for (int j = fict; j < Nx + 1 - fict; j++) {
            u[j] = new_u[j];
        }
        for (int j = fict; j < Nx - fict; j++) {
            rho[j] = new_rho[j];
            P[j] = new_P[j];
            I[j] = new_I[j];
            W[j] = new_W[j];
        }

        t += tau;
        iter++;
        if (iter % iterwrite == 0) {
            title = "CSVs\\";
            title += filename;
            title += "\\Iter=";
            title += std::to_string(iter);
            title += ".csv";
            writeCSV(title, center_grid, u, P, rho, W, t, 1);
        }
    }
    if (t >= t_end) {
        std::cout << "Solve stopped by time";
    }
    else {
        std::cout << "Solve stopped by iterations number";
    }
    title = "CSVs\\";
    title += filename;
    title += "\\Iter=";
    title += std::to_string(iter);
    title += ".csv";
    writeCSV(title, center_grid, u, P, rho, W, t, 1);
}


int main() {
    krest_detonate(1);
}