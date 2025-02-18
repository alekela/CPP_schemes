#include <iostream>
#include<fstream>
#include <string>
#include <vector>
#include <functional>
#include <cmath>


void writeCSV(std::string name, std::vector<double> grid_P, std::vector<double> u, std::vector<double> P, std::vector<double> rho, double t, int fict) {
    std::ofstream outfile(name);
    std::string tmp;
    tmp = "Time,X,Rho,P,U\n";
    outfile << tmp;
    for (int i = fict; i < P.size() - fict; i++) {
        tmp = "";
        tmp += std::to_string(t) + ',';
        tmp += std::to_string(grid_P[i]) + ',';
        tmp += std::to_string(rho[i]) + ',';
        tmp += std::to_string(P[i]) + ',';
        tmp += std::to_string(u[i]) + '\n';
        outfile << tmp;
    }
    outfile.close();
}


void Boundary(int state, std::vector<double>* u, double v0) {
    (*u)[0] = v0;
    (*u)[1] = v0;
    if (state == 0) { // free end
        (*u)[u->size() - 2] = 0;
    }
    else if (state == 1) { // reflection
        (*u)[u->size() - 2] = -(*u)[u->size() - 1];
    }
}


void Artificial_viscosity(int Nx, double* rho, double* vb, double* dP, int state) { //s - mass
    const double nu_0 = 1.0;
    const double mu_0 = 2.0;
    double fict = 1;
    if (state == 1) { // linear
        for (int i = 1; i < Nx - 1; ++i) {
            if ((vb[i + 1] - vb[i]) < 0) {
                dP[i] = -nu_0 * rho[i] * (vb[i + 1] - vb[i]);
            }
            else {
                dP[i] = 0.0;
            }
        }
        dP[0] = dP[1];
        dP[Nx - 1] = dP[Nx - 2];
    }
    else if (state == 2) { // Latter
        for (int i = 1; i < Nx - 1; ++i) {
            if ((vb[i + 1] - vb[i]) < 0) {
                dP[i] = -mu_0 * rho[i] * std::abs(vb[i + 1] - vb[i]) * (vb[i + 1] - vb[i]);
            }
            else {
                dP[i] = 0.0;
            }
        }
        dP[0] = dP[1];
        dP[Nx - 1] = dP[Nx - 2];
    }
}


void krest_detonate() {
    int Nx = 1000;
    int Nt = 5000;
    double x_start = 0;
    double x_end = 1;
    double t_start = 0;
    double t_end = 0.2;
    double tau;
    double h_start = (x_end - x_start) / Nx;
    double gamma = 1.4;
    double mu = 79.3 * 1.e9; //Pa
    double Y0 = 0.3 * 1.e9; //Pa
    double PLAP = 1;
    double P0 = 100000; //Pa
    double K = 156. * 1.e9;
    double rho0 = 7850;
    double Q = 2000000;
    double Vcj = gamma / (gamma + 1) / rho0;
    double u0 = sqrt(2 * (gamma - 1) / (gamma + 1) * Q);
    std::string filename = "KrestDetonate";

    double CFL = 0.5;

    // Boundary conditions
    double u_left, u_right, P_left, P_right, rho_left, rho_right;
    // Test1(&rho_left, &u_left, &P_left, &rho_right, &u_right, &P_right);


    std::vector<double> grid(Nx + 1);
    std::vector<double> center_grid(Nx);
    for (int j = 0; j < Nx + 1; j++) {
        grid[j] = (x_start + h_start * j);
        if (j != Nx) {
            center_grid[j] = (x_start + h_start * j + h_start / 2);
        }
    }

    /*
    // theoretic graphs
    std::vector<double> theoretic_P(Nx), theoretic_rho(Nx), theoretic_u(Nx);
    for (int j = 0; j < Nx; j++) {
        AnalitSolveTest1(center_grid[j], t_end, &theoretic_P[j], &theoretic_u[j], &theoretic_rho[j]);
    }
    writeCSV("CSVs\\Theory_graph.csv", center_grid, theoretic_u, theoretic_P, theoretic_rho, t_end);
    */

    std::vector<double> u(Nx + 1), P(Nx), rho(Nx), mass(Nx), I(Nx), dP_vis(Nx), W(Nx), T(Nx);

    //initial conditions
    for (int j = 0; j < Nx + 1; j++) {
        u[j] = 0;
        if (j != Nx) {
            P[j] = P0;
            rho[j] = rho0;
        }
    }
    u[0] = u0;
    u[1] = u0;

    for (int j = 0; j < Nx; j++) {
        I[j] = 0;
        mass[j] = (h_start * rho[j]);
    }
    std::string title = "CSVs\\";
    title += filename;
    title += "\\Iter=000.csv";
    writeCSV(title, center_grid, u, P, rho, t_start, 1);

    double t = 0;
    int iter = 0;
    int iterwrite = 100;
    while (t < t_end && iter < Nt) {
        tau = (t_end - t_start) / Nt;
        for (int j = 0; j < Nx; j++) {
            // double sz = sound_speed(P[j], rho[j], gamma);
            tau = std::min(tau, CFL * (grid[j + 1] - grid[j]) / (std::abs(u[j]) + 5000.));
        }

        std::vector<double> new_u(Nx + 1), new_P(Nx), new_rho(Nx), new_mass(Nx), new_I(Nx), new_W(Nx), new_T(Nx);

        Artificial_viscosity(Nx, rho.data(), u.data(), dP_vis.data(), 1); // 1 - linear, 2 - Latter

        for (int j = 1; j < Nx + 1; j++) {
            new_u[j] = ((-tau * 2. / (mass[j] + mass[j - 1]) * (P[j] - P[j - 1] + dP_vis[j] - dP_vis[j - 1])) + u[j]);
        }
        Boundary(1, &new_u, u0);

        for (int j = 1; j < Nx; j++) {
            grid[j] = new_u[j] * tau + grid[j];
        }

        for (int j = 0; j < Nx; j++) {
            center_grid[j] = (grid[j] + grid[j + 1]) / 2.;
        }

        for (int j = 1; j < Nx - 1; j++) {
            new_rho[j] = mass[j] / (grid[j + 1] - grid[j]);
        }
        new_rho[0] = new_rho[1];

        for (int j = 1; j < Nx - 1; j++) {
            double first, second, third;
            first = (mass[j] * P[j - 1] + mass[j - 1] * P[j]) / (mass[j] + mass[j - 1]) + 0.5 * (dP_vis[j] + dP_vis[j - 1]);
            second = (mass[j + 1] * P[j + 1] + mass[j] * P[j + 1]) / (mass[j + 1] + mass[j]) + 0.5 * (dP_vis[j] + dP_vis[j + 1]);
            third = (u[j + 1] + u[j]) * (u[j + 1] + u[j]) - (new_u[j + 1] + new_u[j]) * (new_u[j + 1] + new_u[j]);
            new_I[j] = I[j] + tau / mass[j] * (first * new_u[j] - new_u[j + 1] * second) + 1. / 8. * third;
        }

        double tmp;
        for (int j = 1; j < Nx - 1; j++) {
            tmp = 1 - (1 / rho0 - 1 / new_rho[j]) / (1 / rho0 - Vcj);
            if (tmp > W[j] && tmp < 0.9) {
                new_W[j] = 0;
            }
            else {
                new_W[j] = tmp;
            }
        }

        for (int j = 1; j < Nx - 1; j++) {
            if (new_W[j] < 0.99) {
                new_P[j] = P0;
            }
            else {
                new_P[j] = (1 - new_W[j]) * (gamma - 1) * new_rho[j] * (new_I[j] + Q);
            }
        }

        for (int j = 1; j < Nx; j++) {
            u[j] = new_u[j];
        }
        for (int j = 1; j < Nx - 1; j++) {
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
            writeCSV(title, center_grid, u, P, rho, t, 1);
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
    writeCSV(title, center_grid, u, P, rho, t, 1);
}


int main() {
    krest_detonate();
}