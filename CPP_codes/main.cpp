#include <iostream>
#include<fstream>
#include <string>
#include <vector>
#include <map>
#include <functional>


struct Point {
    double x, y, z;
};


double Rho(Point point) {
    if (point.x < 5) {
        return 1;
    }
    else {
        return 0;
    }
}


double U(Point point) {
    if (point.y < 2) {
        return 2;
    }
    else {
        return 0;
    }
}


double P(Point point) {
    if (point.z < 0.5) {
        return 1;
    }
    else {
        return 0;
    }
}


void writeCSV(std::string name, std::vector<double> grid_P, std::vector<double> u, std::vector<double> P, std::vector<double> rho, double t) {
    std::ofstream outfile(name);
    std::string tmp;
    tmp = "Time,X,Rho,P,U\n";
    outfile << tmp;
    for (int i = 0; i < P.size(); i++) {
        tmp = "";
        tmp += std::to_string(t) + ',';
        tmp += std::to_string(grid_P[i]) + ',';
        tmp += std::to_string(rho[i]) + ',';
        tmp += std::to_string(P[i]) + ',';
        tmp += std::to_string((u[i] + u[i + 1]) / 2) + '\n';
        outfile << tmp;
    }
    outfile.close();
}


void old_krest() {
    int Nx = 1000;
    int Nt = 500;
    double x_start = 0;
    double x_end = 1;
    double t_start = 0;
    double t_end = 0.25;
    double tau = (t_end - t_start) / Nt;
    double h = (x_end - x_start) / Nx;
    double gamma = 1.4;


    double u_left = 0;
    double u_right = 0;
    double P_left = 1;
    double P_right = 0.1;
    double rho_left = 1;
    double rho_right = 0.125;


    std::vector<double> grid_u;
    std::vector<double> grid_P_rho;
    for (int j = 0; j < Nx + 1; j++) {
        grid_u.push_back(x_start + h * j);
        if (j != Nx) {
            grid_P_rho.push_back(x_start + h * j + h / 2);
        }
    }

    std::vector<double> u, P, rho, mass;

    //initing
    for (int j = 0; j < Nx + 1; j++) {
        double tmp = x_start + h * j;
        if (tmp > (x_end + x_start) / 2) {
            u.push_back(u_right);
            if (j != Nx) {
                P.push_back(P_right);
                rho.push_back(rho_right);
            }
        }
        else {
            u.push_back(u_left);
            P.push_back(P_left);
            rho.push_back(rho_left);
        }
    }

    writeCSV("first_step.csv", grid_P_rho, u, P, rho, t_start);

    for (int n = 0; n < Nt; n++) {
        std::vector<double> E;
        for (int j = 0; j < Nx; j++) {
            E.push_back(P[j] / rho[j] / (gamma - 1));
        }

        std::vector<double> new_u;
        for (int j = 0; j < Nx + 1; j++) {
            if (j == 0) {
                new_u.push_back(u_left);
            }
            else if (j == Nx) {
                new_u.push_back(u_right);
            }
            else {
                new_u.push_back((-tau / h * (P[j] - P[j - 1])) + u[j]);
            }
        }
        u = new_u;

        std::vector<double> new_rho;
        for (int j = 0; j < Nx; j++) {
            if (j == 0) {
                new_rho.push_back(rho_left);
            }
            else {
                new_rho.push_back(1 / ((u[j + 1] - u[j]) * tau / h + 1 / rho[j]));
            }
        }
        rho = new_rho;

        std::vector<double> new_E;
        for (int j = 0; j < Nx; j++) {
            if (j == 0) {
                new_E.push_back(P_left / rho_left / (gamma - 1));
            }
            else {
                double tmp;
                tmp = 1 + tau / h * (gamma - 1) * rho[j] * (u[j + 1] - u[j]);
                new_E.push_back(E[j] / tmp);
            }
        }

        std::vector<double> new_P;
        for (int j = 0; j < Nx; j++) {
            new_P.push_back(rho[j] * (gamma - 1) * new_E[j]);
        }
        P = new_P;

    }
    writeCSV("last_step.csv", grid_P_rho, u, P, rho, t_end);
}



void Artificial_viscosity_Latter(int Nx, double* rho, double* w, double* vb, double* p) {
    const double mu_0 = 2.0;
    double fict = 1;

    for (int i = 0; i < Nx; ++i) {
        if ((vb[i + 1] - vb[i]) < 0) {
            w[i] = -mu_0 * rho[i] * std::fabs(vb[i + 1] - vb[i]) * (vb[i + 1] - vb[i]);
        }
        else {
            w[i] = 0.0;
        }
    }
    w[0] = w[1];
    w[Nx - 1] = w[Nx - 2];
    for (int i = 0; i < Nx; ++i) {
        p[i] += w[i];
    }
}


void Artificial_linear(int Nx, double* rho, double* w, double* vb, double* p) { //s - mass
    const double nu_0 = 1.0;
    double fict = 1;

    for (int i = 0; i < Nx; ++i) {
        if ((vb[i + 1] - vb[i]) < 0) {
            w[i] = -nu_0 * rho[i] * (vb[i + 1] - vb[i]);
        }
        else {
            w[i] = 0.0;
        }
    }
    w[0] = w[1];
    w[Nx - 1] = w[Nx - 2];
    for (int i = 0; i < Nx; ++i) {
        p[i] += w[i];
    }
}


void new_krest() {
    int Nx = 1000;
    int Nt = 500;
    double x_start = 0;
    double x_end = 1;
    double t_start = 0;
    double t_end = 0.25;
    double tau = (t_end - t_start) / Nt;
    double h = (x_end - x_start) / Nx;
    double gamma = 1.4;

    double CFL = 0.5;


    double u_left = 0;
    double u_right = 0;
    double P_left = 1;
    double P_right = 0.1;
    double rho_left = 1;
    double rho_right = 0.125;


    std::vector<double> grid_u;
    std::vector<double> grid_P_rho;
    for (int j = 0; j < Nx + 1; j++) {
        grid_u.push_back(x_start + h * j);
        if (j != Nx) {
            grid_P_rho.push_back(x_start + h * j + h / 2);
        }

    }

    std::vector<double> u, P, rho, mass, w(Nx);

    //initing
    for (int j = 0; j < Nx + 1; j++) {
        double tmp = x_start + h * j;
        if (tmp > (x_end + x_start) / 2) {
            u.push_back(u_right);
            if (j != Nx) {
                P.push_back(P_right);
                rho.push_back(rho_right);
            }
        }
        else {
            u.push_back(u_left);
            P.push_back(P_left);
            rho.push_back(rho_left);
        }
    }

    std::vector<double> E;
    for (int j = 0; j < Nx; j++) {
        E.push_back(P[j] / rho[j] / (gamma - 1));
    }

    for (int j = 0; j < Nx; j++) {
        mass.push_back(h * rho[j]);
    }

    writeCSV("first_step_new.csv", grid_P_rho, u, P, rho, t_start);

    for (int n = 0; n < Nt; n++) {

        Artificial_viscosity_Latter(Nx, rho.data(), w.data(), u.data(), P.data());
        Artificial_linear(Nx, rho.data(), w.data(), u.data(), P.data());
        tau = 1111;
        for (int i = 0; i < Nx; i++) {
            double sz = sqrt(P[i] / rho[i]);
            if (CFL * h / (std::fabs(u[i]) + sz) < tau)
                tau = CFL * h / (std::fabs(u[i]) + sz);
        }

        std::vector<double> new_u;
        for (int j = 0; j < Nx + 1; j++) {
            if (j == 0) {
                new_u.push_back(0); //stenka
            }
            else if (j == Nx) {
                new_u.push_back(0);
            }
            else {
                new_u.push_back((-tau / 2.0 / (mass[j] + mass[j - 1]) * (P[j] - P[j - 1])) + u[j]);
            }
        }

        u = new_u;
        for (int i = 0; i < Nx /*+1*/; i++) {
            grid_u[i] = u[i + 1] * tau + grid_u[i];
            if (i != Nx) {
                grid_P_rho[i] = (grid_u[i] + grid_u[i + 1]) / 2.0;
            }
        }

        std::vector<double> new_rho;
        for (int j = 0; j < Nx; j++) {
            if (j == 0) {
                new_rho.push_back(rho_left);
            }
            else {
                new_rho.push_back(1 / ((u[j + 1] - u[j]) * tau / mass[j] + 1 / rho[j]));
            }
        }
        rho = new_rho;

        std::vector<double> new_E;
        for (int j = 0; j < Nx; j++) {
            if (j == 0) {
                new_E.push_back(P_left / rho_left / (gamma - 1));
            }
            else {
                double tmp;
                tmp = 1 + tau / mass[j] * (gamma - 1) * rho[j] * (u[j + 1] - u[j]);
                new_E.push_back(E[j] / tmp);
            }
        }

        std::vector<double> new_P;
        for (int j = 0; j < Nx; j++) {
            new_P.push_back(rho[j] * (gamma - 1) * new_E[j]);
        }

        P = new_P;

    }
    writeCSV("last_step_new.csv", grid_P_rho, u, P, rho, t_end);
}


int main() {
    new_krest();
}
