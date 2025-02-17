#include "Tests.cpp"


void Artificial_viscosity(int Nx, double* rho, double* vb, double* dP, int state) { //s - mass
    const double nu_0 = 1.0;
    const double mu_0 = 2.0;
    double fict = 1;
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
                dP[i] = -mu_0 * rho[i] * std::abs(vb[i + 1] - vb[i]) * (vb[i + 1] - vb[i]);
            }
            else {
                dP[i] = 0.0;
            }
        }
    }
}


void krest() {
    int Nx = 400;
    int Nt = 10000;
    double x_start = 0;
    double x_end = 1;
    double t_start = 0;
    double t_end = 0.2;
    double tau;
    double h_start = (x_end - x_start) / Nx;
    double gamma = 1.4;

    double CFL = 0.8;


    double u_left, u_right, P_left, P_right, rho_left, rho_right;
    Test1(&rho_left, &u_left, &P_left, &rho_right, &u_right, &P_right);
    

    std::vector<double> grid(Nx + 1);
    std::vector<double> center_grid(Nx);
    for (int j = 0; j < Nx + 1; j++) {
        grid[j] = (x_start + h_start * j);
        if (j != Nx) {
            center_grid[j] = (x_start + h_start * j + h_start / 2);
        }
    }

    // theoretic graphs
    std::vector<double> theoretic_P(Nx), theoretic_rho(Nx), theoretic_u(Nx);
    for (int j = 0; j < Nx; j++) {
        AnalitSolveTest1(center_grid[j], t_end, &theoretic_P[j], &theoretic_u[j], &theoretic_rho[j]);
    }
    writeCSV("CSVs\\Theory_graph.csv", center_grid, theoretic_u, theoretic_P, theoretic_rho, t_end);

    std::vector<double> u(Nx + 1), P(Nx), rho(Nx), mass(Nx), I(Nx), dP_vis(Nx);

    //initing
    for (int j = 0; j < Nx + 1; j++) {
        double tmp = x_start + h_start * j;
        if (tmp > (x_end + x_start) / 2) {
            u[j] = (u_right);
            if (j != Nx) {
                P[j] = (P_right);
                rho[j] = (rho_right);
            }
        }
        else {
            u[j] = (u_left);
            P[j] = (P_left);
            rho[j] = (rho_left);
        }
    }

    for (int j = 0; j < Nx; j++) {
        I[j] = (P[j] / rho[j] / (gamma - 1));
        mass[j] = (h_start * rho[j]);
    }

    writeCSV("CSVs\\first_step_Krest.csv", center_grid, u, P, rho, t_start);

    double t = 0;
    int iter = 0;

    while (t < t_end && iter < Nt) {
        tau = (t_end - t_start) / Nt;
        for (int j = 0; j < Nx; j++) {
            double sz = sound_speed(P[j], rho[j], gamma);
            tau = std::min(tau, CFL * (grid[j + 1] - grid[j]) / (std::abs(u[j]) + sz));
        }

        std::vector<double> new_u(Nx + 1), new_P(Nx), new_rho(Nx), new_mass(Nx), new_I(Nx);

        Artificial_viscosity(Nx, rho.data(), u.data(), dP_vis.data(), 1); // 1 - linear, 2 - Latter
        for (int j = 1; j < Nx - 1; j++) {
            P[j] += dP_vis[j];
        }

        for (int j = 1; j < Nx; j++) {
            new_u[j] = ((-tau / mass[j] * (P[j] - P[j - 1])) + u[j]);
        }

        for (int j = 1; j < Nx; j++) {
            grid[j] = new_u[j] * tau + grid[j];
        }

        for (int j = 0; j < Nx; j++) {
            center_grid[j] = (grid[j] + grid[j + 1]) / 2.;
        }

        for (int j = 1; j < Nx - 1; j++) {
            new_rho[j] = mass[j] / (grid[j + 1] - grid[j]);
        }
        for (int j = 1; j < Nx - 1; j++) {
            double tmp;
            tmp = 1 + tau / mass[j] * (gamma - 1) * new_rho[j] * (new_u[j + 1] - new_u[j]);
            new_I[j] = (I[j] / tmp);
        }

        for (int j = 1; j < Nx - 1; j++) {
            new_P[j] = (new_rho[j] * (gamma - 1) * new_I[j]);
        }


        for (int j = 1; j < Nx; j++) {
            u[j] = new_u[j];
        }
        for (int j = 1; j < Nx - 1; j++) {
            rho[j] = new_rho[j];
            P[j] = new_P[j];
            I[j] = new_I[j];
        }

        t += tau;
        iter++;
    }
    if (t >= t_end) {
        std::cout << "Solve stopped by time";
    }
    else {
        std::cout << "Solve stopped by iterations number";
    }
    writeCSV("CSVs\\last_step_Krest.csv", center_grid, u, P, rho, t);
}

int main() {
    krest();
}