#include <iostream>
#include<fstream>
#include <string>
#include <vector>
#include <functional>
#include <math.h>


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
        tmp += std::to_string(u[i]) + '\n';
        outfile << tmp;
    }
    outfile.close();
}


double sound_speed(double p, double rho, double gamma) {
    return std::sqrt(gamma * p / rho);
}


void Test1(double* rhoL, double* uL, double* PL, double* rhoR, double* uR, double* PR) {
    *rhoL = 1;
    *rhoR = 0.125;
    *PL = 1;
    *PR = 0.1;
    *uL = 0;
    *uR = 0;
}
void AnalitSolveTest1(double x, double t, double* p, double* u, double* rho) {
    double rhoL, rhoR, PL, PR, uL, uR;
    Test1(&rhoL, &uL, &PL, &rhoR, &uR, &PR);
    double gamma = 1.4;

    double cL = sound_speed(PL, rhoL, gamma);
    double cR = sound_speed(PR, rhoR, gamma);

    double p3 = (PL + PR) / 2.;

    double tmp_PR = PR;
    double tmp_PL = PL;

    for (int i = 0; i < 35; i++) {
        double func = p3 - PL * std::pow((1 - (gamma - 1) / 2. / sound_speed(PL, rhoL, gamma) * (std::sqrt(2 * (p3 - PR) * (p3 - PR) / rhoR / ((gamma - 1) * PR + (gamma + 1) * p3)))), 2 * gamma / (gamma - 1));
        if (func > 0) {
            tmp_PL = p3;
        }
        else {
            tmp_PR = p3;
        }
        p3 = (tmp_PL + tmp_PR) / 2;
    }

    double u3 = std::sqrt(2 * (p3 - PR) * (p3 - PR) / rhoR / ((gamma - 1) * PR + (gamma + 1) * p3));
    double rho3 = rhoR * ((gamma - 1) * PR + (gamma + 1) * p3) / ((gamma + 1) * PR + (gamma - 1) * p3);
    double rho2 = rhoL * std::pow((p3 / PL), 1. / gamma);
    double p2 = p3;
    double u2 = u3;

    double u1 = 2 / (gamma + 1) * (cL + (x - 0.5) / t);
    double rho1 = rhoL * std::pow((1 - (gamma - 1) / 2. / cL * u1), 2. / (gamma - 1));
    double p1 = PL * std::pow((1 - (gamma - 1) / 2. / cL * u1), 2. * gamma / (gamma - 1));

    double XRW1 = 0.5 - cL * t;
    double XI = 0.5 + u3 * t;
    double XSW = 0.5 + (p3 - PR) / rhoR / u3 * t;
    double XRW2 = 0.5 - t * (cL - (gamma + 1) * u3 / 2.);


    if (x < XRW1 && x >= 0) {
        *p = PL;
        *u = uL;
        *rho = rhoL;
    }
    else if (x > XRW1 && x < XRW2) {
        *p = p1;
        *u = u1;
        *rho = rho1;
    }
    else if (x > XRW2 && x < XI) {
        *p = p2;
        *u = u2;
        *rho = rho2;
    }
    else if (x < XSW && x > XI) {
        *p = p3;
        *u = u3;
        *rho = rho3;
    }
    else if (x > XSW && x < 1) {
        *p = PR;
        *u = uR;
        *rho = rhoR;
    }

}
void Test2(double* rhoL, double* uL, double* PL, double* rhoR, double* uR, double* PR) {
    *rhoL = 1;
    *rhoR = 1;
    *PL = 0.4;
    *PR = 0.4;
    *uL = -2;
    *uR = 2;
}
void Test3(double* rhoL, double* uL, double* PL, double* rhoR, double* uR, double* PR) {
    *rhoL = 1;
    *rhoR = 1;
    *PL = 1000;
    *PR = 0.01;
    *uL = 0;
    *uR = 0;
}
void Test4(double* rhoL, double* uL, double* PL, double* rhoR, double* uR, double* PR) {
    *rhoL = 1;
    *rhoR = 1;
    *PL = 0.01;
    *PR = 100;
    *uL = 0;
    *uR = 0;
}
void Test5(double* rhoL, double* uL, double* PL, double* rhoR, double* uR, double* PR) {
    *rhoL = 5.99924;
    *rhoR = 5.99242;
    *PL = 460.894;
    *PR = 46.095;
    *uL = 19.5975;
    *uR = -6.19633;
}


struct InitialValues {
    double x_start = 0;
    double x_end = 1;
    int Nx = 500;
    double h = (x_end - x_start) / Nx;

    double t_start = 0;
    double t_end = 0.2;
    int max_Nt = 1000;

    double gamma = 1.4;
    double CFL = 0.5;

    double u_left, u_right, P_left, P_right, rho_left, rho_right;
};


void InitParams(InitialValues IS, int test_num, std::vector<double>* grid, std::vector<double>* rho, std::vector<double>* u, std::vector<double>* P,
    std::vector<double>* mass, std::vector<double>* Imp, std::vector<double>* rhoe, std::vector<double>* i) 
{
    for (int j = 0; j < IS.Nx + 1; j++) {
        (*grid)[j] = IS.x_start + IS.h * j;
    }
    if (test_num == 1) {
        Test1(&IS.rho_left, &IS.u_left, &IS.P_left, &IS.rho_right, &IS.u_right, &IS.P_right);
    }
    else if (test_num == 2) {
        Test2(&IS.rho_left, &IS.u_left, &IS.P_left, &IS.rho_right, &IS.u_right, &IS.P_right);
    }
    else if (test_num == 3) {
        Test3(&IS.rho_left, &IS.u_left, &IS.P_left, &IS.rho_right, &IS.u_right, &IS.P_right);
    }
    else if (test_num == 4) {
        Test4(&IS.rho_left, &IS.u_left, &IS.P_left, &IS.rho_right, &IS.u_right, &IS.P_right);
    }
    else if (test_num == 5) {
        Test5(&IS.rho_left, &IS.u_left, &IS.P_left, &IS.rho_right, &IS.u_right, &IS.P_right);
    }
    else {
        std::cerr << "Error, wrong test!";
        exit(0);
    }

    for (int j = 0; j < IS.Nx; j++) {
        double tmp = IS.x_start + IS.h * j;
        if (tmp > (IS.x_end + IS.x_start) / 2) {
            (*u)[j] = (IS.u_right);
            (*P)[j] = (IS.P_right);
            (*rho)[j] = (IS.rho_right);
        }
        else {
            (*u)[j] = (IS.u_left);
            (*P)[j] = (IS.P_left);
            (*rho)[j] = (IS.rho_left);
        }
    }

    for (int j = 0; j < IS.Nx; j++) {
        (*i)[j] = ((*P)[j] / (*rho)[j] / (IS.gamma - 1));
        (*rhoe)[j] = (*rho)[j] * ((*i)[j] + (*u)[j] * (*u)[j] / 2.);
        (*mass)[j] = (IS.h * (*rho)[j]);
        (*Imp)[j] = (*rho)[j] * (*u)[j];
    }
}


void Bound(InitialValues& IS, int fict, std::vector<double>* rho, std::vector<double>* Imp, std::vector<double>* rhoe) {
    for (int i = 1; i <= fict; ++i) {
        (*rho)[fict - i] = (*rho)[fict + i - 1];
        (*Imp)[fict - i] = (*Imp)[fict + i - 1];
        (*rhoe)[fict - i] = (*rhoe)[fict + i - 1];
        (*rho)[IS.Nx - fict + i - 1] = (*rho)[IS.Nx - fict - i];
        (*Imp)[IS.Nx - fict + i - 1] = (*Imp)[IS.Nx - fict - i];
        (*rhoe)[IS.Nx - fict + i - 1] = (*rhoe)[IS.Nx - fict - i];
    }
}


void Laks_Friedrichs(InitialValues IS, double A, std::vector<double> u, 
    std::vector<double> flux, std::vector<double>* flux_plus, std::vector<double>* flux_minus) {
    for (int j = 0; j < IS.Nx; j++) {
        (*flux_plus)[j] = 0.5 * (flux[j] + A * u[j]);
        (*flux_minus)[j] = 0.5 * (flux[j] - A * u[j]);
    }
}


void WENO_reconstruction_plus(InitialValues IS, int fict, std::vector<double> flux, std::vector<double>* new_flux) {
    double gammaLeft[3] = { 0.3, 0.6, 0.1 };

    double b[3];
    double P[3];
    double alphas[3];
    double weights[3];
    double eps = 1.e-20;

    for (int j = fict; j < IS.Nx - fict; j++) {
        double um2 = flux[j - 2];
        double um1 = flux[j - 1];
        double u = flux[j];
        double up1 = flux[j + 1];
        double up2 = flux[j + 2];

        b[2] = 13.0 / 12.0 * std::pow((um2 - 2 * um1 + u), 2)
            + 1.0 / 4.0 * std::pow((um2 - 4 * um1 + 3 * u), 2);
        b[1] = 13.0 / 12.0 * std::pow((um1 - 2 * u + up1), 2)
            + 1.0 / 4.0 * std::pow((um1 - up1), 2);
        b[0] = 13.0 / 12.0 * std::pow((u - 2 * up1 + up2), 2)
            + 1.0 / 4.0 * std::pow((3 * u - 4 * up1 + up2), 2);

        P[0] = 1.0 / 3.0 * u + 5.0 / 6.0 * up1 - 1.0 / 6.0 * up2;
        P[1] = -1.0 / 6.0 * um1 + 5.0 / 6.0 * u + 1.0 / 3.0 * up1;
        P[2] = 1.0 / 3.0 * um2 - 7.0 / 6.0 * um1 + 11.0 / 6.0 * u;

        double sum_alpha = 0;
        for (int i = 0; i < 3; i++) {
            alphas[i] = gammaLeft[i] / std::pow((eps + b[i]), 2);
            sum_alpha += alphas[i];
        }

        double tmp = 0;
        for (int i = 0; i < 3; i++) {
            weights[i] = alphas[i] / sum_alpha;
            tmp += weights[i] * P[i];
        }

        (*new_flux)[j] = tmp;
    }
}


void WENO_reconstruction_minus(InitialValues IS, int fict, std::vector<double> flux, std::vector<double>* new_flux) {
    double gammaRight[3] = { 0.1, 0.6, 0.3 };

    double b[3];
    double P[3];
    double alphas[3];
    double weights[3];
    double eps = 1.e-20;

    for (int j = fict; j < IS.Nx - fict; j++) {
        double um2 = flux[j - 2];
        double um1 = flux[j - 1];
        double u = flux[j];
        double up1 = flux[j + 1];
        double up2 = flux[j + 2];

        /* индикаторы гладкости */
        b[2] = 13.0 / 12.0 * std::pow((um2 - 2 * um1 + u), 2)
            + 1.0 / 4.0 * std::pow((um2 - 4 * um1 + 3 * u), 2);
        b[1] = 13.0 / 12.0 * std::pow((um1 - 2 * u + up1), 2)
            + 1.0 / 4.0 * std::pow((um1 - up1), 2);
        b[0] = 13.0 / 12.0 * std::pow((u - 2 * up1 + up2), 2)
            + 1.0 / 4.0 * std::pow((3 * u - 4 * up1 + up2), 2);

        /* шаблоны */
        P[0] = 11.0 / 6.0 * u - 7.0 / 6.0 * up1 + 1.0 / 3.0 * up2;
        P[1] = 1.0 / 3.0 * um1 + 5.0 / 6.0 * u - 1.0 / 6.0 * up1;
        P[2] = -1.0 / 6.0 * um2 + 5.0 / 6.0 * um1 + 1.0 / 3.0 * u;

        double sum_alpha = 0;
        for (int i = 0; i < 3; i++) {
            alphas[i] = gammaRight[i] / std::pow((eps + b[i]), 2);
            sum_alpha += alphas[i];
        }

        double tmp = 0;
        for (int i = 0; i < 3; i++) {
            weights[i] = alphas[i] / sum_alpha;
            tmp += weights[i] * P[i];
        }

        (*new_flux)[j] = tmp;
    }
}


void WENO4(InitialValues IS, int fict, double h, double tau, std::vector<std::vector<double>> U, std::vector<std::vector<double>>* L) {
    std::vector<double> flux_rho(IS.Nx), flux_Imp(IS.Nx), flux_e(IS.Nx);
    std::vector<double> u(IS.Nx), P(IS.Nx), rho(IS.Nx), mass(IS.Nx), rhoe(IS.Nx), Imp(IS.Nx), i(IS.Nx);

    for (int j = 0; j < IS.Nx; j++) {
        rho[j] = U[0][j];
        Imp[j] = U[1][j];
        rhoe[j] = U[2][j];
        mass[j] = rho[j] * h;
        u[j] = Imp[j] / rho[j];
        i[j] = rhoe[j] / rho[j] - u[j] * u[j] / 2.;
        P[j] = i[j] * rho[j] * (IS.gamma - 1);
    }

    for (int j = 0; j < IS.Nx; j++) {
        flux_rho[j] = u[j] * rho[j];
        flux_Imp[j] = u[j] * u[j] * rho[j] + P[j];
        flux_e[j] = (rhoe[j] + P[j]) * u[j];
    }

    double A_tmp = 0.;
    double A = -100000;
    for (int j = fict; j < IS.Nx - fict; j++) {
        A_tmp = std::abs(u[j]) + sound_speed(P[j], rho[j], IS.gamma);
        if (A_tmp > A) {
            A = A_tmp;
        }
    }

    std::vector<double> flux_rho_plus(IS.Nx), flux_rho_minus(IS.Nx), flux_Imp_plus(IS.Nx),
        flux_Imp_minus(IS.Nx), flux_e_plus(IS.Nx), flux_e_minus(IS.Nx);
    Laks_Friedrichs(IS, A, rho, flux_rho, &flux_rho_plus, &flux_rho_minus);
    Laks_Friedrichs(IS, A, Imp, flux_Imp, &flux_Imp_plus, &flux_Imp_minus);
    Laks_Friedrichs(IS, A, rhoe, flux_e, &flux_e_plus, &flux_e_minus);

    std::vector<double> flux_rhoR(IS.Nx), flux_rhoL(IS.Nx), flux_ImpR(IS.Nx),
        flux_ImpL(IS.Nx), flux_eR(IS.Nx), flux_eL(IS.Nx);

    WENO_reconstruction_plus(IS, fict, flux_rho_plus, &flux_rhoR);
    WENO_reconstruction_plus(IS, fict, flux_Imp_plus, &flux_ImpR);
    WENO_reconstruction_plus(IS, fict, flux_e_plus, &flux_eR);
    WENO_reconstruction_minus(IS, fict, flux_rho_minus, &flux_rhoL);
    WENO_reconstruction_minus(IS, fict, flux_Imp_minus, &flux_ImpL);
    WENO_reconstruction_minus(IS, fict, flux_e_minus, &flux_eL);

    Bound(IS, fict, &flux_rhoR, &flux_ImpR, &flux_eR);
    Bound(IS, fict, &flux_rhoL, &flux_ImpL, &flux_eL);

    for (int j = fict - 1; j < IS.Nx - fict + 1; j++) {
        (*L)[0][j] = ((flux_rhoR[j] + flux_rhoL[j + 1]) - (flux_rhoR[j - 1] + flux_rhoL[j])) / h;
        (*L)[1][j] = ((flux_ImpR[j] + flux_ImpL[j + 1]) - (flux_ImpR[j - 1] + flux_ImpL[j])) / h;
        (*L)[2][j] = ((flux_eR[j] + flux_eL[j + 1]) - (flux_eR[j - 1] + flux_eL[j])) / h;
    }
}


void WENO(int fict_num) {
    InitialValues IS;
    double tau;
    double h = (IS.x_end - IS.x_start) / IS.Nx;

    int fict = fict_num;

    // initing
    std::vector<double> grid(IS.Nx + 1), u(IS.Nx), P(IS.Nx), rho(IS.Nx), mass(IS.Nx), e(IS.Nx), Imp(IS.Nx), i(IS.Nx);
    InitParams(IS, 1, &grid, &rho, &u, &P, &mass, &Imp, &e, &i);

    std::string outname = "CSVs\\Initial_state.csv";
    writeCSV(outname, grid, u, P, rho, IS.t_start);

    double t = IS.t_start;
    int iter = 0;

    int RK = 4;

    std::vector<std::vector<std::vector<double>>> U(RK, std::vector<std::vector<double>>(3, std::vector<double>(IS.Nx)));
    std::vector<std::vector<double>> rk = { {1, 0, 0, 0},
                                            {1, 0, 0, 1},
                                            {0.75, 0.25, 0, 0.25},
                                            { 1.0 / 3.0, 0, 2.0 / 3.0, 2.0 / 3.0} };

    while (t < IS.t_end && iter < IS.max_Nt) {
        tau = (IS.t_end - IS.t_start) / IS.max_Nt;
        for (int j = fict; j < IS.Nx - fict; j++) {
            double sz = sound_speed(P[j], rho[j], IS.gamma);
            tau = std::min(IS.CFL * h / (std::abs(u[j]) + sz), tau);
        }

        for (int j = 0; j < IS.Nx; j++) {
            U[0][0][j] = rho[j];
            U[0][1][j] = Imp[j];
            U[0][2][j] = e[j];
        }

        for (int s = 1; s < RK; s++) {
            std::vector<std::vector<double>> L(3, std::vector<double>(IS.Nx));
            WENO4(IS, fict, h, tau, U[s - 1], &L);
            for (int k = 0; k < 3; k++) {
                for (int j = fict; j < IS.Nx - fict; j++) {
                    U[s][k][j] = rk[s][0] * U[0][k][j] + rk[s][1] * U[1][k][j] + rk[s][2] * U[2][k][j] - rk[s][3] * tau * L[k][j];
                }
            }
            Bound(IS, fict, &U[s][0], &U[s][1], &U[s][2]);
        }

        for (int j = 0; j < IS.Nx; j++) {
            rho[j] = U[3][0][j];
            Imp[j] = U[3][1][j];
            e[j] = U[3][2][j];
            mass[j] = rho[j] * h;
            u[j] = Imp[j] / rho[j];
            i[j] = e[j] / rho[j] - u[j] * u[j] / 2.;
            P[j] = i[j] * rho[j] * (IS.gamma - 1);
        }
        iter++;
        t += tau;
    }
    if (t >= IS.t_end) {
        std::cout << "Solve stopped by time";
    }
    else {
        std::cout << "Solve stopped by iterations number";
    }
    outname = "CSVs\\last_step_WENO.csv";
    writeCSV(outname, grid, u, P, rho, t);

    // theoretic graphs
    std::vector<double> theoretic_P(IS.Nx), theoretic_rho(IS.Nx), theoretic_u(IS.Nx);
    for (int j = 0; j < IS.Nx; j++) {
        AnalitSolveTest1(j * h, t, &theoretic_P[j], &theoretic_u[j], &theoretic_rho[j]);
    }
    writeCSV("CSVs\\Theory_graph.csv", grid, theoretic_u, theoretic_P, theoretic_rho, t);
}

int main() {
    WENO(2);
}
