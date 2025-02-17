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
    *PR = 0.3;
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


void calculate_flux(std::vector<double> rho, std::vector<double> P, std::vector<double> u, std::vector<double> e, int index, double* fluxImp, double* fluxRho, double* fluxE) {
    (*fluxRho) = rho[index] * u[index] - rho[index - 1] * u[index - 1];
    (*fluxImp) = rho[index] * u[index] * u[index] + P[index] - rho[index - 1] * u[index - 1] * u[index - 1] - P[index - 1];
    (*fluxE) = (P[index] + e[index]) * u[index] - (P[index - 1] + e[index - 1]) * u[index - 1];
}


double calc_phi_c(double Q, std::vector<double> u, std::vector<double> Du, int index) {
    double left, right, middle;
    int s = (u[index + 1] - u[index] > 0) ? 1 : -1;
    middle = std::abs(Q * (Du[index + 1] - Du[index]));
    left = s * (Du[index] - Du[index - 1]);
    right = s * (Du[index + 2] - Du[index + 1]);
    double ans = s * std::max(0., std::min({ left, middle, right }));
    return ans;
}


void diffusion_antidiffusion(InitialValues IS, int fict, double Q, std::vector<double> u, std::vector<double>* new_u) {
    std::vector<double> Du(IS.Nx);
    for (int j = fict; j < IS.Nx - fict; j++) {
        Du[j] = u[j] + Q * (u[j + 1] - u[j] - (u[j] - u[j - 1]));
    }
    for (int j = fict; j < IS.Nx - fict; j++) {
        (*new_u)[j] = Du[j] - (calc_phi_c(Q, u, Du, j) - calc_phi_c(Q, u, Du, j - 1));
    }
}


void diff_adiff_Polina(InitialValues IS, int fict, double Q, std::vector<double> h, std::vector<double>* new_h) {
    double s, um, up, u, sudm, sudp, phi_p, phi_m;
    for (int i = fict; i < IS.Nx - fict; i++) {
        s = (h[i + 1] - h[i] > 0) ? 1 : -1;
        um = 0.5 * (h[i - 1] + h[i - 2]);
        u = 0.5 * (h[i] + h[i - 1]);
        up = 0.5 * (h[i + 1] + h[i]);

        sudm = u + Q * (up - 2 * u + um);
        um = 0.5 * (h[i + 1] + h[i]);
        u = 0.5 * (h[i + 2] + h[i + 1]);
        up = 0.5 * (h[i + 3] + h[i + 2]);
        sudp = u + Q * (up - 2 * u + um);
        phi_p = s * std::max(0.0, std::min({ s * sudm, std::fabs(Q * (h[i + 1] - h[i])), s * sudp }));

        s = (h[i] - h[i - 1] > 0) ? 1 : -1;

        um = 0.5 * (h[i - 2] + h[i - 3]);

        u = 0.5 * (h[i - 1] + h[i - 2]);

        up = 0.5 * (h[i] + h[i - 1]);
        sudm = u + Q * (up - 2 * u + um);
        um = 0.5 * (h[i] + h[i - 1]);
        u = 0.5 * (h[i + 1] + h[i]);
        up = 0.5 * (h[i + 2] + h[i + 1]);
        sudp = u + Q * (up - 2 * u + um);
        phi_m = s * std::max(0.0, std::min({ s * sudm, std::fabs(Q * (h[i] - h[i - 1])), s * sudp }));

        (*new_h)[i] = sudm - (phi_p - phi_m);
    }

}


void MakKormak(int fict_num, int diff) {
    InitialValues IS;
    double tau;
    double h = (IS.x_end - IS.x_start) / IS.Nx;

    int fict = fict_num;

    // initing
    std::vector<double> grid(IS.Nx + 1), u(IS.Nx), P(IS.Nx), rho(IS.Nx), mass(IS.Nx), rhoe(IS.Nx), Imp(IS.Nx), i(IS.Nx);
    InitParams(IS, 1, &grid, &rho, &u, &P, &mass, &Imp, &rhoe, &i);
    std::vector<double> new_u(IS.Nx), new_rho(IS.Nx), new_P(IS.Nx), new_rhoe(IS.Nx), new_mass(IS.Nx), new_Imp(IS.Nx), new_i(IS.Nx);
    
    // write initial state
    std::string outname = "CSVs\\Initial_state.csv";
    writeCSV(outname, grid, u, P, rho, IS.t_start);

    // theoretic graphs
    std::vector<double> theoretic_P(IS.Nx), theoretic_rho(IS.Nx), theoretic_u(IS.Nx);
    for (int j = 0; j < IS.Nx; j++) {
        AnalitSolveTest1(j * h, IS.t_end, &theoretic_P[j], &theoretic_u[j], &theoretic_rho[j]);
    }
    writeCSV("CSVs\\Theory_graph.csv", grid, theoretic_u, theoretic_P, theoretic_rho, IS.t_end);


    // init cycle for solver
    double t = IS.t_start;
    int iter = 0;

    // initing specially for Kormak
    std::vector<double> temp_rho(IS.Nx), temp_Imp(IS.Nx), temp_rhoe(IS.Nx);
    int Kormak_steps = 2;
    // for diffusion
    double Q = 0.2;

    while (t < IS.t_end && iter < IS.max_Nt) {
        tau = (IS.t_end - IS.t_start) / IS.max_Nt;
        for (int j = 0; j < IS.Nx; j++) {
            double sz = sound_speed(P[j], rho[j], IS.gamma);
            tau = std::min(IS.CFL * h / (std::abs(u[j]) + sz), tau);
        }

        new_u = u;
        new_rho = rho;
        new_P = P;
        new_rhoe = rhoe;
        new_mass = mass;
        new_Imp = Imp;
        new_i = i;

        for (int k = 0; k < Kormak_steps; k++) {
            for (int j = fict; j < IS.Nx - fict; j++) {
                // bad initing here, but understandable
                double fluxImp, fluxRho, fluxE;
                calculate_flux(new_rho, new_P, new_u, new_rhoe, j + k, &fluxImp, &fluxRho, &fluxE);
                if (k == 0) {
                    temp_rho[j] = new_rho[j] - fluxRho * tau / h;
                    temp_Imp[j] = new_Imp[j] - fluxImp * tau / h;
                    temp_rhoe[j] = new_rhoe[j] - fluxE * tau / h;
                }
                else if (k == 1) {
                    temp_rho[j] = 0.5 * (rho[j] + new_rho[j]) - 0.5 * fluxRho * tau / h;
                    temp_Imp[j] = 0.5 * (Imp[j] + new_Imp[j]) - 0.5 * fluxImp * tau / h;
                    temp_rhoe[j] = 0.5 * (rhoe[j] + new_rhoe[j]) - 0.5 * fluxE * tau / h;
                }
            }
            for (int j = fict; j < IS.Nx - fict; j++) {
                new_rho[j] = temp_rho[j];
                new_Imp[j] = temp_Imp[j];
                new_rhoe[j] = temp_rhoe[j];
                new_mass[j] = new_rho[j] * h;
                new_u[j] = new_Imp[j] / new_rho[j];
                new_i[j] = new_rhoe[j] / new_rho[j] - new_u[j] * new_u[j] / 2.;
                new_P[j] = new_i[j] * new_rho[j] * (IS.gamma - 1);
            }
        }

        if (diff == 1) {
            diffusion_antidiffusion(IS, fict, Q, new_rho, &rho);
            diffusion_antidiffusion(IS, fict, Q, new_Imp, &Imp);
            diffusion_antidiffusion(IS, fict, Q, new_rhoe, &rhoe);
            //diff_adiff_Polina(IS, fict, Q, new_rho, &rho);
            //diff_adiff_Polina(IS, fict, Q, new_Imp, &Imp);
            //diff_adiff_Polina(IS, fict, Q, new_rhoe, &rhoe);
            for (int j = 0; j < IS.Nx; j++) {
                mass[j] = rho[j] * h;
                u[j] = Imp[j] / rho[j];
                i[j] = rhoe[j] / rho[j] - u[j] * u[j] / 2.;
                P[j] = i[j] * rho[j] * (IS.gamma - 1);
            }
        }
        else {
            u = new_u;
            rho = new_rho;
            P = new_P;
            rhoe = new_rhoe;
            mass = new_mass;
            Imp = new_Imp;
            i = new_i;
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
    outname = "CSVs\\last_step_MakKormak.csv";
    writeCSV(outname, grid, u, P, rho, t);
}

int main() {
    MakKormak(2, 0);
}