#include <iostream>
#include<fstream>
#include <string>
#include <vector>
#include <functional>


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
    int fict;
};


struct Sol {
    std::vector<double> grid;
    std::vector<double> u, P, rho, mass, rhoe, Imp;
};


void writeCSV_new(std::string name, Sol sol, double t) {
    std::ofstream outfile(name);
    std::string tmp;
    tmp = "Time,X,Rho,P,U,E\n";
    outfile << tmp;
    for (int i = 0; i < sol.grid.size() - 1; i++) {
        tmp = "";
        tmp += std::to_string(t) + ',';
        tmp += std::to_string((sol.grid[i] + sol.grid[i + 1]) / 2.) + ',';
        tmp += std::to_string(sol.rho[i]) + ',';
        tmp += std::to_string(sol.P[i]) + ',';
        tmp += std::to_string(sol.u[i]) + ',';
        tmp += std::to_string(sol.rhoe[i] * (sol.grid[i + 1] - sol.grid[i])) + '\n';
        outfile << tmp;
    }
    outfile.close();
}


void InitParams(InitialValues IS, Sol* sol) {
    for (int j = 0; j < IS.Nx + 1; j++) {
        sol->grid.push_back(IS.x_start + IS.h * j);
    }

    Test1(&IS.rho_left, &IS.u_left, &IS.P_left, &IS.rho_right, &IS.u_right, &IS.P_right);

    for (int j = 0; j < IS.Nx; j++) {
        double tmp = IS.x_start + IS.h * j;
        if (tmp > (IS.x_end + IS.x_start) / 2) {
            sol->u.push_back(IS.u_right);
            sol->P.push_back(IS.P_right);
            sol->rho.push_back(IS.rho_right);
        }
        else {
            sol->u.push_back(IS.u_left);
            sol->P.push_back(IS.P_left);
            sol->rho.push_back(IS.rho_left);
        }
    }

    for (int j = 0; j < IS.Nx; j++) {
        sol->rhoe.push_back(sol->rho[j] * (((sol->P[j] / sol->rho[j] / (IS.gamma - 1)) + sol->u[j] * sol->u[j] / 2.)));
        sol->mass.push_back(IS.h * sol->rho[j]);
        sol->Imp.push_back(sol->rho[j] * sol->u[j]);
    }
}


void Bound(InitialValues& IS, std::vector<double>* rho, std::vector<double>* Imp, std::vector<double>* rhoe) {
    for (int i = 1; i <= IS.fict; ++i) {
        (*rho)[IS.fict - i] = (*rho)[IS.fict + i - 1];
        (*Imp)[IS.fict - i] = (*Imp)[IS.fict + i - 1];
        (*rhoe)[IS.fict - i] = (*rhoe)[IS.fict + i - 1];
        (*rho)[IS.Nx - IS.fict + i - 1] = (*rho)[IS.Nx - IS.fict - i];
        (*Imp)[IS.Nx - IS.fict + i - 1] = (*Imp)[IS.Nx - IS.fict - i];
        (*rhoe)[IS.Nx - IS.fict + i - 1] = (*rhoe)[IS.Nx - IS.fict - i];
    }
}


void cons_to_noncons(InitialValues IS, Sol* sol) {
    for (int j = 0; j < sol->grid.size() - 1; j++) {
        sol->u[j] = sol->Imp[j] / sol->rho[j];
        sol->mass[j] = sol->rho[j] * (sol->grid[j + 1] - sol->grid[j]);
        sol->P[j] = (sol->rhoe[j] / sol->rho[j] - sol->u[j] * sol->u[j] / 2.) * (IS.gamma - 1) * sol->rho[j];
    }
}


void cons_to_noncons(InitialValues IS, std::vector<double> rho, std::vector<double> Imp, 
    std::vector<double> rhoe, std::vector<double>* u, std::vector<double>* mass, std::vector<double>* P) {
    for (int j = 0; j < IS.Nx; j++) {
        (*u)[j] = Imp[j] / rho[j];
        (*mass)[j] = rho[j] * IS.h;
        (*P)[j] = (rhoe[j] / rho[j] - (*u)[j] * (*u)[j] / 2.) * (IS.gamma - 1) * rho[j];
    }
}


void noncons_to_cons(InitialValues IS, Sol* sol) {
    for (int j = 0; j < sol->grid.size() - 1; j++) {
        sol->rho[j] = sol->mass[j] / (sol->grid[j + 1] - sol->grid[j]);
        sol->Imp[j] = sol->u[j] * sol->rho[j];
        sol->rhoe[j] = (sol->rho[j] * (((sol->P[j] / sol->rho[j] / (IS.gamma - 1)) + sol->u[j] * sol->u[j] / 2.)));
    }
}


void noncons_to_cons(InitialValues IS, std::vector<double>* rho, std::vector<double>* Imp,
    std::vector<double>* rhoe, std::vector<double> u, std::vector<double> mass, std::vector<double> P) {
    for (int j = 0; j < IS.Nx; j++) {
        (*rho)[j] = mass[j] / IS.h;
        (*Imp)[j] = u[j] * (*rho)[j];
        (*rhoe)[j] = ((*rho)[j] * (((P[j] / (*rho)[j] / (IS.gamma - 1)) + u[j] * u[j] / 2.)));
    }
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


void Artificial_viscosity_linear(int Nx, double* rho, double* w, double* vb, double* p) { //s - mass
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


void krest() {
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

    writeCSV("CSVs\\first_step_Krest.csv", grid_P_rho, u, P, rho, t_start);

    for (int n = 0; n < Nt; n++) {

        Artificial_viscosity_Latter(Nx, rho.data(), w.data(), u.data(), P.data());
        Artificial_viscosity_linear(Nx, rho.data(), w.data(), u.data(), P.data());
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
    writeCSV("CSVs\\last_step_Krest.csv", grid_P_rho, u, P, rho, t_end);
}


double shock_wave_func(double analit_p, double p, double rho, double gamma) {
    return (analit_p - p) / std::sqrt(rho / 2 * ((gamma + 1) * analit_p + (gamma - 1) * p));
}


double der_shock_wave_func(double analit_p, double p, double rho, double gamma) {
    return 1 / std::sqrt(2 * rho) * ((gamma + 1) * analit_p + (3 * gamma - 1) * p) / pow(((gamma + 1) * analit_p + (gamma - 1) * p), 1.5);
}


double rare_wave_func(double analit_p, double p, double rho, double gamma) {
    return 2. / (gamma - 1) * sound_speed(p, rho, gamma) * (std::pow((analit_p / p), (gamma - 1.) / 2. / gamma) - 1.);
}


double der_rare_wave_func(double analit_p, double p, double rho, double gamma) {
    double FD = sound_speed(p, rho, gamma) * gamma / p * std::pow((analit_p / p), (gamma - 1) / gamma / 2. - 1);
    return FD;
}


double shock_wave_rho(double analit_p, double p, double rho, double gamma) {
    return rho * ((gamma + 1) * analit_p + (gamma - 1) * p) / ((gamma - 1) * analit_p + (gamma + 1) * p);
}


double rare_wave_rho(double analit_p, double p, double rho, double gamma) {
    return rho * std::pow((analit_p / p), 1 / gamma);
}


void Newton_find_P(std::function<double(double, double, double, double)> fL,
    std::function<double(double, double, double, double)> fR,
    std::function<double(double, double, double, double)> dfL,
    std::function<double(double, double, double, double)> dfR,
    double uL, double uR, double rhoL, double rhoR, double gamma, double pL, double pR, double* p) {
    double add = (fL(*p, pL, rhoL, gamma) + fR(*p, pR, rhoR, gamma) - (uL - uR)) / (dfL(*p, pL, rhoL, gamma) + dfR(*p, pR, rhoR, gamma));
    for (int i = 0; i < 30; i++) {
        (*p) -= add;
        add = (fL(*p, pL, rhoL, gamma) + fR(*p, pR, rhoR, gamma) - (uL - uR)) / (dfL(*p, pL, rhoL, gamma) + dfR(*p, pR, rhoR, gamma));
    }
}


void find_in_which_part(double* PM, double* UM, double* rhoM, double PL, double PR, double uL, double uR, double rhoL, double rhoR, double gamma, double S) {
    if (S <= *UM) {
        if (*PM <= PL) {
            // слева волна разрежения 
            if (S <= (uL - sound_speed(PL, rhoL, gamma))) {
                // слева вне распада разрыва
                *rhoM = rhoL;
                *UM = uL;
                *PM = PL;
            }
            else {
                if (S > (*UM - sound_speed(PL, rhoL, gamma) * std::pow((*PM / PL), (gamma - 1) / 2. / gamma))) {
                    // слева от контактного разрыва
                    *rhoM = rhoL * std::pow((*PM / PL), 1. / gamma);
                }
                else {
                    // слева внутри волны разрежения
                    *UM = 2. / (gamma + 1.) * (sound_speed(PL, rhoL, gamma) + (gamma - 1.) / 2. * uL + S);
                    double C = 2. / (gamma + 1.) * (sound_speed(PL, rhoL, gamma) + (gamma - 1.) / 2. * (uL - S));
                    *rhoM = rhoL * std::pow((C / sound_speed(PL, rhoL, gamma)), 2. / (gamma - 1));
                    *PM = PL * std::pow((C / sound_speed(PL, rhoL, gamma)), 2. * gamma / (gamma - 1));
                }
            }
        }
        else {
            // слева ударная волна
            if (S <= uL - sound_speed(PL, rhoL, gamma) * std::sqrt(*PM / PL * (gamma + 1) / 2. / gamma + (gamma - 1) / 2. / gamma)) {
                // слева вне распада разрыва
                *rhoM = rhoL;
                *UM = uL;
                *PM = PL;
            }
            else {
                // слева от контактного разрыва
                *rhoM = rhoL * (*PM / PL + (gamma - 1) / (gamma + 1)) / (*PM / PL * (gamma - 1) / (gamma + 1) + 1.0);
            }
        }
    }
    else {
        if (*PM > PR) {
            // справа ударная волна
            double PMR = *PM / PR;
            double SR = uR + sound_speed(PR, rhoR, gamma) * std::sqrt(PMR * (gamma + 1) / 2. / gamma + (gamma - 1) / 2. / gamma);
            if (S >= SR) {
                // справа вне распада разрыва
                *rhoM = rhoR;
                *UM = uR;
                *PM = PR;
            }
            else {
                // справа от контактного разрыва
                double G6 = (gamma - 1) / (gamma + 1);
                *rhoM = rhoR * (PMR + G6) / (PMR * G6 + 1.);
            }
        }
        else {
            // справа волна разрежения
            double SHR = uR + sound_speed(PR, rhoR, gamma);
            if (S >= SHR) {
                // справа вне распада разрыва
                *rhoM = rhoR;
                *UM = uR;
                *PM = PR;
            }
            else {
                double CMR = sound_speed(PR, rhoR, gamma) * std::pow((*PM / PR), (gamma - 1) / 2. / gamma);
                double STR = *UM + CMR;
                if (S <= STR) {
                    // справа от контактного разрыва
                    *rhoM = rhoR * std::pow((*PM / PR), 1. / gamma);
                }
                else {
                    // внутри правой волны разрежения
                    *UM = 2. / (gamma + 1.) * (-sound_speed(PR, rhoR, gamma) + (gamma - 1.) / 2. * uR + S);
                    double C = 2. / (gamma + 1.) * (sound_speed(PR, rhoR, gamma) - (gamma - 1.) / 2. * (uR - S));
                    *rhoM = rhoR * std::pow((C / sound_speed(PR, rhoR, gamma)), 2. / (gamma - 1));
                    *PM = PR * std::pow((C / sound_speed(PR, rhoR, gamma)), 2. * gamma / (gamma - 1));
                }
            }
        }
    }
}


void Rieman_solver(double PL, double uL, double rhoL, double PR, double uR, double rhoR, double* Pans, double* uans, double* rhoans, double gamma, double S) {
    // temporary
    double cR = std::sqrt(gamma * PR / rhoR);
    double cL = std::sqrt(gamma * PL / rhoL);
    (*Pans) = (PL * rhoR * cR + PR * rhoL * cL + (uL - uR) * rhoL * rhoR * cL * cR) / (rhoL * cL + rhoR * cR);
    std::function<double(double, double, double, double)> funcL, funcR, derfuncL, derfuncR;
    if (2 / (gamma - 1) * (cL + cR) <= (uR - uL)) {
        std::cout << "Vacuum!!!" << std::endl;
        exit(0);
    }
    if (*Pans < PL) {
        funcL = rare_wave_func;
        derfuncL = der_rare_wave_func;
    }
    else {
        funcL = shock_wave_func;
        derfuncL = der_shock_wave_func;
    }
    if (*Pans < PR) {
        funcR = rare_wave_func;
        derfuncR = der_rare_wave_func;
    }
    else {
        funcR = shock_wave_func;
        derfuncR = der_shock_wave_func;
    }

    Newton_find_P(funcL, funcR, derfuncL, derfuncR, uL, uR, rhoL, rhoR, gamma, PL, PR, Pans);

    (*uans) = uR + funcR(*Pans, PR, rhoR, gamma);
    S = 0;
    find_in_which_part(Pans, uans, rhoans, PL, PR, uL, uR, rhoL, rhoR, gamma, S);
}


double minmod(double a, double b, int type) {
    if (type == 1) {
        if (std::abs(a) < std::abs(b)) {
            return a;
        }
        else {
            return b;
        }
    }
    else if (type == 2) {
        double c = (a + b) / 2.;
        if (std::abs(a) < std::abs(b) && std::abs(a) < std::abs(c)) {
            return a;
        }
        else if (std::abs(b) < std::abs(a) && std::abs(b) < std::abs(c)) {
            return b;
        }
        else {
            return c;
        }
    }
    else if (type == 3) {
        if (a * a < a * b) {
            return a;
        }
        else if (b * b < a * b) {
            return b;
        }
        else if (a * b < 0) {
            return 0;
        }
    }
}


void BoundaryValuesForRiemanSolver(std::vector<double> P, std::vector<double> u, std::vector<double> rho,
    int index, double* PL, double* PR, double* uL, double* uR, double* rhoL, double* rhoR, int type) {
    if (type == 1) {
        *PL = P[index];
        *PR = P[index + 1];
        *uL = u[index];
        *uR = u[index + 1];
        *rhoL = rho[index];
        *rhoR = rho[index + 1];
    }
    else if (type == 2) {
        int minmod_type = 2;
        double kL_P = minmod((P[index] - P[index - 1]), (P[index + 1] - P[index]), minmod_type);
        *PL = P[index] + kL_P / 2;
        double kR_P = minmod((P[index + 1] - P[index]), (P[index + 2] - P[index + 1]), minmod_type);
        *PR = P[index + 1] - kR_P / 2;

        double kL_u = minmod((u[index] - u[index - 1]), (u[index + 1] - u[index]), minmod_type);
        *uL = u[index] + kL_u / 2;
        double kR_u = minmod((u[index + 1] - u[index]), (u[index + 2] - u[index + 1]), minmod_type);
        *uR = u[index + 1] - kR_u / 2;

        double kL_rho = minmod((rho[index] - rho[index - 1]), (rho[index + 1] - rho[index]), minmod_type);
        *rhoL = rho[index] + kL_rho / 2;
        double kR_rho = minmod((rho[index + 1] - rho[index]), (rho[index + 2] - rho[index + 1]), minmod_type);
        *rhoR = rho[index + 1] - kR_rho / 2;
    }
}


void Godunov_Kolgan_solver(int fict_num) {
    int Nx = 1000;
    int max_Nt = 1000;
    double x_start = 0;
    double x_end = 1;
    double t_start = 0;
    double t_end = 0.2;
    double tau;
    double h = (x_end - x_start) / Nx;
    double gamma = 1.4;

    int fict = fict_num;

    double CFL = 0.5;
    double u_left, u_right, P_left, P_right, rho_left, rho_right;

    Test1(&rho_left, &u_left, &P_left, &rho_right, &u_right, &P_right);

    std::vector<double> grid;
    for (int j = 0; j < Nx; j++) {
        grid.push_back(x_start + h * j);
    }

    std::vector<double> u, P, rho, mass, w(Nx), e;

    // theoretic graphs
    std::vector<double> theoretic_P(Nx), theoretic_rho(Nx), theoretic_u(Nx);
    for (int j = 0; j < Nx; j++) {
        AnalitSolveTest1(j * h, t_end, &theoretic_P[j], &theoretic_u[j], &theoretic_rho[j]);
    }
    writeCSV("CSVs\\Theory_graph.csv", grid, theoretic_u, theoretic_P, theoretic_rho, t_end);

    //initing
    for (int j = 0; j < Nx; j++) {
        double tmp = x_start + h * j;
        if (tmp > (x_end + x_start) / 2) {
            u.push_back(u_right);
            P.push_back(P_right);
            rho.push_back(rho_right);
        }
        else {
            u.push_back(u_left);
            P.push_back(P_left);
            rho.push_back(rho_left);
        }
    }

    for (int j = 0; j < Nx; j++) {
        e.push_back(P[j] / rho[j] / (gamma - 1));
    }

    for (int j = 0; j < Nx; j++) {
        mass.push_back(h * rho[j]);
    }
    std::string outname = "CSVs\\step0.csv";
    writeCSV(outname, grid, u, P, rho, t_start);

    double t = t_start;
    int iter = 0;
    while (t < t_end && iter < max_Nt) {
        tau = (t_end - t_start) / max_Nt;
        for (int j = 0; j < Nx; j++) {
            double sz = sqrt(gamma * P[j] / rho[j]);
            tau = std::min(CFL * h / (std::abs(u[j]) + sz), tau);
        }

        std::vector<double> new_u(Nx), new_rho(Nx), new_P(Nx), new_e(Nx), new_mass(Nx), add_mass(Nx), add_imp(Nx), add_E(Nx);
        std::vector<double> analit_P(Nx + 1), analit_u(Nx + 1), analit_rho(Nx + 1), analit_e(Nx + 1);

        for (int j = 0; j < fict - 1; j++) {
            analit_P[j] = P_left;
            analit_u[j] = u_left;
            analit_rho[j] = rho_left;
        }

        for (int j = 0; j < fict; j++) {
            analit_P[Nx - j] = P_right;
            analit_u[Nx - j] = u_right;
            analit_rho[Nx - j] = rho_right;
        }

        for (int j = fict - 1; j < Nx - fict; j++) {
            double PL, uL, rhoL, PR, uR, rhoR, Pans, uans, rhoans;

            // 1 для Годунова, 2 для Годунова-Колгана
            BoundaryValuesForRiemanSolver(P, u, rho, j, &PL, &PR, &uL, &uR, &rhoL, &rhoR, fict_num);
            double S = 0.5 * h / tau;
            Rieman_solver(PL, uL, rhoL, PR, uR, rhoR, &Pans, &uans, &rhoans, gamma, S);
            analit_P[j + 1] = Pans;
            analit_u[j + 1] = uans;
            analit_rho[j + 1] = rhoans;
            analit_e[j + 1] = Pans / rhoans / (gamma - 1) + uans * uans / 2.;
        }

        for (int j = 0; j < Nx; j++) {
            if (j > fict - 1 && j < Nx - fict) {
                add_mass[j] = (analit_rho[j] * analit_u[j] - analit_rho[j + 1] * analit_u[j + 1]) * tau;
                add_imp[j] = tau * (analit_rho[j] * analit_u[j] * analit_u[j] - analit_rho[j + 1] * analit_u[j + 1] * analit_u[j + 1] + analit_P[j] - analit_P[j + 1]);
                add_E[j] = tau * (analit_rho[j] * analit_u[j] * analit_e[j] - analit_rho[j + 1] * analit_u[j + 1] * analit_e[j + 1] + analit_P[j] * analit_u[j] - analit_P[j + 1] * analit_u[j + 1]);

            }
            else {
                add_mass[j] = 0;
                add_imp[j] = 0;
                add_E[j] = 0;
            }
        }

        for (int j = 0; j < Nx; j++) {
            new_rho[j] = rho[j] + add_mass[j] / h;
            new_mass[j] = h * new_rho[j];
            new_u[j] = u[j] * rho[j] / new_rho[j] + add_imp[j] / new_rho[j] / h;
            new_e[j] = e[j] * rho[j] / new_rho[j] + add_E[j] / new_rho[j] / h;
            new_P[j] = (new_e[j] - new_u[j] * new_u[j] / 2.) * (gamma - 1) * new_rho[j];
        }
        rho = new_rho;
        mass = new_mass;
        u = new_u;
        e = new_e;
        P = new_P;

        iter++;
        t += tau;
    }
    if (t >= t_end) {
        std::cout << "Solve stopped by time";
    }
    else {
        std::cout << "Solve stopped by iterations number";
    }
    outname = "CSVs\\last_step_Godunov.csv";
    writeCSV(outname, grid, u, P, rho, t);
}



void calculate_flux(std::vector<double> rho, std::vector<double> P, std::vector<double> u, std::vector<double> e, int index, double* fluxImp, double* fluxRho, double* fluxE) {
    (*fluxRho) = rho[index] * u[index] - rho[index - 1] * u[index - 1];
    (*fluxImp) = rho[index] * u[index] * u[index] + P[index] - rho[index - 1] * u[index - 1] * u[index - 1] - P[index - 1];
    (*fluxE) = P[index] * u[index] + u[index] * e[index] - (P[index - 1] * u[index - 1] + u[index - 1] * e[index - 1]);
}

/*
void MakKormak(int fict_num) {
    InitialValues IS;
    double tau;
    double h = (IS.x_end - IS.x_start) / IS.Nx;

    int fict = fict_num;

    std::vector<double> grid;
    for (int j = 0; j < IS.Nx; j++) {
        grid.push_back(IS.x_start + h * j);
    }

    // theoretic graphs
    std::vector<double> theoretic_P(IS.Nx), theoretic_rho(IS.Nx), theoretic_u(IS.Nx);
    for (int j = 0; j < IS.Nx; j++) {
        AnalitSolveTest1(j * h, IS.t_end, &theoretic_P[j], &theoretic_u[j], &theoretic_rho[j]);
    }
    writeCSV("CSVs\\Theory_graph.csv", grid, theoretic_u, theoretic_P, theoretic_rho, IS.t_end);

    // initing
    std::vector<double> u(IS.Nx), P(IS.Nx), rho(IS.Nx), mass(IS.Nx), e(IS.Nx), Imp(IS.Nx), i(IS.Nx);
    InitParams(IS, &rho, &u, &P, &mass, &Imp, &e, &i);
    std::vector<double> new_u(IS.Nx), new_rho(IS.Nx), new_P(IS.Nx), new_e(IS.Nx), new_mass(IS.Nx), new_Imp(IS.Nx), new_i(IS.Nx);
    
    // write initial state
    std::string outname = "CSVs\\Initial_state.csv";
    writeCSV(outname, grid, u, P, rho, IS.t_start);

    // init cycle for solver
    double t = IS.t_start;
    int iter = 0;

    // initing specially for Kormak
    std::vector<double> temp_rho(IS.Nx), temp_Imp(IS.Nx), temp_e(IS.Nx);
    int Kormak_steps = 2;

    while (t < IS.t_end && iter < IS.max_Nt) {
        tau = (IS.t_end - IS.t_start) / IS.max_Nt;
        for (int j = 0; j < IS.Nx; j++) {
            double sz = sound_speed(P[j], rho[j], IS.gamma);
            tau = std::min(IS.CFL * h / (std::abs(u[j]) + sz), tau);
        }

        new_u = u;
        new_rho = rho;
        new_P = P;
        new_e = e;
        new_mass = mass;
        new_Imp = Imp;
        new_i = i;

        for (int k = 0; k < Kormak_steps; k++) {
            for (int j = fict; j < IS.Nx - fict; j++) {
                // bad initing here, but understandable
                double fluxImp, fluxRho, fluxE;
                calculate_flux(new_rho, new_P, new_u, new_e, j + k, &fluxImp, &fluxRho, &fluxE);
                if (k == 0) {
                    temp_rho[j] = new_rho[j] - fluxRho * tau / h;
                    temp_Imp[j] = new_Imp[j] - fluxImp * tau / h;
                    temp_e[j] = new_e[j] - fluxE * tau / h;
                }
                else if (k == 1) {
                    temp_rho[j] = 0.5 * (rho[j] + new_rho[j]) - 0.5 * fluxRho * tau / h;
                    temp_Imp[j] = 0.5 * (Imp[j] + new_Imp[j]) - 0.5 * fluxImp * tau / h;
                    temp_e[j] = 0.5 * (e[j] + new_e[j]) - 0.5 * fluxE * tau / h;
                }
            }
            for (int j = fict; j < IS.Nx - fict; j++) {
                new_rho[j] = temp_rho[j];
                new_Imp[j] = temp_Imp[j];
                new_e[j] = temp_e[j];
                new_mass[j] = new_rho[j] * h;
                new_u[j] = new_Imp[j] / new_rho[j];
                new_i[j] = new_e[j] / new_rho[j] - new_u[j] * new_u[j] / 2.;
                new_P[j] = new_i[j] * new_rho[j] * (IS.gamma - 1);
            }
        }

        u = new_u;
        rho = new_rho;
        P = new_P;
        e = new_e;
        mass = new_mass;
        Imp = new_Imp;
        i = new_i;

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
*/

void Laks_Friedrichs(InitialValues IS, double A, std::vector<double> u, 
    std::vector<double> flux, std::vector<double>* flux_plus, std::vector<double>* flux_minus) {
    for (int j = 0; j < IS.Nx; j++) {
        (*flux_plus)[j] = 0.5 * (flux[j] + A * u[j]);
        (*flux_minus)[j] = 0.5 * (flux[j] - A * u[j]);
    }
}


void WENO_reconstruction_plus(InitialValues IS, std::vector<double> flux, std::vector<double>* new_flux) {
    double gammaLeft[3] = { 0.3, 0.6, 0.1 };

    double b[3];
    double P[3];
    double alphas[3];
    double weights[3];
    double eps = 1.e-20;

    for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
        double um2 = flux[j - 2];
        double um1 = flux[j - 1];
        double u = flux[j];
        double up1 = flux[j + 1];
        double up2 = flux[j + 2];

        b[2] = 13.0 / 12.0 * pow((um2 - 2 * um1 + u), 2)
            + 1.0 / 4.0 * pow((um2 - 4 * um1 + 3 * u), 2);
        b[1] = 13.0 / 12.0 * pow((um1 - 2 * u + up1), 2)
            + 1.0 / 4.0 * pow((um1 - up1), 2);
        b[0] = 13.0 / 12.0 * pow((u - 2 * up1 + up2), 2)
            + 1.0 / 4.0 * pow((3 * u - 4 * up1 + up2), 2);

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


void WENO_reconstruction_minus(InitialValues IS, std::vector<double> flux, std::vector<double>* new_flux) {
    double gammaRight[3] = { 0.1, 0.6, 0.3 };

    double b[3];
    double P[3];
    double alphas[3];
    double weights[3];
    double eps = 1.e-20;

    for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
        double um2 = flux[j - 2];
        double um1 = flux[j - 1];
        double u = flux[j];
        double up1 = flux[j + 1];
        double up2 = flux[j + 2];

        /* индикаторы гладкости */
        b[2] = 13.0 / 12.0 * pow((um2 - 2 * um1 + u), 2)
            + 1.0 / 4.0 * pow((um2 - 4 * um1 + 3 * u), 2);
        b[1] = 13.0 / 12.0 * pow((um1 - 2 * u + up1), 2)
            + 1.0 / 4.0 * pow((um1 - up1), 2);
        b[0] = 13.0 / 12.0 * pow((u - 2 * up1 + up2), 2)
            + 1.0 / 4.0 * pow((3 * u - 4 * up1 + up2), 2);

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


void WENO4(InitialValues IS, double tau, std::vector<std::vector<double>> U, std::vector<std::vector<double>>* L) {
    std::vector<double> flux_rho(IS.Nx), flux_Imp(IS.Nx), flux_e(IS.Nx);
    std::vector<double> u(IS.Nx), P(IS.Nx), rho(IS.Nx), mass(IS.Nx), rhoe(IS.Nx), Imp(IS.Nx);

    for (int j = 0; j < IS.Nx; j++) {
        rho[j] = U[0][j];
        Imp[j] = U[1][j];
        rhoe[j] = U[2][j];
        cons_to_noncons(IS, rho, Imp, rhoe, &u, &mass, &P);
    }

    for (int j = 0; j < IS.Nx; j++) {
        flux_rho[j] = u[j] * rho[j];
        flux_Imp[j] = u[j] * u[j] * rho[j] + P[j];
        flux_e[j] = (rhoe[j] + P[j]) * u[j];
    }

    double A_tmp = 0.;
    double A = -100000;
    for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
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

    WENO_reconstruction_plus(IS, flux_rho_plus, &flux_rhoR);
    WENO_reconstruction_plus(IS, flux_Imp_plus, &flux_ImpR);
    WENO_reconstruction_plus(IS, flux_e_plus, &flux_eR);
    WENO_reconstruction_minus(IS, flux_rho_minus, &flux_rhoL);
    WENO_reconstruction_minus(IS, flux_Imp_minus, &flux_ImpL);
    WENO_reconstruction_minus(IS, flux_e_minus, &flux_eL);

    Bound(IS, &flux_rhoR, &flux_ImpR, &flux_eR);
    Bound(IS, &flux_rhoL, &flux_ImpL, &flux_eL);

    for (int j = IS.fict - 1; j < IS.Nx - IS.fict + 1; j++) {
        (*L)[0][j] = ((flux_rhoR[j] + flux_rhoL[j + 1]) - (flux_rhoR[j - 1] + flux_rhoL[j]));
        (*L)[1][j] = ((flux_ImpR[j] + flux_ImpL[j + 1]) - (flux_ImpR[j - 1] + flux_ImpL[j]));
        (*L)[2][j] = ((flux_eR[j] + flux_eL[j + 1]) - (flux_eR[j - 1] + flux_eL[j]));
    }
}


void WENO(int fict_num=2) {
    InitialValues IS;
    double tau;
    IS.fict = fict_num;

    // initing
    Sol sol;
    InitParams(IS, &sol);

    std::string outname = "CSVs\\Initial_state.csv";
    writeCSV_new(outname, sol, IS.t_start);

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
        for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
            double sz = sound_speed(sol.P[j], sol.rho[j], IS.gamma);
            tau = std::min(IS.CFL * IS.h / (std::abs(sol.u[j]) + sz), tau);
        }
        //WENO method

        for (int j = 0; j < IS.Nx; j++) {
            U[0][0][j] = sol.rho[j];
            U[0][1][j] = sol.Imp[j];
            U[0][2][j] = sol.rhoe[j];
        }

        for (int s = 1; s < RK; s++) {
            std::vector<std::vector<double>> L(3, std::vector<double>(IS.Nx));
            WENO4(IS, tau, U[s - 1], &L);
            for (int k = 0; k < 3; k++) {
                for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
                    U[s][k][j] = rk[s][0] * U[0][k][j] + rk[s][1] * U[1][k][j] + rk[s][2] * U[2][k][j] - rk[s][3] * tau / IS.h * L[k][j];
                }
            }
            Bound(IS, &U[s][0], &U[s][1], &U[s][2]);
        }

        for (int j = 0; j < IS.Nx; j++) {
            sol.rho[j] = U[3][0][j];
            sol.Imp[j] = U[3][1][j];
            sol.rhoe[j] = U[3][2][j];
            cons_to_noncons(IS, &sol);
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
    writeCSV_new(outname, sol, t);

    // theoretic graphs
    std::vector<double> grid(IS.Nx), theoretic_P(IS.Nx), theoretic_rho(IS.Nx), theoretic_u(IS.Nx);
    for (int j = 0; j < IS.Nx; j++) {
        grid[j] = (IS.x_start + IS.h * j + IS.h / 2.);
        AnalitSolveTest1(j * IS.h, t, &theoretic_P[j], &theoretic_u[j], &theoretic_rho[j]);
    }
    writeCSV("CSVs\\Theory_graph.csv", grid, theoretic_u, theoretic_P, theoretic_rho, t);
}

int main() {
    // fict_cells: for Godunov 1, for Kolgan 2, for MakKormak 1, for WENO 2
    // MakKormak(1);
    WENO(2);
    //Godunov_Kolgan_solver(1);
}
