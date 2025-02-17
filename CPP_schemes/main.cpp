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
    *rhoR = 1;
    *PL = 10;
    *PR = 1;
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
    int Nx = 400;
    double h = (x_end - x_start) / Nx;

    double t_start = 0;
    double t_end = 0.25;
    int max_Nt = 1000;

    double gamma = 1.4;
    double CFL = 0.5;

    double u_left, u_right, P_left, P_right, rho_left, rho_right;
};


struct Sol {
    std::vector<double> grid;
    std::vector<double> u, P, rho, mass, e, Imp;
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


void cons_to_noncons(InitialValues IS, std::vector<double> rho, std::vector<double> Imp,
    std::vector<double> rhoe, std::vector<double>* u, std::vector<double>* mass, std::vector<double>* P) {
    for (int j = 0; j < IS.Nx; j++) {
        (*u)[j] = Imp[j] / rho[j];
        (*mass)[j] = rho[j] * IS.h;
        (*P)[j] = (rhoe[j] / rho[j] - (*u)[j] * (*u)[j] / 2.) * (IS.gamma - 1) * rho[j];
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
    InitialValues IS;

    int fict = fict_num;
    double tau;

    std::vector<double> grid(IS.Nx + 1), u(IS.Nx), P(IS.Nx), rho(IS.Nx), mass(IS.Nx), Imp(IS.Nx), rhoe(IS.Nx), i(IS.Nx);

    InitParams(IS, 1, &grid, &rho, &u, &P, &mass, &Imp, &rhoe, &i);

    // write initial state
    std::string outname = "CSVs\\Initial_state.csv";
    writeCSV(outname, grid, u, P, rho, IS.t_start);

    // theoretic graphs
    std::vector<double> theoretic_P(IS.Nx), theoretic_rho(IS.Nx), theoretic_u(IS.Nx);
    for (int j = 0; j < IS.Nx; j++) {
        AnalitSolveTest1(j * IS.h, IS.t_end, &theoretic_P[j], &theoretic_u[j], &theoretic_rho[j]);
    }
    writeCSV("CSVs\\Theory_graph.csv", grid, theoretic_u, theoretic_P, theoretic_rho, IS.t_end);

    double t = IS.t_start;
    int iter = 0;
    while (t < IS.t_end && iter < IS.max_Nt) {
        tau = (IS.t_end - IS.t_start) / IS.max_Nt;
        for (int j = 0; j < IS.Nx; j++) {
            double sz = sqrt(IS.gamma * P[j] / rho[j]);
            tau = std::min(IS.CFL * IS.h / (std::abs(u[j]) + sz), tau);
        }

        std::vector<double> new_u(IS.Nx), new_rho(IS.Nx), new_P(IS.Nx), new_rhoe(IS.Nx), new_mass(IS.Nx), add_mass(IS.Nx), add_imp(IS.Nx), add_E(IS.Nx);
        std::vector<double> analit_P(IS.Nx + 1), analit_u(IS.Nx + 1), analit_rho(IS.Nx + 1), analit_rhoe(IS.Nx + 1);


        for (int j = fict - 1; j < IS.Nx - fict; j++) {
            double PL, uL, rhoL, PR, uR, rhoR, Pans, uans, rhoans;

            // 1 для Годунова, 2 для Годунова-Колгана
            BoundaryValuesForRiemanSolver(P, u, rho, j, &PL, &PR, &uL, &uR, &rhoL, &rhoR, fict);
            double S = 0.5 * IS.h / tau;
            Rieman_solver(PL, uL, rhoL, PR, uR, rhoR, &Pans, &uans, &rhoans, IS.gamma, S);
            analit_P[j + 1] = Pans;
            analit_u[j + 1] = uans;
            analit_rho[j + 1] = rhoans;
            analit_rhoe[j + 1] = Pans / (IS.gamma - 1) + rhoans * uans * uans / 2.;
        }

        for (int j = 0; j < fict - 1; j++) {
            analit_P[j] = IS.P_left;
            analit_u[j] = IS.u_left;
            analit_rho[j] = IS.rho_left;
        }

        for (int j = 0; j < fict; j++) {
            analit_P[IS.Nx - j] = IS.P_right;
            analit_u[IS.Nx - j] = IS.u_right;
            analit_rho[IS.Nx - j] = IS.rho_right;
        }

        for (int j = fict; j < IS.Nx - fict; j++) {
            add_mass[j] = (analit_rho[j] * analit_u[j] - analit_rho[j + 1] * analit_u[j + 1]) * tau;
            add_imp[j] = tau * (analit_rho[j] * analit_u[j] * analit_u[j] - analit_rho[j + 1] * analit_u[j + 1] * analit_u[j + 1] + analit_P[j] - analit_P[j + 1]);
            add_E[j] = tau * (analit_u[j] * (analit_rhoe[j] + analit_P[j]) - analit_u[j + 1] * (analit_rhoe[j + 1] + analit_P[j + 1]));
        }

        for (int j = 0; j < IS.Nx; j++) {
            new_rho[j] = rho[j] + add_mass[j] / IS.h;
            new_mass[j] = IS.h * new_rho[j];
            new_u[j] = u[j] * rho[j] / new_rho[j] + add_imp[j] / new_rho[j] / IS.h;
            new_rhoe[j] = rhoe[j] + add_E[j] / IS.h;
            new_P[j] = (new_rhoe[j] - new_rho[j] * new_u[j] * new_u[j] / 2.) * (IS.gamma - 1);
        }
        rho = new_rho;
        mass = new_mass;
        u = new_u;
        rhoe = new_rhoe;
        P = new_P;

        iter++;
        t += tau;
    }
    if (t >= IS.t_end) {
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
    // fict_cells: for Godunov 1, for Kolgan 2, for MakKormak 1, for WENO 2
    //MakKormak(3, 1);
    // WENO(2);
    Godunov_Kolgan_solver(1);
    //krest();
}
