#include "Tests.cpp"


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
        AnalitSolveTest1(j * h, &theoretic_P[j], &theoretic_u[j], &theoretic_rho[j]);
    }
    writeCSV("Theory_t=0_25.csv", grid, theoretic_u, theoretic_P, theoretic_rho, t_end);

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
    std::string outname = "Steps_Godunov/step0.csv";
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
    outname = "last_step_Godunov.csv";
    writeCSV(outname, grid, u, P, rho, t);
}