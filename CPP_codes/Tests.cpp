#ifndef TESTS
#define TESTS

#include "Base.cpp"

void Test1(double* rhoL, double* uL, double* PL, double* rhoR, double* uR, double* PR) {
    *rhoL = 1;
    *rhoR = 0.125;
    *PL = 1;
    *PR = 0.2;
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

#endif