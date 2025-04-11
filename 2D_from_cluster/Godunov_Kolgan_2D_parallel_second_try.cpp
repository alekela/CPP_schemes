#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <experimental/filesystem>
#include <stdio.h>
#include <stdlib.h>
// #include <mpi.h>
#include <unordered_map>

namespace fs = std::experimental::filesystem;
typedef std::vector<double> vec1d;
typedef std::vector<vec1d> vec2d;
typedef std::vector<vec2d> vec3d;

struct InitialState {
    double gamma;
    double x_start;
    double x_end;
    double y_start;
    double y_end;
    double hx, hy;

    int nx;
    int ny;
    int Nx;
    int Ny;
    double t_end;
    double CFL;

    int b_left;
    int b_right;
    int b_up;
    int b_down;
    int initial;
    int write_interval;
    int max_iter;

    int s_type;
    int mod_type;

    int fict;
};


void change_params(InitialState& IS, int rank, int px, int py) {
    int start_x, stop_x;
    int start_y, stop_y;
    int rank_x, rank_y;
    rank_x = rank % px;
    rank_y = rank / py;
    start_x = IS.nx / px * rank_x;
    stop_x = start_x + IS.nx / px;
    start_y = IS.ny / py * rank_y;
    stop_y = start_y + IS.ny / py;
    if (IS.nx % px) {
        if (rank_x >= px - IS.nx % px) {
            stop_x += IS.nx % px - (px - 1 - rank_x);
            start_x += IS.nx % px - 1 - (px - 1 - rank_x);
        }
    }
    if (IS.ny % py) {
        if (rank_y >= py - IS.ny % py) {
            stop_y += IS.ny % py - (py - 1 - rank_y);
            start_y += IS.ny % py - 1 - (py - 1 - rank_y);
        }
    }
    IS.nx = stop_x - start_x;
    IS.ny = stop_y - start_y;
    IS.x_start = start_x * IS.hx;
    IS.x_end = stop_x * IS.hx;
    IS.y_start = start_y * IS.hy;
    IS.y_end = stop_y * IS.hy;
    IS.Nx = IS.nx + 2 * IS.fict;
    IS.Ny = IS.ny + 2 * IS.fict;

}

//constant parameters
void init_params(InitialState& IS) {
    IS.gamma = 1.4;
    IS.x_start = 0.;
    IS.x_end = 1.;
    IS.y_start = 0.;
    IS.y_end = 1.;

    IS.nx = 10;
    IS.ny = 5;
    IS.t_end = 0.15;
    IS.CFL = 0.5;
    IS.hx = (IS.x_end - IS.x_start) / IS.nx;
    IS.hy = (IS.y_end - IS.y_start) / IS.ny;

    IS.b_left = 1;
    IS.b_right = 1;
    IS.b_up = 1;
    IS.b_down = 1;
    IS.initial = 2;
    IS.write_interval = 10;
    IS.max_iter = 1000;

    IS.s_type = 3; // 1 - Godunov, 2 - Kolgan, 3 - Rodionov, 4 - WENO
    IS.mod_type = 3;
    if (IS.s_type == 1) IS.fict = 1;
    else if (IS.s_type == 2 || IS.s_type == 3) IS.fict = 2;
    else if (IS.s_type == 4) IS.fict = 4;
    else {
        IS.s_type = 1;
        IS.fict = 1;
    }

    IS.Nx = IS.nx + 2 * IS.fict;
    IS.Ny = IS.ny + 2 * IS.fict;
}


//convert conservative variables to nonconservative
void cons_to_noncons(InitialState& params, double& p, double& vx, double& vy, double& r, double& m, double& impx, double& impy, double& e) {
    p = (params.gamma - 1.0) * (e - 0.5 * (pow(impx, 2.0) + pow(impy, 2.0)) / m);
    if (p < 0.) {
        p = 1.e-6;
    }
    vx = impx / m;
    vy = impy / m;
    r = m;
}

void cons_to_noncons(InitialState& params, vec2d& p, vec2d& vx, vec2d& vy, vec2d& r, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    for (int i = 0; i < params.Ny; ++i) {
        for (int j = 0; j < params.Nx; ++j) {
            cons_to_noncons(params, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
        }
    }
}

//convert nonconservative variables to conservative
void noncons_to_cons(InitialState& params, double& p, double& vx, double& vy, double& r, double& m, double& impx, double& impy, double& e) {
    m = r;
    impx = r * vx;
    impy = r * vy;
    e = 0.5 * r * (pow(vx, 2.0) + pow(vy, 2.0)) + p / (params.gamma - 1.0);
}

//calculate internal energy
void calc_ei(InitialState& params, vec2d& p, vec2d& r, vec2d& ei) {
    for (int i = 0; i < params.Ny; ++i) {
        for (int j = 0; j < params.Nx; ++j) {
            ei[i][j] = p[i][j] / ((params.gamma - 1.) * r[i][j]);
        }
    }
}

//calculate sound velocity
double calc_c(double gamma, double p, double r) {
    return std::sqrt(gamma * p / r);
}

//build x grid
void x_grid(InitialState& params, vec1d& xc, vec1d& x) {
    for (int i = 0; i < params.Nx + 1; ++i) {
        x[i] = params.x_start + (i - params.fict) * params.hx;
    }
    for (int i = 0; i < params.Nx; ++i) {
        xc[i] = 0.5 * (x[i] + x[i + 1]);
    }
}

//build y grid
void y_grid(InitialState& params, vec1d& yc, vec1d& y) {
    for (int i = 0; i < params.Ny + 1; ++i) {
        y[i] = params.y_start + (i - params.fict) * params.hy;
    }
    for (int i = 0; i < params.Ny; ++i) {
        yc[i] = 0.5 * (y[i] + y[i + 1]);
    }
}

//boundary Y
void Boundary_y(InitialState& params, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    for (int i = 1; i <= params.fict; ++i) {
        for (int j = 0; j < params.Nx; ++j) {
            if (params.b_left == 0) {   //wall
                m[params.fict - i][j] = m[params.fict + i - 1][j];
                impx[params.fict - i][j] = -impx[params.fict + i - 1][j];
                impy[params.fict - i][j] = impy[params.fict + i - 1][j];
                e[params.fict - i][j] = e[params.fict + i - 1][j];
            }
            else if (params.b_left == 1) {   //free boundary
                m[params.fict - i][j] = m[params.fict + i - 1][j];
                impx[params.fict - i][j] = impx[params.fict + i - 1][j];
                impy[params.fict - i][j] = impy[params.fict + i - 1][j];
                e[params.fict - i][j] = e[params.fict + i - 1][j];
            }
            if (params.b_right == 0) {   //wall
                m[params.Ny - params.fict + i - 1][j] = m[params.Ny - params.fict - i][j];
                impx[params.Ny - params.fict + i - 1][j] = -impx[params.Ny - params.fict - i][j];
                impy[params.Ny - params.fict + i - 1][j] = impy[params.Ny - params.fict - i][j];
                e[params.Ny - params.fict + i - 1][j] = e[params.Ny - params.fict - i][j];
            }
            else if (params.b_right == 1) {   //free boundary
                m[params.Ny - params.fict + i - 1][j] = m[params.Ny - params.fict - i][j];
                impx[params.Ny - params.fict + i - 1][j] = impx[params.Ny - params.fict - i][j];
                impy[params.Ny - params.fict + i - 1][j] = impy[params.Ny - params.fict - i][j];
                e[params.Ny - params.fict + i - 1][j] = e[params.Ny - params.fict - i][j];
            }
        }
    }
}

//boundary X
void Boundary_x(InitialState& params, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    for (int j = 1; j <= params.fict; ++j) {
        for (int i = 0; i < params.Ny; ++i) {
            if (params.b_down == 0) {   //wall
                m[i][params.fict - j] = m[i][params.fict + j - 1];
                impx[i][params.fict - j] = impx[i][params.fict + j - 1];
                impy[i][params.fict - j] = -impy[i][params.fict + j - 1];
                e[i][params.fict - j] = e[i][params.fict + j - 1];
            }
            else if (params.b_down == 1) {   //free boundery
                m[i][params.fict - j] = m[i][params.fict + j - 1];
                impx[i][params.fict - j] = impx[i][params.fict + j - 1];
                impy[i][params.fict - j] = impy[i][params.fict + j - 1];
                e[i][params.fict - j] = e[i][params.fict + j - 1];
            }
            if (params.b_up == 0) {   //wall
                m[i][params.Nx - params.fict + j - 1] = m[i][params.Nx - params.fict - j];
                impx[i][params.Nx - params.fict + j - 1] = impx[i][params.Nx - params.fict - j];
                impy[i][params.Nx - params.fict + j - 1] = -impy[i][params.Nx - params.fict - j];
                e[i][params.Nx - params.fict + j - 1] = e[i][params.Nx - params.fict - j];
            }
            else if (params.b_up == 1) {   //free boundery
                m[i][params.Nx - params.fict + j - 1] = m[i][params.Nx - params.fict - j];
                impx[i][params.Nx - params.fict + j - 1] = impx[i][params.Nx - params.fict - j];
                impy[i][params.Nx - params.fict + j - 1] = impy[i][params.Nx - params.fict - j];
                e[i][params.Nx - params.fict + j - 1] = e[i][params.Nx - params.fict - j];
            }
        }
    }
}

//caluclate t_start step
double get_dt(InitialState& params, vec1d& x, vec1d& y, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    double new_step = 1.e6;
    double p, vx, vy, r, c;
    double c_step;
    double CFL = params.CFL;
    for (int i = params.fict; i < params.Ny - params.fict; ++i) {
        for (int j = params.fict; j < params.Nx - params.fict; ++j) {
            cons_to_noncons(params, p, vx, vy, r, m[i][j], impx[i][j], impy[i][j], e[i][j]);
            c = calc_c(params.gamma, p, r);
            c_step = std::min(CFL * (x[j + 1] - x[j]) / (std::fabs(vx) + c), CFL * (y[i + 1] - y[i]) / (std::fabs(vy) + c));
            if (c_step < new_step) {
                new_step = c_step;
            }
        }
    }
    return new_step;
}

/*
void exchange_dt(InitialState& params, double& dt, vec1d& x, vec1d& y, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e, int myrank, int size, MPI_Comm comm) {
    MPI_Status Status;
    if (myrank != 0) {
        MPI_Send(&dt, 1, MPI_DOUBLE, 0, myrank, comm);
        MPI_Recv(&dt, 1, MPI_DOUBLE, 0, myrank, comm, &Status);
    }
    if (myrank == 0) {
        double temp;
        for (int i = 1; i < size; ++i) {
            MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, comm, &Status);
            if (temp < dt) dt = temp;
        }
        for (int i = 1; i < size; ++i) {
            MPI_Send(&dt, 1, MPI_DOUBLE, i, i, comm);
        }
    }
}
*/
//initialize 0 t_start layer

void Soda(InitialState& params, vec1d& xc, vec1d& yc, vec2d& p, vec2d& vx, vec2d& vy, vec2d& r, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    double p1 = 1.0;
    double r1 = 1.0;
    double p2 = 0.1;
    double r2 = 0.125;

    for (int i = 0; i < params.Ny; ++i) {
        for (int j = 0; j < params.Nx; ++j) {
            if (xc[j] <= 0.5) {
                r[i][j] = r1;
                p[i][j] = p1;
                vx[i][j] = 0.;
                vy[i][j] = 0.;
            }
            else {
                r[i][j] = r2;
                p[i][j] = p2;
                vx[i][j] = 0.;
                vy[i][j] = 0.;
            }
            noncons_to_cons(params, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
        }
    }
}

void Gelmgolc(InitialState& params, vec1d& xc, vec1d& yc, vec2d& p, vec2d& vx, vec2d& vy, vec2d& r, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    double p1 = 1.0;
    double r1 = 1.0;
    double vx1 = 0.;
    double vy1 = 0.5;
    double p2 = 1.;
    double r2 = 2.;
    double vx2 = 0.;
    double vy2 = -0.5;

    for (int i = 0; i < params.Ny; ++i) {
        for (int j = 0; j < params.Nx; ++j) {
            if (xc[j] <= 0.5) {
                r[i][j] = r1;
                p[i][j] = p1;
                vx[i][j] = vx1;
                vy[i][j] = vy1;
            }
            else {
                r[i][j] = r2;
                p[i][j] = p2;
                vx[i][j] = vx2;
                vy[i][j] = vy2;
            }
            noncons_to_cons(params, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
        }
    }
}

void Bubble(InitialState& params, vec1d& xc, vec1d& yc, vec2d& p, vec2d& vx, vec2d& vy, vec2d& r, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    double R = 0.2;
    double p1, vx1, vy1, r1, p2, vx2, vy2, r2;
    p1 = 1.;
    vx1 = 0.5;
    vy1 = 0.;
    r1 = 1.;

    p2 = 1.;
    vx2 = 0.;
    vy2 = 0;
    r2 = 0.125;

    for (int i = 0; i < params.Ny; ++i) {
        for (int j = 0; j < params.Nx; ++j) {
            if (pow((xc[j] - 0.5), 2.) + pow((yc[i] - 0.5), 2.) < pow(R, 2.)) {
                p[i][j] = p1;
                vx[i][j] = vx1;
                vy[i][j] = vy1;
                r[i][j] = r1;
            }
            else {
                p[i][j] = p2;
                vx[i][j] = vx2;
                vy[i][j] = vy2;
                r[i][j] = r2;
            }
            noncons_to_cons(params, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
        }
    }
}

void Explosion(InitialState& params, vec1d& xc, vec1d& yc, vec2d& p, vec2d& vx, vec2d& vy, vec2d& r, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    double R = 0.2;
    double p1, vx1, vy1, r1, p2, vx2, vy2, r2;
    p1 = 1.;
    vx1 = 0.;
    vy1 = 0.;
    r1 = 1.;

    p2 = 0.1;
    vx2 = 0.;
    vy2 = 0.;
    r2 = 0.125;

    for (int i = 0; i < params.Ny; ++i) {
        for (int j = 0; j < params.Nx; ++j) {
            if (pow((xc[j] - 0.5), 2.) + pow((yc[i] - 0.5), 2.) < pow(R, 2.)) {
                p[i][j] = p1;
                vx[i][j] = vx1;
                vy[i][j] = vy1;
                r[i][j] = r1;
            }
            else {
                p[i][j] = p2;
                vx[i][j] = vx2;
                vy[i][j] = vy2;
                r[i][j] = r2;
            }
            noncons_to_cons(params, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
        }
    }
}

void Init(InitialState& params, vec1d& xc, vec1d& yc, vec2d& p, vec2d& vx, vec2d& vy, vec2d& r, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    switch (params.initial) {
    case 1:
        Bubble(params, xc, yc, p, vx, vy, r, m, impx, impy, e);
        break;
    case 2:
        Soda(params, xc, yc, p, vx, vy, r, m, impx, impy, e);
        break;
    case 3:
        Explosion(params, xc, yc, p, vx, vy, r, m, impx, impy, e);
        break;
    case 4:
        Gelmgolc(params, xc, yc, p, vx, vy, r, m, impx, impy, e);
        break;
    default:
        Bubble(params, xc, yc, p, vx, vy, r, m, impx, impy, e);
        break;
    }
}

//save results in the directory out_dir
void writeCSV(std::string filename, InitialState& IS, int iter, vec1d& xc, vec1d& yc, vec2d& ux, vec2d& uy, vec2d& P, vec2d& rho, double _time) {
    std::string name = filename;
    name += "/Iter=";
    name += std::to_string(iter);
    name += ".csv";
    if (!fs::exists(filename)) {
        fs::create_directory(filename);
    }
    std::ofstream csv;
    csv.open(name);
    csv << "Time: " << _time << std::endl;
    csv << "X,Y,Rho,P,Ux,Uy\n";
    for (int i = IS.fict; i < IS.Ny - IS.fict; ++i) {
        for (int j = IS.fict; j < IS.Nx - IS.fict; ++j) {
            csv << xc[j] << "," << yc[i] << "," << rho[i][j] << "," << P[i][j] << "," << ux[i][j] << "," << uy[i][j] << "\n";
        }
    }
    csv.close();
}

/*
void writeCSV_p(std::string filename, InitialState& IS, int iter, vec1d& xc, vec1d& yc, vec2d& ux, vec2d& uy, vec2d& P, vec2d& rho,
    double _time, int myrank, MPI_Comm comm) {
    if (myrank == 0) {
        if (!fs::exists(filename)) {
            fs::create_directory(filename);
        }
    }
    MPI_Barrier(comm);

    std::string name = "CSVs\\";
    name += filename;
    name += "\\Iter=";
    name += std::to_string(iter);
    name += ".csv";

    std::ostringstream buffer;

    if (myrank == 0) {
        buffer << (std::to_string(_time) + "\n") << "X,Y,Rho,P,Ux,Uy\n";
    }
    for (int i = IS.fict; i < IS.Ny - IS.fict; ++i) {
        for (int j = IS.fict; j < IS.Nx - IS.fict; ++j) {
            buffer << xc[i] << "," << yc[i] << "," << rho[i][j] << ',' << P[i][j] << ',' << ux[i][j] << ',' << uy[i][j] << "\n";
        }
    }

    std::string local_data = buffer.str();
    int local_size = local_data.size();

    int offset = 0;
    MPI_Exscan(&local_size, &offset, 1, MPI_INT, MPI_SUM, comm);
    MPI_File fh;
    MPI_File_open(comm, name.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write_at_all(fh, offset, local_data.c_str(), local_size, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
}
*/

double Newton_guess_P(InitialState& params, double pl, double vl, double rl, double cl, double pr, double vr, double rr, double cr) {
    double p_lin; //linear approximation
    double p_min, p_max; //min max pressure
    double p_ratio;
    double p_ratio_max = 2.0;
    double tl, tr;
    double g1 = (params.gamma - 1.0) / 2. / params.gamma;
    double g2 = (params.gamma - 1.0) / (params.gamma + 1.0);
    double g3 = 2.0 / (params.gamma + 1.0);

    double eps = 1.e-6;

    p_lin = std::max(eps, 0.5 * (pl + pr) + 0.5 * (vl - vr) * 0.25 * (rl + rr) * (cl + cr)); //Toro P.299 (9.20, 9.21)
    p_min = std::min(pl, pr);
    p_max = std::max(pl, pr);
    p_ratio = p_max / p_min;

    if ((p_ratio <= p_ratio_max) && ((p_min < p_lin && p_lin < p_max))) {
        return p_lin;
    }
    else {
        if (p_lin < p_min) {
            //two-rarefaction
            return pow(((cl + cr - 0.5 * (params.gamma - 1.0) * (vr - vl)) / (cl / pow(pl, g1) + cr / pow(pr, g1))), 1. / g1); //Toro P.301 (9.32)
        }
        else {
            // two-shock approximation
            tl = std::sqrt(g3 / rl / (g2 * pl + p_lin));
            tr = std::sqrt(g3 / rr / (g2 * pr + p_lin));
            return (tl * pl + tr * pr - (vr - vl)) / (tl + tr); //Toro P.303 (9.42)
        }
    }
}


//to evaluate the pressure functions FL and FR in Riemann solver
void calc_f(InitialState& IS, double& F, double& FD, double p_curr, double p, double v, double r, double c) {
    double g1 = (IS.gamma - 1.0) / 2.0 / IS.gamma;
    double g2 = 0.5 * (IS.gamma + 1.0) / IS.gamma;
    double g4 = 2.0 / (IS.gamma - 1.0);
    double g5 = 2.0 / (IS.gamma + 1.0);
    double g6 = (IS.gamma - 1.0) / (IS.gamma + 1.0);
    double p_ratio = p_curr / p;
    double qrt;

    if (p_curr <= p) {
        //rarefaction
        F = g4 * c * (pow(p_ratio, g1) - 1.0); //Toro P.301 (9.30)
        FD = (1.0 / r / c) * pow(p_ratio, -g2); // dF/dp
    }
    else {
        //shock wave
        qrt = sqrt((g5 / r) / (g6 * p + p_curr));
        F = (p_curr - p) * qrt; //Toro P.301 (9.30)
        FD = (1.0 - 0.5 * (p_curr - p) / (g6 * p + p_curr)) * qrt; // dF/dp
    }
}

//to sample the solution throughout the wave pattern. Sampling is performed in terms of the ’speed’ S = X/T.
void find_in_which_part(InitialState& params, double pl, double vl, double rl, double cl, double pr, double vr, double rr, double cr, double p_cont, double v_cont, double s, double& p_res, double& v_res, double& r_res) {

    double g1 = (params.gamma - 1.0) / 2. / params.gamma;
    double g2 = (params.gamma + 1.0) / 2. / params.gamma;
    double g3 = 2.0 * params.gamma / (params.gamma - 1.0);
    double g4 = 2.0 / (params.gamma - 1.0);
    double g5 = 2.0 / (params.gamma + 1.0);
    double g6 = (params.gamma - 1.0) / (params.gamma + 1.0);
    double g7 = (params.gamma - 1.0) / 2.;
    double g8 = params.gamma - 1.0;

    //left shock wave
    double shl, stl;
    double sl;
    //right shock wave
    double shr, str;
    double sr;

    double cml, cmr; //sound velocity left and right of contact discontinuity
    double crf; //sound velocity in rarefaction
    double p_ratio;
    double r, v, p; //sample values

    if (s <= v_cont) {
        //Sampling point lies to the left of the contact discontinuity
        if (p_cont <= pl) {
            //Left rarefaction
            shl = vl - cl;
            if (s <= shl) {
                //слева от КР
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                cml = cl * pow(p_cont / pl, g1);
                stl = v_cont - cml;
                if (s > stl) {
                    //слева от КР
                    r = rl * pow(p_cont / pl, 1.0 / params.gamma); //Toro P.134 (4.53)
                    v = v_cont;
                    p = p_cont;
                }
                else {
                    //inside left rarefaction
                    //Toro P.135 (4.56)
                    v = g5 * (cl + g7 * vl + s);
                    crf = g5 * (cl + g7 * (vl - s));
                    r = rl * pow(crf / cl, g4);
                    p = pl * pow(crf / cl, g3);
                }
            }
        }
        else {
            //left shock
            p_ratio = p_cont / pl;
            sl = vl - cl * std::sqrt(g2 * p_ratio + g1);
            if (s <= sl) {
                //слева от КР
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                //за левой УВ
                r = rl * (p_ratio + g6) / (p_ratio * g6 + 1.0); //Toro P.133 (4.50)
                v = v_cont;
                p = p_cont;
            }
        }
    }
    else {
        if (p_cont > pr) { //Sampling point lies to the left of the contact discontinuity
            //right shock
            p_ratio = p_cont / pr;
            sr = vr + cr * std::sqrt(g2 * p_ratio + g1);
            if (s >= sr) {
                //справа от КР
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                //за правой УВ
                r = rr * (p_ratio + g6) / (p_ratio * g6 + 1.0); //Toro P.135 (4.57)
                v = v_cont;
                p = p_cont;
            }
        }
        else {
            //right rarefaction
            shr = vr + cr;
            if (s >= shr) {
                //справа от КР
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                cmr = cr * pow(p_cont / pr, g1);
                str = v_cont + cmr;
                if (s <= str) {
                    //справа от КР
                    r = rr * pow(p_cont / pr, 1.0 / params.gamma); //Toro P.136 (4.60)
                    v = v_cont;
                    p = p_cont;
                }
                else {
                    //inside right rarefaction
                    //Toro P.136 (4.63)
                    v = g5 * (-cr + g7 * vr + s);
                    crf = g5 * (cr - g7 * (vr - s));
                    r = rr * pow(crf / cr, g4);
                    p = pr * pow(crf / cr, g3);
                }
            }
        }
    }
    //result values
    r_res = r;
    v_res = v;
    p_res = p;

}

//Toro P.155. Subroutine STARPU
void Newton_contact(InitialState& params, double pl, double vl, double rl, double cl, double pr, double vr, double rr, double cr, double& p_cont, double& v_cont) {
    double p_old; //previous iteration
    double fl, fr;
    double fl_diff, fr_diff; //derivative of f

    int step = 0;
    int max_iter = 500;

    double change; //convergence parameter
    double eps = 1.e-6; //accuracy

    if (2. * (cl + cr) / (params.gamma - 1.) <= (vr - vl)) { //Toro P.127 (4.40)
        std::cerr << "Vacuum was generated" << std::endl;
    }

    //initial pressure for iterational algorithm
    p_old = Newton_guess_P(params, pl, vl, rl, cl, pr, vr, rr, cr);
    if (p_old < 0.) {
        p_old = eps;
        std::cerr << "Initial pressure < 0" << std::endl;
    }

    //itterational algorithm
    do {
        calc_f(params, fl, fl_diff, p_old, pl, vl, rl, cl);
        calc_f(params, fr, fr_diff, p_old, pr, vr, rr, cr);
        p_cont = p_old - (fl + fr + vr - vl) / (fl_diff + fr_diff);
        change = 2.0 * std::fabs((p_cont - p_old) / (p_cont + p_old));
        step++;
        if (step > max_iter) {
            std::cerr << "Newton algorhitm doesn't converge with max iterations" << std::endl;
        }
        if (p_cont < 0.0) {
            p_cont = eps;
            std::cerr << "Contact pressure < 0" << std::endl;
        }
        p_old = p_cont;
    } while (change > eps);
    v_cont = 0.5 * (vl + vr + fr - fl); //Toro P.125
}


void Riemann_solver(InitialState& params, double rl, double vl, double pl, double rr, double vr, double pr, double& p_res, double& v_res, double& r_res) {
    double p_cont, v_cont;
    double p, v, r;
    double cr, cl;

    cl = calc_c(params.gamma, pl, rl);
    cr = calc_c(params.gamma, pr, rr);

    if (2. * (cl + cr) / (params.gamma - 1.) <= (vr - vl)) {
        std::cerr << "Vacuum was generated" << std::endl;
    }

    Newton_contact(params, pl, vl, rl, cl, pr, vr, rr, cr, p_cont, v_cont); // contact pressure and velocity

    find_in_which_part(params, pl, vl, rl, cl, pr, vr, rr, cr, p_cont, v_cont, 0.0, p, v, r); //solution selection

    p_res = p;
    v_res = v;
    r_res = r;
}


void Godunov_flux_x(InitialState& params, double p, double vx, double vy, double r, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
    double m, impx, impy, e;
    noncons_to_cons(params, p, vx, vy, r, m, impx, impy, e);

    Fm = r * vx;
    Fimpx = Fm * vx + p;
    Fimpy = Fm * vy;
    Fe = (p + e) * vx;
}


void Godunov_flux_y(InitialState& params, double p, double vx, double vy, double r, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
    double m, impx, impy, e;
    noncons_to_cons(params, p, vx, vy, r, m, impx, impy, e);

    Fm = r * vy;
    Fimpx = Fm * vx;
    Fimpy = Fm * vy + p;
    Fe = (p + e) * vy;
}


void cons_flux_x(InitialState& params, double m, double impx, double impy, double e, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
    double p, vx, vy, r;
    cons_to_noncons(params, p, vx, vy, r, m, impx, impy, e);

    Fm = r * vx;
    Fimpx = Fm * vx + p;
    Fimpy = Fm * vy;
    Fe = (p + e) * vx;
}


void cons_flux_y(InitialState& params, double m, double impx, double impy, double e, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
    double p, vx, vy, r;
    cons_to_noncons(params, p, vx, vy, r, m, impx, impy, e);

    Fm = r * vy;
    Fimpx = Fm * vx;
    Fimpy = Fm * vy + p;
    Fe = (p + e) * vy;
}


void Godunov_solver_x(InitialState& params, double ml, double impxl, double impyl, double el, double mr, double impxr, double impyr, double er, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
    double p, vx, vy, r;
    double pl, vxl, vyl, rl;
    double pr, vxr, vyr, rr;

    cons_to_noncons(params, pl, vxl, vyl, rl, ml, impxl, impyl, el);
    cons_to_noncons(params, pr, vxr, vyr, rr, mr, impxr, impyr, er);

    Riemann_solver(params, rl, vxl, pl, rr, vxr, pr, p, vx, r);

    if (vx >= 0)
        vy = impyl / ml;
    else
        vy = impyr / mr;

    Godunov_flux_x(params, p, vx, vy, r, Fm, Fimpx, Fimpy, Fe);
}


void Godunov_solver_y(InitialState& params, double md, double impxd, double impyd, double ed, double mu, double impxu, double impyu, double eu, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
    double p, vx, vy, r;
    double pd, vxd, vyd, rd;
    double pu, vxu, vyu, ru;

    cons_to_noncons(params, pd, vxd, vyd, rd, md, impxd, impyd, ed);
    cons_to_noncons(params, pu, vxu, vyu, ru, mu, impxu, impyu, eu);

    Riemann_solver(params, rd, vyd, pd, ru, vyu, pu, p, vy, r);

    if (vy >= 0)
        vx = impxd / md;
    else
        vx = impxu / mu;

    Godunov_flux_y(params, p, vx, vy, r, Fm, Fimpx, Fimpy, Fe);
}

double minmod(double a, double b, int func) {
    if (func == 0)
        return (abs(a) < abs(b) ? a : b);
    else if (func == 1)
    {
        double c = (a + b) / 2;
        return ((abs(a) < abs(b) && abs(a) < abs(c)) ? a : (abs(b) < abs(a) && abs(b) < abs(c) ? b : c));
    }
    else
        return (a * b < 0 ? 0 : (a * a < a * b ? a : b));
}


//Godunov-Kolgan reconstruction on the left cell boundary
void Kolgan_left(double qi, double qi_1, double qi_2, double qi1, double xi, double xi_1, double xi_2, double xi1, double& ql, double& qr, int mod_type) {
    double dif_l, dif_r, dif_c; //linear coeffs

    dif_l = (qi - qi_1) / (xi - xi_1);
    dif_r = (qi1 - qi) / (xi1 - xi);

    dif_c = minmod(dif_l, dif_r, mod_type);

    qr = -dif_c * (xi - xi_1) / 2. + qi;

    dif_r = dif_l;
    dif_l = (qi_1 - qi_2) / (xi_1 - xi_2);

    dif_c = minmod(dif_l, dif_r, mod_type);

    ql = dif_c * (xi - xi_1) / 2. + qi_1;
}

//Godunov-Kolgan reconstruction on the right cell boundary
void Kolgan_right(double qi, double qi_1, double qi1, double qi2, double xi, double xi_1, double xi1, double xi2, double& ql, double& qr, int mod_type) {
    double dif_l, dif_r, dif_c; //linear coeffs

    dif_l = (qi - qi_1) / (xi - xi_1);
    dif_r = (qi1 - qi) / (xi1 - xi);

    dif_c = minmod(dif_l, dif_r, mod_type);

    ql = dif_c * (xi1 - xi) / 2. + qi;

    dif_l = dif_r;
    dif_r = (qi2 - qi1) / (xi2 - xi1);

    dif_c = minmod(dif_l, dif_r, mod_type);

    qr = -dif_c * (xi1 - xi) / 2. + qi1;
}


void calc_Fx(InitialState& params, vec3d& u, vec3d& F) {
    double r, vx, vy, p, e;
    for (int i = 0; i < params.Ny; i++) {
        for (int j = 0; j < params.Nx; j++) {
            r = u[0][i][j];
            vx = u[1][i][j] / r;
            vy = u[2][i][j] / r;
            e = u[3][i][j];
            p = (params.gamma - 1.0) * (e - 0.5 * r * (vx * vx + vy * vy));
            F[0][i][j] = vx * r;
            F[1][i][j] = vx * vx * r + p;
            F[2][i][j] = vx * vy * r;
            F[3][i][j] = vx * (e + p);
        }
    }
}

void calc_Fy(InitialState& params, vec3d& u, vec3d& F) {
    double r, vx, vy, p, e;
    for (int i = 0; i < params.Ny; i++) {
        for (int j = 0; j < params.Nx; j++) {
            r = u[0][i][j];
            vx = u[1][i][j] / r;
            vy = u[2][i][j] / r;
            e = u[3][i][j];
            p = (params.gamma - 1.0) * (e - 0.5 * r * (vx * vx + vy * vy));
            F[0][i][j] = vy * r;
            F[1][i][j] = vy * vx * r;
            F[2][i][j] = vy * vy * r + p;
            F[3][i][j] = vy * (e + p);
        }
    }
}

double calc_Ax(InitialState& params, vec3d& u) {
    double new_step = 0.;
    double p, vx, vy, r, e, c;
    double a_step;
    for (int i = params.fict; i < params.Ny - params.fict; ++i) {
        for (int j = params.fict; j < params.Nx - params.fict; ++j) {
            r = u[0][i][j];
            vx = u[1][i][j] / r;
            vy = u[2][i][j] / r;
            e = u[3][i][j];
            p = (params.gamma - 1.0) * (e - 0.5 * r * (vx * vx + vy * vy));
            c = calc_c(params.gamma, p, r);
            a_step = abs(vx) + c;
            if (a_step > new_step) {
                new_step = a_step;
            }
        }
    }
    return new_step;
}

double calc_Ay(InitialState& params, vec3d& u) {
    double new_step = 0.;
    double p, vx, vy, r, e, c;
    double a_step;
    for (int i = params.fict; i < params.Ny - params.fict; ++i) {
        for (int j = params.fict; j < params.Nx - params.fict; ++j) {
            r = u[0][i][j];
            vx = u[1][i][j] / r;
            vy = u[2][i][j] / r;
            e = u[3][i][j];
            p = (params.gamma - 1.0) * (e - 0.5 * r * (vx * vx + vy * vy));
            c = calc_c(params.gamma, p, r);
            a_step = abs(vy) + c;
            if (a_step > new_step) {
                new_step = a_step;
            }
        }
    }
    return new_step;
}

double WENO_reconstruction_minus(InitialState& params, double um2, double um1, double u, double up1, double up2) {
    double eps = 1.e-40;
    double a_sum = 0.;
    vec1d gammaRight = { 0.1, 0.6, 0.3 };
    vec1d a(3);
    vec1d b(3);
    vec1d P(3);
    vec1d w(3);

    b[2] = 13.0 / 12.0 * pow((um2 - 2 * um1 + u), 2) + 1.0 / 4.0 * pow((um2 - 4 * um1 + 3 * u), 2);
    b[1] = 13.0 / 12.0 * pow((um1 - 2 * u + up1), 2) + 1.0 / 4.0 * pow((um1 - up1), 2);
    b[0] = 13.0 / 12.0 * pow((u - 2 * up1 + up2), 2) + 1.0 / 4.0 * pow((3 * u - 4 * up1 + up2), 2);

    P[0] = 11.0 / 6.0 * u - 7.0 / 6.0 * up1 + 1.0 / 3.0 * up2;
    P[1] = 1.0 / 3.0 * um1 + 5.0 / 6.0 * u - 1.0 / 6.0 * up1;
    P[2] = -1.0 / 6.0 * um2 + 5.0 / 6.0 * um1 + 1.0 / 3.0 * u;

    for (int j = 0; j < 3; j++) {
        a[j] = gammaRight[j] / pow((b[j] + eps), 2.0);
        a_sum += a[j];
    }

    for (int j = 0; j < 3; j++) {
        w[j] = a[j] / a_sum;
    }
    return w[0] * P[0] + w[1] * P[1] + w[2] * P[2];
}

double WENO_reconstruction_plus(InitialState& params, double um2, double um1, double u, double up1, double up2) {
    double eps = 1.e-40;
    double a_sum = 0.;
    vec1d gammaLeft = { 0.3, 0.6, 0.1 };
    vec1d a(3);
    vec1d b(3);
    vec1d P(3);
    vec1d w(3);

    b[2] = 13.0 / 12.0 * pow((um2 - 2 * um1 + u), 2) + 1.0 / 4.0 * pow((um2 - 4 * um1 + 3 * u), 2);
    b[1] = 13.0 / 12.0 * pow((um1 - 2 * u + up1), 2) + 1.0 / 4.0 * pow((um1 - up1), 2);
    b[0] = 13.0 / 12.0 * pow((u - 2 * up1 + up2), 2) + 1.0 / 4.0 * pow((3 * u - 4 * up1 + up2), 2);

    P[0] = 1.0 / 3.0 * u + 5.0 / 6.0 * up1 - 1.0 / 6.0 * up2;
    P[1] = -1.0 / 6.0 * um1 + 5.0 / 6.0 * u + 1.0 / 3.0 * up1;
    P[2] = 1.0 / 3.0 * um2 - 7.0 / 6.0 * um1 + 11.0 / 6.0 * u;

    for (int j = 0; j < 3; j++) {
        a[j] = gammaLeft[j] / pow((b[j] + eps), 2.0);
        a_sum += a[j];
    }

    for (int j = 0; j < 3; j++) {
        w[j] = a[j] / a_sum;
    }
    return w[0] * P[0] + w[1] * P[1] + w[2] * P[2];
}

void WENO5(InitialState& IS, vec3d& u, vec3d& L, double dx, double dy) {
    vec3d Fx(4, vec2d(IS.Ny, vec1d(IS.Nx)));
    vec3d Fy(4, vec2d(IS.Ny, vec1d(IS.Nx)));
    vec3d Fx_plus(4, vec2d(IS.Ny, vec1d(IS.Nx)));
    vec3d Fy_plus(4, vec2d(IS.Ny, vec1d(IS.Nx)));
    vec3d Fx_minus(4, vec2d(IS.Ny, vec1d(IS.Nx)));
    vec3d Fy_minus(4, vec2d(IS.Ny, vec1d(IS.Nx)));


    calc_Fx(IS, u, Fx);
    calc_Fy(IS, u, Fy);
    double Ax = calc_Ax(IS, u);
    double Ay = calc_Ay(IS, u);

    for (int k = 0; k < 4; k++) {
        for (int i = 0; i < IS.Ny; i++) {
            for (int j = 0; j < IS.Nx; j++) {
                Fx_plus[k][i][j] = 0.5 * (Fx[k][i][j] + Ax * u[k][i][j]);
                Fx_minus[k][i][j] = 0.5 * (Fx[k][i][j] - Ax * u[k][i][j]);
                Fy_plus[k][i][j] = 0.5 * (Fy[k][i][j] + Ay * u[k][i][j]);
                Fy_minus[k][i][j] = 0.5 * (Fy[k][i][j] - Ay * u[k][i][j]);
            }
        }
    }

    vec3d FlowL(4, vec2d(IS.Ny, vec1d(IS.Nx)));
    vec3d FlowR(4, vec2d(IS.Ny, vec1d(IS.Nx)));
    for (int k = 0; k < 4; k++) {
        for (int i = IS.fict; i < IS.Ny - IS.fict; i++) {
            for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
                FlowL[k][i][j] = WENO_reconstruction_minus(IS, Fx_minus[k][i - 2][j], Fx_minus[k][i - 1][j], Fx_minus[k][i][j], Fx_minus[k][i + 1][j], Fx_minus[k][i + 2][j]);
                FlowR[k][i][j] = WENO_reconstruction_plus(IS, Fx_plus[k][i - 2][j], Fx_plus[k][i - 1][j], Fx_plus[k][i][j], Fx_plus[k][i + 1][j], Fx_plus[k][i + 2][j]);
            }
        }
    }
    Boundary_x(IS, FlowL[0], FlowL[1], FlowL[2], FlowL[3]);
    Boundary_x(IS, FlowR[0], FlowR[1], FlowR[2], FlowR[3]);

    vec3d FlowD(4, vec2d(IS.Ny, vec1d(IS.Nx)));
    vec3d FlowU(4, vec2d(IS.Ny, vec1d(IS.Nx)));
    for (int k = 0; k < 4; k++) {
        for (int i = IS.fict; i < IS.Ny - IS.fict; i++) {
            for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
                FlowD[k][i][j] = WENO_reconstruction_minus(IS, Fy_minus[k][i][j - 2], Fy_minus[k][i][j - 1], Fy_minus[k][i][j], Fy_minus[k][i][j + 1], Fy_minus[k][i][j + 2]);
                FlowU[k][i][j] = WENO_reconstruction_plus(IS, Fy_plus[k][i][j - 2], Fy_plus[k][i][j - 1], Fy_plus[k][i][j], Fy_plus[k][i][j + 1], Fy_plus[k][i][j + 2]);
            }
        }
    }
    Boundary_y(IS, FlowD[0], FlowD[1], FlowD[2], FlowD[3]);
    Boundary_y(IS, FlowU[0], FlowU[1], FlowU[2], FlowU[3]);

    for (int k = 0; k < 4; k++) {
        for (int i = IS.fict; i < IS.Ny - IS.fict; i++) {
            for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
                L[k][i][j] = ((FlowR[k][i][j] + FlowL[k][i + 1][j]) - (FlowR[k][i - 1][j] + FlowL[k][i][j])) / dx + \
                    ((FlowU[k][i][j] + FlowD[k][i][j + 1]) - (FlowU[k][i][j - 1] + FlowD[k][i][j])) / dy;
            }
        }
    }
}

/*
void mpi_get_px_py(int p, int Nx, int Ny, int* px, int* py) {
    if (p == 1) {
        (*px) = 1;
        (*py) = 1;
    }
    else {
        int best_px, best_py;
        int best_p;
        int mera = -1;
        for (int i = 2; i < p + 1; i++) {
            int first, second;
            for (int j = 1; j < sqrt(i) + 1; j++) {
                if (i % j == 0) {
                    first = j;
                    second = i / j;
                }
            }
            if (i - abs(second - first) > mera) {
                mera = i - abs(second - first);
                best_px = first;
                best_py = second;
                best_p = i;
            }
        }
        if ((Nx % best_px) * (Ny % best_py) < (Nx % best_py) * (Ny % best_px)) {
            (*px) = best_px;
            (*py) = best_py;
        }
        else {
            (*px) = best_py;
            (*py) = best_px;
        }
    }
}


void sending(InitialState IS, vec2d& P, vec2d& ux, vec2d& uy, vec2d& rho, int* neighbors, int rank, MPI_Status* Status) {
    int length;
    int recv_len;
    vec1d send_arr;
    for (int d = 0; d < 4; d++) {
        if (neighbors[d] > 0) {
            // packing data to send
            if (d == 0) {
                length = 4 * (IS.Ny - 2 * IS.fict) * IS.fict;
                send_arr.reserve(length);
                for (int i = IS.fict; i < IS.Ny - IS.fict; i++) {
                    for (int j = IS.fict; j < IS.fict + IS.fict; j++) {
                        for (int k = 0; k < 4; k++) {
                            if (k == 0) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = P[i][j];
                            if (k == 1) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = ux[i][j];
                            if (k == 2) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = uy[i][j];
                            if (k == 3) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = rho[i][j];
                        }
                    }
                }
            }
            else if (d == 1) {
                length = 4 * (IS.Nx - 2 * IS.fict) * IS.fict;
                send_arr.reserve(length);
                for (int i = IS.Ny - 2 * IS.fict; i < IS.Ny - IS.fict; i++) {
                    for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
                        for (int k = 0; k < 4; k++) {
                            if (k == 0) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = P[i][j];
                            if (k == 1) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = ux[i][j];
                            if (k == 2) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = uy[i][j];
                            if (k == 3) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = rho[i][j];
                        }
                    }
                }
            }
            else if (d == 2) {
                length = 4 * (IS.Ny - 2 * IS.fict) * IS.fict;
                send_arr.reserve(length);
                for (int i = IS.fict; i < IS.Ny - IS.fict; i++) {
                    for (int j = IS.Nx - 2 * IS.fict; j < IS.Nx - IS.fict; j++) {
                        for (int k = 0; k < 4; k++) {
                            if (k == 0) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = P[i][j];
                            if (k == 1) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = ux[i][j];
                            if (k == 2) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = uy[i][j];
                            if (k == 3) send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k] = rho[i][j];
                        }
                    }
                }
            }
            else if (d == 3) {
                length = 4 * (IS.Nx - 2 * IS.fict) * IS.fict;
                send_arr.reserve(length);
                for (int i = IS.fict; i < IS.fict + IS.fict; i++) {
                    for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
                        for (int k = 0; k < 4; k++) {
                            if (k == 0) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = P[i][j];
                            if (k == 1) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = ux[i][j];
                            if (k == 2) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = uy[i][j];
                            if (k == 3) send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k] = rho[i][j];
                        }
                    }
                }
            }
            MPI_Send(&(send_arr[0]), length, MPI_DOUBLE, neighbors[d], neighbors[d], MPI_COMM_WORLD);
            MPI_Recv(&(recv_len), 1, MPI_INT, neighbors[d], rank, MPI_COMM_WORLD, Status);
            std::vector<double> recv_arr;
            MPI_Recv(&(recv_arr[0]), recv_len, MPI_DOUBLE, neighbors[d], rank, MPI_COMM_WORLD, Status);
            // unpacking received data
            if (d == 0) {
                if (recv_len != 4 * (IS.Ny - 2 * IS.fict) * IS.fict) {
                    std::cerr << "Error in sending" << std::endl;
                }
                for (int i = IS.fict; i < IS.Ny - IS.fict; i++) {
                    for (int j = IS.Nx - IS.fict; j < IS.Nx; j++) {
                        for (int k = 0; k < 4; k++) {
                            if (k == 0) P[i][j] = recv_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
                            if (k == 1) ux[i][j] = recv_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
                            if (k == 2) uy[i][j] = recv_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
                            if (k == 3) rho[i][j] = recv_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
                        }
                    }
                }
            }
            else if (d == 1) {
                if (recv_len != 4 * (IS.Nx - 2 * IS.fict) * IS.fict) {
                    std::cerr << "Error in sending" << std::endl;
                }
                for (int i = 0; i < IS.fict; i++) {
                    for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
                        for (int k = 0; k < 4; k++) {
                            if (k == 0) P[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
                            if (k == 1) ux[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
                            if (k == 2) uy[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
                            if (k == 3) rho[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
                        }
                    }
                }
            }
            else if (d == 2) {
                if (recv_len != 4 * (IS.Ny - 2 * IS.fict) * IS.fict) {
                    std::cerr << "Error in sending" << std::endl;
                }
                for (int i = IS.fict; i < IS.Ny - IS.fict; i++) {
                    for (int j = 0; j < IS.fict; j++) {
                        for (int k = 0; k < 4; k++) {
                            if (k == 0) P[i][j] = send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
                            if (k == 1) ux[i][j] = send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
                            if (k == 2) uy[i][j] = send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
                            if (k == 3) rho[i][j] = send_arr[i * (IS.Ny - 2 * IS.fict) + j * IS.fict + k];
                        }
                    }
                }
            }
            else if (d == 3) {
                if (recv_len != 4 * (IS.Nx - 2 * IS.fict) * IS.fict) {
                    std::cerr << "Error in sending" << std::endl;
                }
                for (int i = IS.Ny - IS.fict; i < IS.Ny; i++) {
                    for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
                        for (int k = 0; k < 4; k++) {
                            if (k == 0) P[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
                            if (k == 1) ux[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
                            if (k == 2) uy[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
                            if (k == 3) rho[i][j] = send_arr[i * IS.fict + j * (IS.Nx - 2 * IS.fict) + k];
                        }
                    }
                }
            }
        }
    }
}

void send_vec2(vec2d& data, int rows, int cols, int receiver, int tag, MPI_Comm comm) {
    vec1d flat_data(rows * cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            flat_data[i * cols + j] = data[i][j];
        }
    }

    MPI_Send(flat_data.data(), rows * cols, MPI_DOUBLE, receiver, tag, comm);
}

vec2d recv_vec2(int rows, int cols, int sender, int tag, MPI_Comm comm, MPI_Status& Status) {
    vec1d flat_data(rows * cols);
    MPI_Recv(flat_data.data(), rows * cols, MPI_DOUBLE, sender, tag, comm, &Status);

    vec2d result(rows, vec1d(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i][j] = flat_data[i * cols + j];
        }
    }
    return result;
}

void exchange_data(InitialState& IS, int myrank, int px, int py, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e, vec2d& buff_x, vec2d& buff_y, MPI_Comm comm, MPI_Status& Status) {
    // X sending
    if ((myrank % px) % 2 != 0) {
        if ((myrank % px) > 0) {
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Ny + 2 * IS.fict; j++) {
                    buff_x[i][j] = m[i + IS.fict][j];
                    buff_x[i + IS.fict][j] = impx[i + IS.fict][j];
                    buff_x[i + 2 * IS.fict][j] = impy[i + IS.fict][j];
                    buff_x[i + 3 * IS.fict][j] = e[i + IS.fict][j];
                }
            }
            send_vec2(buff_x, 4 * IS.fict, (IS.Ny + 2 * IS.fict), myrank - 1, 0, comm);
        }
        if ((myrank % px) > 0) {
            buff_x = recv_vec2(4 * IS.fict, IS.Ny + 2 * IS.fict, myrank - 1, 0, comm, Status);
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Ny + 2 * IS.fict; j++) {
                    m[i][j] = buff_x[i][j];
                    impx[i][j] = buff_x[i + IS.fict][j];
                    impy[i][j] = buff_x[i + 2 * IS.fict][j];
                    e[i][j] = buff_x[i + 3 * IS.fict][j];
                }
            }
        }
        if ((myrank % px) < px - 1) {
            buff_x = recv_vec2(4 * IS.fict, IS.Ny + 2 * IS.fict, myrank + 1, 0, comm, Status);
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Ny + 2 * IS.fict; j++) {
                    m[IS.Nx + IS.fict + i][j] = buff_x[i][j];
                    impx[IS.Nx + IS.fict + i][j] = buff_x[i + IS.fict][j];
                    impy[IS.Nx + IS.fict + i][j] = buff_x[i + 2 * IS.fict][j];
                    e[IS.Nx + IS.fict + i][j] = buff_x[i + 3 * IS.fict][j];
                }
            }
        }
        if ((myrank % px) < px - 1) {
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Ny + 2 * IS.fict; j++) {
                    buff_x[i][j] = m[IS.Nx + i][j];
                    buff_x[i + IS.fict][j] = impx[IS.Nx + i][j];
                    buff_x[i + 2 * IS.fict][j] = impy[IS.Nx + i][j];
                    buff_x[i + 3 * IS.fict][j] = e[IS.Nx + i][j];
                }
            }
            send_vec2(buff_x, 4 * IS.fict, (IS.Ny + 2 * IS.fict), myrank + 1, 0, comm);
        }
    }
    if ((myrank % px) % 2 == 0) {
        if ((myrank % px) < px - 1) {
            buff_x = recv_vec2(4 * IS.fict, IS.Ny + 2 * IS.fict, myrank + 1, 0, comm, Status);
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Ny + 2 * IS.fict; j++) {
                    m[IS.Nx + IS.fict + i][j] = buff_x[i][j];
                    impx[IS.Nx + IS.fict + i][j] = buff_x[i + IS.fict][j];
                    impy[IS.Nx + IS.fict + i][j] = buff_x[i + 2 * IS.fict][j];
                    e[IS.Nx + IS.fict + i][j] = buff_x[i + 3 * IS.fict][j];
                }
            }
        }
        if ((myrank % px) < px - 1) {
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Ny + 2 * IS.fict; j++) {
                    buff_x[i][j] = m[IS.Nx + i][j];
                    buff_x[i + IS.fict][j] = impx[IS.Nx + i][j];
                    buff_x[i + 2 * IS.fict][j] = impy[IS.Nx + i][j];
                    buff_x[i + 3 * IS.fict][j] = e[IS.Nx + i][j];
                }
            }
            send_vec2(buff_x, 4 * IS.fict, (IS.Ny + 2 * IS.fict), myrank + 1, 0, comm);
        }
        if ((myrank % px) > 0) {
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Ny + 2 * IS.fict; j++) {
                    buff_x[i][j] = m[i + IS.fict][j];
                    buff_x[i + IS.fict][j] = impx[i + IS.fict][j];
                    buff_x[i + 2 * IS.fict][j] = impy[i + IS.fict][j];
                    buff_x[i + 3 * IS.fict][j] = e[i + IS.fict][j];
                }
            }
            send_vec2(buff_x, 4 * IS.fict, (IS.Ny + 2 * IS.fict), myrank - 1, 0, comm);
        }
        if ((myrank % px) > 0) {
            buff_x = recv_vec2(4 * IS.fict, IS.Ny + 2 * IS.fict, myrank - 1, 0, comm, Status);
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Ny + 2 * IS.fict; j++) {
                    m[i][j] = buff_x[i][j];
                    impx[i][j] = buff_x[i + IS.fict][j];
                    impy[i][j] = buff_x[i + 2 * IS.fict][j];
                    e[i][j] = buff_x[i + 3 * IS.fict][j];
                }
            }
        }
    }

    // Y sending
    if ((myrank / px) % 2 != 0) {
        if ((myrank / px) > 0) {
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Nx + 2 * IS.fict; j++) {
                    buff_y[i][j] = m[j][i + IS.fict];
                    buff_y[i + IS.fict][j] = impx[j][i + IS.fict];
                    buff_y[i + 2 * IS.fict][j] = impy[j][i + IS.fict];
                    buff_y[i + 3 * IS.fict][j] = e[j][i + IS.fict];
                }
            }
            send_vec2(buff_y, 4 * IS.fict, IS.Nx + 2 * IS.fict, myrank - px, 1, comm);
        }
        if ((myrank / px) > 0) {
            buff_y = recv_vec2(4 * IS.fict, IS.Nx + 2 * IS.fict, myrank - px, 1, comm, Status);
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Nx + 2 * IS.fict; j++) {
                    m[j][i] = buff_y[i][j];
                    impx[j][i] = buff_y[i + IS.fict][j];
                    impy[j][i] = buff_y[i + 2 * IS.fict][j];
                    e[j][i] = buff_y[i + 3 * IS.fict][j];
                }
            }
        }
        if ((myrank / px) < py - 1) {
            buff_y = recv_vec2(4 * IS.fict, IS.Nx + 2 * IS.fict, myrank + px, 1, comm, Status);
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Nx + 2 * IS.fict; j++) {
                    m[j][IS.Ny + IS.fict + i] = buff_y[i][j];
                    impx[j][IS.Ny + IS.fict + i] = buff_y[i + IS.fict][j];
                    impy[j][IS.Ny + IS.fict + i] = buff_y[i + 2 * IS.fict][j];
                    e[j][IS.Ny + IS.fict + i] = buff_y[i + 3 * IS.fict][j];
                }
            }
        }
        if ((myrank / px) < py - 1) {
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Nx + 2 * IS.fict; j++) {
                    buff_y[i][j] = m[j][IS.Ny + i];
                    buff_y[i + IS.fict][j] = impx[j][IS.Ny + i];
                    buff_y[i + 2 * IS.fict][j] = impy[j][IS.Ny + i];
                    buff_y[i + 3 * IS.fict][j] = e[j][IS.Ny + i];
                }
            }
            send_vec2(buff_y, 4 * IS.fict, IS.Nx + 2 * IS.fict, myrank + px, 1, comm);
        }
    }
    if ((myrank / px) % 2 == 0) {
        if ((myrank / px) < py - 1) {
            buff_y = recv_vec2(4 * IS.fict, IS.Nx + 2 * IS.fict, myrank + px, 1, comm, Status);
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Nx + 2 * IS.fict; j++) {
                    m[j][IS.Ny + IS.fict + i] = buff_y[i][j];
                    impx[j][IS.Ny + IS.fict + i] = buff_y[i + IS.fict][j];
                    impy[j][IS.Ny + IS.fict + i] = buff_y[i + 2 * IS.fict][j];
                    e[j][IS.Ny + IS.fict + i] = buff_y[i + 3 * IS.fict][j];
                }
            }
        }
        if ((myrank / px) < py - 1) {
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Nx + 2 * IS.fict; j++) {
                    buff_y[i][j] = m[j][IS.Ny + i];
                    buff_y[i + IS.fict][j] = impx[j][IS.Ny + i];
                    buff_y[i + 2 * IS.fict][j] = impy[j][IS.Ny + i];
                    buff_y[i + 3 * IS.fict][j] = e[j][IS.Ny + i];
                }
            }
            send_vec2(buff_y, 4 * IS.fict, IS.Nx + 2 * IS.fict, myrank + px, 1, comm);
        }
        if ((myrank / px) > 0) {
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Nx + 2 * IS.fict; j++) {
                    buff_y[i][j] = m[j][i + IS.fict];
                    buff_y[i + IS.fict][j] = impx[j][i + IS.fict];
                    buff_y[i + 2 * IS.fict][j] = impy[j][i + IS.fict];
                    buff_y[i + 3 * IS.fict][j] = e[j][i + IS.fict];
                }
            }
            send_vec2(buff_y, 4 * IS.fict, IS.Nx + 2 * IS.fict, myrank - px, 1, comm);
        }
        if ((myrank / px) > 0) {
            buff_y = recv_vec2(4 * IS.fict, IS.Nx + 2 * IS.fict, myrank - px, 1, comm, Status);
            for (int i = 0; i < IS.fict; i++) {
                for (int j = 0; j < IS.Nx + 2 * IS.fict; j++) {
                    m[j][i] = buff_y[i][j];
                    impx[j][i] = buff_y[i + IS.fict][j];
                    impy[j][i] = buff_y[i + 2 * IS.fict][j];
                    e[j][i] = buff_y[i + 3 * IS.fict][j];
                }
            }
        }
    }
}
*/

void Solve_one_step(InitialState IS, double* t_start, vec1d& xc, vec1d& yc,
    vec2d& P, vec2d& ux, vec2d& uy, vec2d& rho, vec2d& mass, vec2d& Imp_x, vec2d& Imp_y, vec2d& e, std::vector<vec3d>& U, vec2d rk, vec3d& L) {

    double dt = get_dt(IS, xc, yc, mass, Imp_x, Imp_y, e);
    *t_start += dt;

    vec2d m_new(IS.Ny, vec1d(IS.Nx));
    vec2d impx_new(IS.Ny, vec1d(IS.Nx));
    vec2d impy_new(IS.Ny, vec1d(IS.Nx));
    vec2d e_new(IS.Ny, vec1d(IS.Nx));

    // for Kolgan/Radionov
    vec2d m_kolganl(IS.Ny, vec1d(IS.Nx));
    vec2d impx_kolganl(IS.Ny, vec1d(IS.Nx));
    vec2d impy_kolganl(IS.Ny, vec1d(IS.Nx));
    vec2d e_kolganl(IS.Ny, vec1d(IS.Nx));
    vec2d m_kolganr(IS.Ny, vec1d(IS.Nx));
    vec2d impx_kolganr(IS.Ny, vec1d(IS.Nx));
    vec2d impy_kolganr(IS.Ny, vec1d(IS.Nx));
    vec2d e_kolganr(IS.Ny, vec1d(IS.Nx));

    vec2d m_kolgand(IS.Ny, vec1d(IS.Nx));
    vec2d impx_kolgand(IS.Ny, vec1d(IS.Nx));
    vec2d impy_kolgand(IS.Ny, vec1d(IS.Nx));
    vec2d e_kolgand(IS.Ny, vec1d(IS.Nx));
    vec2d m_kolganu(IS.Ny, vec1d(IS.Nx));
    vec2d impx_kolganu(IS.Ny, vec1d(IS.Nx));
    vec2d impy_kolganu(IS.Ny, vec1d(IS.Nx));
    vec2d e_kolganu(IS.Ny, vec1d(IS.Nx));

    vec2d m05(IS.Ny, vec1d(IS.Nx));
    vec2d impx05(IS.Ny, vec1d(IS.Nx));
    vec2d impy05(IS.Ny, vec1d(IS.Nx));
    vec2d e05(IS.Ny, vec1d(IS.Nx));

    double dm_x, dimpx_x, dimpy_x, de_x;
    double dm_y, dimpx_y, dimpy_y, de_y;
    double Fml, Fimpxl, Fimpyl, Fel, Fmr, Fimpxr, Fimpyr, Fer;  //cell x fluxes
    double Fmd, Fimpxd, Fimpyd, Fed, Fmu, Fimpxu, Fimpyu, Feu;  //cell y fluxes

    double ml, impxl, impyl, el, mr, impxr, impyr, er; // left and right cell values of cons
    double md, impxd, impyd, ed, mu, impxu, impyu, eu; // down and up cell values of cons

    if (IS.s_type != 4) {
        if (IS.s_type == 2 || IS.s_type == 3) {
            for (int i = IS.fict; i < IS.Ny - IS.fict + 1; ++i) {
                for (int j = IS.fict; j < IS.Nx - IS.fict + 1; ++j) {
                    Kolgan_left(mass[i][j], mass[i - 1][j], mass[i - 2][j], mass[i + 1][j], xc[j], xc[j - 1], xc[j - 2], xc[j + 1], m_kolganl[i][j], m_kolganr[i][j], IS.mod_type);
                    Kolgan_left(Imp_x[i][j], Imp_x[i - 1][j], Imp_x[i - 2][j], Imp_x[i + 1][j], xc[j], xc[j - 1], xc[j - 2], xc[j + 1], impx_kolganl[i][j], impx_kolganr[i][j], IS.mod_type);
                    Kolgan_left(Imp_y[i][j], Imp_y[i - 1][j], Imp_y[i - 2][j], Imp_y[i + 1][j], xc[j], xc[j - 1], xc[j - 2], xc[j + 1], impy_kolganl[i][j], impy_kolganr[i][j], IS.mod_type);
                    Kolgan_left(e[i][j], e[i - 1][j], e[i - 2][j], e[i + 1][j], xc[j], xc[j - 1], xc[j - 2], xc[j + 1], e_kolganl[i][j], e_kolganr[i][j], IS.mod_type);

                    if (i < IS.Ny - IS.fict) {
                        Kolgan_right(mass[i][j], mass[i - 1][j], mass[i + 1][j], mass[i + 2][j], xc[j], xc[j - 1], xc[j + 1], xc[j + 2], m_kolganl[i + 1][j], m_kolganr[i + 1][j], IS.mod_type);
                        Kolgan_right(Imp_x[i][j], Imp_x[i - 1][j], Imp_x[i + 1][j], Imp_x[i + 2][j], xc[j], xc[j - 1], xc[j + 1], xc[j + 2], impx_kolganl[i + 1][j], impx_kolganr[i + 1][j], IS.mod_type);
                        Kolgan_right(Imp_y[i][j], Imp_y[i - 1][j], Imp_y[i + 1][j], Imp_y[i + 2][j], xc[j], xc[j - 1], xc[j + 1], xc[j + 2], impy_kolganl[i + 1][j], impy_kolganr[i + 1][j], IS.mod_type);
                        Kolgan_right(e[i][j], e[i - 1][j], e[i + 1][j], e[i + 2][j], xc[j], xc[j - 1], xc[j + 1], xc[j + 2], e_kolganl[i + 1][j], e_kolganr[i + 1][j], IS.mod_type);
                    }

                    Kolgan_left(mass[i][j], mass[i][j - 1], mass[i][j - 2], mass[i][j + 1], yc[i], yc[i - 1], yc[i - 2], yc[i + 1], m_kolgand[i][j], m_kolganu[i][j], IS.mod_type);
                    Kolgan_left(Imp_x[i][j], Imp_x[i][j - 1], Imp_x[i][j - 2], Imp_x[i][j + 1], yc[i], yc[i - 1], yc[i - 2], yc[i + 1], impx_kolgand[i][j], impx_kolganu[i][j], IS.mod_type);
                    Kolgan_left(Imp_y[i][j], Imp_y[i][j - 1], Imp_y[i][j - 2], Imp_y[i][j + 1], yc[i], yc[i - 1], yc[i - 2], yc[i + 1], impy_kolgand[i][j], impy_kolganu[i][j], IS.mod_type);
                    Kolgan_left(e[i][j], e[i][j - 1], e[i][j - 2], e[i][j + 1], yc[i], yc[i - 1], yc[i - 2], yc[i + 1], e_kolgand[i][j], e_kolganu[i][j], IS.mod_type);

                    if (j < IS.Nx - IS.fict) {
                        Kolgan_right(mass[i][j], mass[i][j - 1], mass[i][j + 1], mass[i][j + 2], yc[i], yc[i - 1], yc[i + 1], yc[i + 2], m_kolgand[i][j + 1], m_kolganu[i][j + 1], IS.mod_type);
                        Kolgan_right(Imp_x[i][j], Imp_x[i][j - 1], Imp_x[i][j + 1], Imp_x[i][j + 2], yc[i], yc[i - 1], yc[i + 1], yc[i + 2], impx_kolgand[i][j + 1], impx_kolganu[i][j + 1], IS.mod_type);
                        Kolgan_right(Imp_y[i][j], Imp_y[i][j - 1], Imp_y[i][j + 1], Imp_y[i][j + 2], yc[i], yc[i - 1], yc[i + 1], yc[i + 2], impy_kolgand[i][j + 1], impy_kolganu[i][j + 1], IS.mod_type);
                        Kolgan_right(e[i][j], e[i][j - 1], e[i][j + 1], e[i][j + 2], yc[i], yc[i - 1], yc[i + 1], yc[i + 2], e_kolgand[i][j + 1], e_kolganr[i][j + 1], IS.mod_type);
                    }
                }
            }
        }
        if (IS.s_type == 3) {
            for (int i = IS.fict; i < IS.Ny - IS.fict; ++i) {
                for (int j = IS.fict; j < IS.Nx - IS.fict; ++j) {
                    // Predictor
                    cons_flux_y(IS, m_kolganr[i][j], impx_kolganr[i][j], impy_kolganr[i][j], e_kolganr[i][j], Fml, Fimpxl, Fimpyl, Fel);
                    cons_flux_y(IS, m_kolganl[i + 1][j], impx_kolganl[i + 1][j], impy_kolganl[i + 1][j], e_kolganl[i + 1][j], Fmr, Fimpxr, Fimpyr, Fer);
                    cons_flux_x(IS, m_kolganu[i][j], impx_kolganu[i][j], impy_kolganu[i][j], e_kolganu[i][j], Fmd, Fimpxd, Fimpyd, Fed);
                    cons_flux_x(IS, m_kolgand[i][j + 1], impx_kolgand[i][j + 1], impy_kolgand[i][j + 1], e_kolgand[i][j + 1], Fmu, Fimpxu, Fimpyu, Feu);
                    m05[i][j] = mass[i][j] - dt * ((Fmr - Fml) / (xc[j] - xc[j - 1]) + (Fmu - Fmd) / (yc[i] - yc[ - 1]));
                    m05[i][j] = mass[i][j] - dt * ((Fmr - Fml) / (xc[j] - xc[j - 1]) + (Fmu - Fmd) / (yc[i] - yc[ - 1]));
                    impx05[i][j] = Imp_x[i][j] - dt * ((Fimpxr - Fimpxl) / (xc[j] - xc[j - 1]) + (Fimpxu - Fimpxd) / (yc[i] - yc[i - 1]));
                    impy05[i][j] = Imp_y[i][j] - dt * ((Fimpyr - Fimpyl) / (xc[j] - xc[j - 1]) + (Fimpyu - Fimpyd) / (yc[i] - yc[i - 1]));
                    e05[i][j] = e[i][j] - dt * ((Fer - Fel) / (xc[j] - xc[j - 1]) + (Feu - Fed) / (yc[i] - yc[i - 1]));
                }
            }
            for (int i = IS.fict; i < IS.Ny - IS.fict + 1; ++i) {
                for (int j = IS.fict; j < IS.Nx - IS.fict + 1; ++j) {
                    // Corrector x
                    dm_x = (m_kolganl[i + 1][j] - m_kolganr[i][j]);
                    dimpx_x = (impx_kolganl[i + 1][j] - impx_kolganr[i][j]);
                    dimpy_x = (impy_kolganl[i + 1][j] - impy_kolganr[i][j]);
                    de_x = (e_kolganl[i + 1][j] - e_kolganr[i][j]);

                    m_kolganr[i][j] = 0.5 * (mass[i][j] + m05[i][j]) - 0.5 * dm_x;
                    impx_kolganr[i][j] = 0.5 * (Imp_x[i][j] + impx05[i][j]) - 0.5 * dimpx_x;
                    impy_kolganr[i][j] = 0.5 * (Imp_y[i][j] + impy05[i][j]) - 0.5 * dimpy_x;
                    e_kolganr[i][j] = 0.5 * (e[i][j] + e05[i][j]) - 0.5 * de_x;
                    m_kolganl[i + 1][j] = 0.5 * (mass[i][j] + m05[i][j]) + 0.5 * dm_x;
                    impx_kolganl[i + 1][j] = 0.5 * (Imp_x[i][j] + impx05[i][j]) + 0.5 * dimpx_x;
                    impy_kolganl[i + 1][j] = 0.5 * (Imp_y[i][j] + impy05[i][j]) + 0.5 * dimpy_x;
                    e_kolganl[i + 1][j] = 0.5 * (e[i][j] + e05[i][j]) + 0.5 * de_x;

                    // Corrector y
                    dm_y = (m_kolgand[i][j + 1] - m_kolganu[i][j]);
                    dimpx_y = (impx_kolgand[i][j + 1] - impx_kolganu[i][j]);
                    dimpy_y = (impy_kolgand[i][j + 1] - impy_kolganu[i][j]);
                    de_y = (e_kolgand[i][j + 1] - e_kolganu[i][j]);

                    m_kolganu[i][j] = 0.5 * (mass[i][j] + m05[i][j]) - 0.5 * dm_y;
                    impx_kolganu[i][j] = 0.5 * (Imp_x[i][j] + impx05[i][j]) - 0.5 * dimpx_y;
                    impy_kolganu[i][j] = 0.5 * (Imp_y[i][j] + impy05[i][j]) - 0.5 * dimpy_y;
                    e_kolganu[i][j] = 0.5 * (e[i][j] + e05[i][j]) - 0.5 * de_y;
                    m_kolgand[i][j + 1] = 0.5 * (mass[i][j] + m05[i][j]) + 0.5 * dm_y;
                    impx_kolgand[i][j + 1] = 0.5 * (Imp_x[i][j] + impx05[i][j]) + 0.5 * dimpx_y;
                    impy_kolgand[i][j + 1] = 0.5 * (Imp_y[i][j] + impy05[i][j]) + 0.5 * dimpy_y;
                    e_kolgand[i][j + 1] = 0.5 * (e[i][j] + e05[i][j]) + 0.5 * de_y;
                }
            }
        }
        for (int i = IS.fict; i < IS.Ny - IS.fict; ++i) {
            for (int j = IS.fict; j < IS.Nx - IS.fict; ++j) {
                if (IS.s_type == 1) { // X
                    ml = mass[i][j - 1];
                    impxl = Imp_x[i][j - 1];
                    impyl = Imp_y[i][j - 1];
                    el = e[i][j - 1];
                    mr = mass[i][j];
                    impxr = Imp_x[i][j];
                    impyr = Imp_y[i][j];
                    er = e[i][j];
                }
                else if (IS.s_type == 2 || IS.s_type == 3) {
                    ml = m_kolganl[i][j];
                    mr = m_kolganr[i][j];
                    impxl = impx_kolganl[i][j];
                    impxr = impx_kolganr[i][j];
                    impyl = impy_kolganl[i][j];
                    impyr = impy_kolganr[i][j];
                    el = e_kolganl[i][j];
                    er = e_kolganr[i][j];
                }
                Godunov_solver_x(IS, ml, impxl, impyl, el, mr, impxr, impyr, er, Fml, Fimpxl, Fimpyl, Fel);

                if (IS.s_type == 1) { // X
                    ml = mass[i][j];
                    impxl = Imp_x[i][j];
                    impyl = Imp_y[i][j];
                    el = e[i][j];
                    mr = mass[i][j + 1];
                    impxr = Imp_x[i][j + 1];
                    impyr = Imp_y[i][j + 1];
                    er = e[i][j + 1];
                }
                else if (IS.s_type == 2 || IS.s_type == 3) {
                    ml = m_kolganl[i][j + 1];
                    mr = m_kolganr[i][j + 1];
                    impxl = impx_kolganl[i][j + 1];
                    impxr = impx_kolganr[i][j + 1];
                    impyl = impy_kolganl[i][j + 1];
                    impyr = impy_kolganr[i][j + 1];
                    el = e_kolganl[i][j + 1];
                    er = e_kolganr[i][j + 1];
                }
                Godunov_solver_x(IS, ml, impxl, impyl, el, mr, impxr, impyr, er, Fmr, Fimpxr, Fimpyr, Fer);

                if (IS.s_type == 1) {  // Y
                    md = mass[i - 1][j];
                    impxd = Imp_x[i - 1][j];
                    impyd = Imp_y[i - 1][j];
                    ed = e[i - 1][j];
                    mu = mass[i][j];
                    impxu = Imp_x[i][j];
                    impyu = Imp_y[i][j];
                    eu = e[i][j];
                }
                else if (IS.s_type == 2 || IS.s_type == 3) {
                    md = m_kolgand[i][j];
                    mu = m_kolganu[i][j];
                    impxd = impx_kolgand[i][j];
                    impxu = impx_kolganu[i][j];
                    impyd = impy_kolgand[i][j];
                    impyu = impy_kolganu[i][j];
                    ed = e_kolgand[i][j];
                    eu = e_kolganu[i][j];
                }
                Godunov_solver_y(IS, md, impxd, impyd, ed, mu, impxu, impyu, eu, Fmd, Fimpxd, Fimpyd, Fed);

                if (IS.s_type == 1) { // Y
                    md = mass[i][j];
                    impxd = Imp_x[i][j];
                    impyd = Imp_y[i][j];
                    ed = e[i][j];
                    mu = mass[i + 1][j];
                    impxu = Imp_x[i + 1][j];
                    impyu = Imp_y[i + 1][j];
                    eu = e[i + 1][j];
                }
                else if (IS.s_type == 2 || IS.s_type == 3) {
                    md = m_kolgand[i + 1][j];
                    mu = m_kolganu[i + 1][j];
                    impxd = impx_kolgand[i + 1][j];
                    impxu = impx_kolganu[i + 1][j];
                    impyd = impy_kolgand[i + 1][j];
                    impyu = impy_kolganu[i + 1][j];
                    ed = e_kolgand[i + 1][j];
                    eu = e_kolganu[i + 1][j];
                }
                Godunov_solver_y(IS, md, impxd, impyd, ed, mu, impxu, impyu, eu, Fmu, Fimpxu, Fimpyu, Feu);

                if (IS.s_type != 0) {

                    m_new[i][j] = mass[i][j] - dt * ((Fmr - Fml) / IS.hx + (Fmu - Fmd) / IS.hy);
                    impx_new[i][j] = Imp_x[i][j] - dt * ((Fimpxr - Fimpxl) / IS.hx + (Fimpxu - Fimpxd) / IS.hy);
                    impy_new[i][j] = Imp_y[i][j] - dt * ((Fimpyr - Fimpyl) / IS.hx + (Fimpyu - Fimpyd) / IS.hy);
                    e_new[i][j] = e[i][j] - dt * ((Fer - Fel) / IS.hx + (Feu - Fed) / IS.hy);
                }
            }
        }
        Boundary_x(IS, m_new, impx_new, impy_new, e_new);
        Boundary_y(IS, m_new, impx_new, impy_new, e_new);


        for (int i = 0; i < IS.Ny; i++) {
            for (int j = 0; j < IS.Nx; j++) {
                Imp_x[i][j] = impx_new[i][j];
                Imp_y[i][j] = impy_new[i][j];
                mass[i][j] = m_new[i][j];
                e[i][j] = e_new[i][j];
                cons_to_noncons(IS, P[i][j], ux[i][j], uy[i][j], rho[i][j], mass[i][j], Imp_x[i][j], Imp_y[i][j], e[i][j]);
            }
        }
    }
    else if (IS.s_type == 4) {
        for (int i = 0; i < IS.Ny; i++) {
            for (int j = 0; j < IS.Nx; j++) {
                U[0][0][i][j] = mass[i][j];
                U[0][1][i][j] = Imp_x[i][j];
                U[0][2][i][j] = Imp_y[i][j];
                U[0][3][i][j] = e[i][j];
            }
        }

        for (int s = 1; s < 4; s++) {
            WENO5(IS, U[s - 1], L, IS.hx, IS.hy);
            for (int k = 0; k < 4; k++) {
                for (int i = IS.fict; i < IS.Ny - IS.fict; i++) {
                    for (int j = IS.fict; j < IS.Nx - IS.fict; j++) {
                        U[s][k][i][j] = rk[s][0] * U[0][k][i][j] + rk[s][1] * U[1][k][i][j] + rk[s][2] * U[2][k][i][j] - rk[s][3] * dt * L[k][i][j];
                    }
                }
            }
            Boundary_x(IS, U[s][0], U[s][1], U[s][2], U[s][3]);
            Boundary_y(IS, U[s][0], U[s][1], U[s][2], U[s][3]);
        }

        for (int i = 0; i < IS.Ny; i++) {
            for (int j = 0; j < IS.Nx; j++) {
                mass[i][j] = U[3][0][i][j];
                Imp_x[i][j] = U[3][1][i][j];
                Imp_y[i][j] = U[3][2][i][j];
                e[i][j] = U[3][3][i][j];
                cons_to_noncons(IS, P[i][j], ux[i][j], uy[i][j], rho[i][j], mass[i][j], Imp_x[i][j], Imp_y[i][j], e[i][j]);
            }
        }
    }
    else {
        std::cerr << "Invalid s_type" << std::endl;
    }
}


void Solve(InitialState& IS, std::string filename) {
    vec1d x(IS.Nx + 1);
    vec1d xc(IS.Nx);
    vec1d y(IS.Ny + 1);
    vec1d yc(IS.Ny);

    vec2d P(IS.Ny, vec1d(IS.Nx));
    vec2d ux(IS.Ny, vec1d(IS.Nx));
    vec2d uy(IS.Ny, vec1d(IS.Nx));
    vec2d rho(IS.Ny, vec1d(IS.Nx));
    vec2d mass(IS.Ny, vec1d(IS.Nx));
    vec2d Imp_x(IS.Ny, vec1d(IS.Nx));
    vec2d Imp_y(IS.Ny, vec1d(IS.Nx));
    vec2d e(IS.Ny, vec1d(IS.Nx));
    vec2d ei(IS.Ny, vec1d(IS.Nx));

    // for WENO
    std::vector<vec3d> U(4, vec3d(4, vec2d(IS.Ny, vec1d(IS.Nx))));
    std::vector<vec1d> rk = { {1., 0., 0., 0.}, {1., 0., 0., 1.}, {0.75, 0.25, 0., 0.25}, {1. / 3., 0., 2. / 3., 2. / 3.} };
    vec3d L(4, vec2d(IS.Ny, vec1d(IS.Nx)));

    int step = 0;
    double t_start = 0.;

    x_grid(IS, xc, x);
    y_grid(IS, yc, y);
    Init(IS, xc, yc, P, ux, uy, rho, mass, Imp_x, Imp_y, e);
    Boundary_x(IS, mass, Imp_x, Imp_y, e);
    Boundary_y(IS, mass, Imp_x, Imp_y, e);
    calc_ei(IS, P, rho, ei);

    if (fs::exists(filename)) {
        fs::remove_all(filename);
    }
    writeCSV(filename, IS, step, xc, yc, ux, uy, P, rho, t_start);

    while (t_start < IS.t_end) {
        Solve_one_step(IS, &t_start, xc, yc, P, ux, uy, rho, mass, Imp_x, Imp_y, e, U, rk, L);
        step++;
        if (step % IS.write_interval == 0) {
            calc_ei(IS, P, rho, ei);
            writeCSV(filename, IS, step, xc, yc, ux, uy, P, rho, t_start);
        }
    }
}


void Solve_p(InitialState& IS, std::string filename, int size, int myrank) {
    vec1d x(IS.Nx + 1);
    vec1d xc(IS.Nx);
    vec1d y(IS.Ny + 1);
    vec1d yc(IS.Ny);

    vec2d P(IS.Ny, vec1d(IS.Nx));
    vec2d ux(IS.Ny, vec1d(IS.Nx));
    vec2d uy(IS.Ny, vec1d(IS.Nx));
    vec2d rho(IS.Ny, vec1d(IS.Nx));
    vec2d mass(IS.Ny, vec1d(IS.Nx));
    vec2d Imp_x(IS.Ny, vec1d(IS.Nx));
    vec2d Imp_y(IS.Ny, vec1d(IS.Nx));
    vec2d e(IS.Ny, vec1d(IS.Nx));
    vec2d ei(IS.Ny, vec1d(IS.Nx));

    // for WENO
    std::vector<vec3d> U(4, vec3d(4, vec2d(IS.Ny, vec1d(IS.Nx))));
    std::vector<vec1d> rk = { {1., 0., 0., 0.}, {1., 0., 0., 1.}, {0.75, 0.25, 0., 0.25}, {1. / 3., 0., 2. / 3., 2. / 3.} };
    vec3d L(4, vec2d(IS.Ny, vec1d(IS.Nx)));

    int step = 0;
    double t_start = 0.;

    x_grid(IS, xc, x);
    y_grid(IS, yc, y);
    Init(IS, xc, yc, P, ux, uy, rho, mass, Imp_x, Imp_y, e);
    Boundary_x(IS, mass, Imp_x, Imp_y, e);
    Boundary_y(IS, mass, Imp_x, Imp_y, e);
    calc_ei(IS, P, rho, ei);

    if (fs::exists(filename)) {
        fs::remove_all(filename);
    }
    writeCSV(filename, IS, step, xc, yc, ux, uy, P, rho, t_start);

    while (t_start < IS.t_end) {
        Solve_one_step(IS, &t_start, xc, yc, P, ux, uy, rho, mass, Imp_x, Imp_y, e, U, rk, L);
        step++;
        if (step % IS.write_interval == 0) {
            calc_ei(IS, P, rho, ei);
            writeCSV(filename, IS, step, xc, yc, ux, uy, P, rho, t_start);
        }
    }
}

int main(int argc, char* argv[]) {
    int myrank, size;
    // MPI_Init(&argc, &argv);
    // MPI_Comm_size(MPI_COMM_WORLD, &size);
    // MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    myrank = 0;
    double begin, end;

    InitialState IS;
    std::string filename = "Godunov_2D";

    /*
    read_params(IS, "G.txt");
    MPI_Barrier(MPI_COMM_WORLD);
    begin = MPI_Wtime();
    Multiproc_solve(IS, filename, myrank, size, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    if(myrank == 0) printf("t_start parallel: %.5f\n\n", end - begin);

    read_params(IS, "GK.txt");
    MPI_Barrier(MPI_COMM_WORLD);
    begin = MPI_Wtime();
    Multiproc_solve(IS, filename, myrank, size, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    if(myrank == 0) printf("t_start parallel: %.5f\n\n", end - begin);

    read_params(IS, "GKR.txt");
    MPI_Barrier(MPI_COMM_WORLD);
    begin = MPI_Wtime();
    Multiproc_solve(IS, filename, myrank, size, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    if(myrank == 0) printf("t_start parallel: %.5f\n\n", end - begin);

    read_params(IS, "W.txt");
    MPI_Barrier(MPI_COMM_WORLD);
    begin = MPI_Wtime();
    Multiproc_solve(IS, filename, myrank, size, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    if(myrank == 0) printf("t_start parallel: %.5f\n\n", end - begin);
    */

    if (myrank == 0) {
        // begin = MPI_Wtime();
	    init_params(IS);
        // read_params(IS, "W.txt");
        Solve(IS, filename);
        // end = MPI_Wtime();
        // printf("time sequently: %.5f\n\n", end - begin);
    }

    // MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Finalize();
    return 0;
}