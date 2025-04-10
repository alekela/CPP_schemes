#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unordered_map>


//namespace fs = std::experimental::filesystem;
namespace fs = std::filesystem;
typedef std::vector<double> vec1d;
typedef std::vector<vec1d> vec2d;
typedef std::vector<vec2d> vec3d;


struct InitialState {
    double x_start;
    double x_end;
    double y_start;
    double y_end;
    double gamma;
    
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


InitialState change_params(InitialState IS, int Nx, int Ny, double x_start, double x_end, double y_start, double y_end, int myrank, int px, int py){
    InitialState result;
    result = IS;
    result.Nx = Nx;
    result.Ny = Ny;
    result.x_start = x_start;
    result.y_start = y_start;
    result.x_end = x_end;
    result.y_end = y_end;
    if(myrank % px != 0) result.b_left = -1;
    if(myrank % px != px - 1) result.b_right = -1;
    if(myrank / px != 0) result.b_down = -1;
    if(myrank / px != py - 1) result.b_up = -1;
    return result;
}

//constant parameters
void init_params(InitialState& IS) {
    IS.gamma = 1.4;
    IS.x_start = 0.;
    IS.x_end = 1.;
    IS.y_start = 0.;
    IS.y_end = 1.;

    IS.Nx = 50;
    IS.Ny = 50;
    IS.t_end = 2;
    IS.CFL = 0.5;  

    IS.b_left = 1;
    IS.b_right = 1;
    IS.b_up = 0;
    IS.b_down = 0;
    IS.initial = 2;
    IS.write_interval = 10;
    IS.max_iter = 1000;
    
    IS.s_type = 1;
    IS.mod_type = 3;

    if (IS.s_type == 1) IS.fict = 1;
    else if (IS.s_type == 2 || IS.s_type == 3) IS.fict = 2;
    else if (IS.s_type == 4) IS.fict = 3;
    else {
        std::cerr << "Invalid s_type in parameters. Change to 1 (Godunov)" << std::endl;
        IS.fict = 1;
    }
}

//parameters from file
void read_params(InitialState& IS, std::string file_name) {
    std::ifstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << file_name << std::endl;
        return;
    }

    std::unordered_map<std::string, std::string> text;
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, ' ') && std::getline(iss, value, ';')) {
            text[key] = value;
        }
    }
    file.close();
                                                
    IS.gamma = std::stod(text["g"]);
    IS.x_start = std::stod(text["x_start"]);
    IS.x_end = std::stod(text["x_end"]);
    IS.y_start = std::stod(text["y_start"]);
    IS.y_end = std::stod(text["y_end"]);

    IS.Nx = std::stoi(text["Nx"]);
    IS.Ny = std::stoi(text["Ny"]);
    IS.t_end = std::stod(text["t_end"]);
    IS.CFL = std::stod(text["CFL"]);

    IS.b_left = std::stoi(text["b_left"]);
    IS.b_right = std::stoi(text["b_right"]);
    IS.b_up = std::stoi(text["b_up"]);
    IS.b_down = std::stoi(text["b_down"]);

    IS.initial = std::stoi(text["initial"]);
    IS.write_interval = std::stoi(text["write_interval"]);
    IS.max_iter = std::stoi(text["max_iter"]);

    IS.s_type = std::stoi(text["s_type"]);
    IS.mod_type = std::stoi(text["mod_type"]);

    if (IS.mod_type > 4) {
        std::cerr << "Invalid mod_type in parameters. Change to default setting - 3 (Osher minmod)" << std::endl;
        IS.mod_type = 3;
    }

    if (IS.s_type == 1) IS.fict = 1;
    else if (IS.s_type == 2 || IS.s_type == 3) IS.fict = 2;
    else if (IS.s_type == 4) IS.fict = 3;
    else {
        std::cerr << "Invalid s_type in parameters. Change to 1 (Godunov)" << std::endl;
        IS.fict = 1;
    }
}

//convert conservative variables to nonconservative
void cons_to_noncons(InitialState& IS, double& p, double& vx, double& vy, double& r, double& m, double& impx, double& impy, double& e) {
    p = (IS.gamma - 1.0) * (e - 0.5 * (pow(impx, 2.0) + pow(impy, 2.0)) / m);
    if(p < 0.){
	    p = 1.e-6;
    }
    vx = impx / m;
    vy = impy / m;
    r = m;
}

void cons_to_noncons(InitialState& IS, vec2d& p, vec2d& vx, vec2d& vy, vec2d& r, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    for(int i = 0; i < IS.Nx + 2 * IS.fict; ++i) {
        for(int j = 0; j < IS.Ny + 2 * IS.fict; ++j) {
            cons_to_noncons(IS, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
        }
    }
}

//convert nonconservative variables to conservative
void noncons_to_cons(InitialState& IS, double& p, double& vx, double& vy, double& r, double& m, double& impx, double& impy, double& e) {
    m = r;
    impx = r * vx;
    impy = r * vy;
    e = 0.5 * r * (pow(vx, 2.0) + pow(vy, 2.0)) + p / (IS.gamma - 1.0);
}

//calculate internal energy
void calc_ei(InitialState& IS, vec2d& p, vec2d& r, vec2d& ei){
    for(int i = 0; i < IS.Nx + 2 * IS.fict; ++i){
        for(int j = 0; j < IS.Ny + 2 * IS.fict; ++j){
            ei[i][j] = p[i][j] / ((IS.gamma - 1.) * r[i][j]);
        }
    }
}

//calculate sound velocity
double calc_c(double gamma, double p, double r) {
    return std::sqrt(gamma * p / r);
}

//build x grid
void x_grid(InitialState& IS, vec1d& xc, vec1d& x) {
    double dx;
    dx = (IS.x_end - IS.x_start) / IS.Nx;
    for (int i = 0; i < IS.Nx + 1 + 2*IS.fict; ++i) {
        x[i] = IS.x_start + (i - IS.fict) * dx;
    }
    for (int i = 0; i < IS.Nx + 2 * IS.fict; ++i) {
        xc[i] = 0.5 * (x[i] + x[i + 1]);
    }
}

//build y grid
void y_grid(InitialState& IS, vec1d& yc, vec1d& y) {
    double dy;
    dy = (IS.y_end - IS.y_start) / IS.Ny;
    for (int i = 0; i < IS.Ny + 1 + 2*IS.fict; ++i) {
        y[i] = IS.y_start + (i - IS.fict) * dy;
    }
    for (int i = 0; i < IS.Ny + 2 * IS.fict; ++i) {
        yc[i] = 0.5 * (y[i] + y[i + 1]);
    }
}

//boundary X
void Boundary_x(InitialState& IS, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    for(int i = 1; i <= IS.fict; ++i){
        for(int j = 0; j < IS.Ny + 2 * IS.fict; ++j){
            if (IS.b_left == 0) {   //wall
                m[IS.fict - i][j] = m[IS.fict + i - 1][j];
                impx[IS.fict - i][j] = -impx[IS.fict + i - 1][j];
                impy[IS.fict - i][j] = impy[IS.fict + i - 1][j];
                e[IS.fict - i][j] = e[IS.fict + i - 1][j];
            }
            else if (IS.b_left == 1) {   //free boundery
                m[IS.fict - i][j] = m[IS.fict + i - 1][j];
                impx[IS.fict - i][j] = impx[IS.fict + i - 1][j];
                impy[IS.fict - i][j] = impy[IS.fict + i - 1][j];
                e[IS.fict - i][j] = e[IS.fict + i - 1][j];
            }
            if (IS.b_right == 0) {   //wall
                m[IS.Nx + IS.fict + i - 1][j] = m[IS.Nx + IS.fict - i][j];
                impx[IS.Nx + IS.fict + i - 1][j] = -impx[IS.Nx + IS.fict - i][j];
                impy[IS.Nx + IS.fict + i - 1][j] = impy[IS.Nx + IS.fict - i][j];
                e[IS.Nx + IS.fict + i - 1][j] = e[IS.Nx + IS.fict - i][j];
            }
            else if (IS.b_right == 1) {   //free boundery
                m[IS.Nx + IS.fict + i - 1][j] = m[IS.Nx + IS.fict - i][j];
                impx[IS.Nx + IS.fict + i - 1][j] = impx[IS.Nx + IS.fict - i][j];
                impy[IS.Nx + IS.fict + i - 1][j] = impy[IS.Nx + IS.fict - i][j];
                e[IS.Nx + IS.fict + i - 1][j] = e[IS.Nx + IS.fict - i][j];
            }
        }
    }
}

//boundary Y
void Boundary_y(InitialState& IS, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    for(int j = 1; j <= IS.fict; ++j){
        for(int i = 0; i < IS.Nx + 2 * IS.fict; ++i){
            if (IS.b_down == 0) {   //wall
                m[i][IS.fict - j] = m[i][IS.fict + j - 1];
                impx[i][IS.fict - j] = impx[i][IS.fict + j - 1];
                impy[i][IS.fict - j] = -impy[i][IS.fict + j - 1];
                e[i][IS.fict - j] = e[i][IS.fict + j - 1];
            }
            else if (IS.b_down == 1) {   //free boundery
                m[i][IS.fict - j] = m[i][IS.fict + j - 1];
                impx[i][IS.fict - j] = impx[i][IS.fict + j - 1];
                impy[i][IS.fict - j] = impy[i][IS.fict + j - 1];
                e[i][IS.fict - j] = e[i][IS.fict + j - 1];
            }
            if (IS.b_up == 0) {   //wall
                m[i][IS.Ny + IS.fict + j - 1] = m[i][IS.Ny + IS.fict - j];
                impx[i][IS.Ny + IS.fict + j - 1] = impx[i][IS.Ny + IS.fict - j];
                impy[i][IS.Ny + IS.fict + j - 1] = -impy[i][IS.Ny + IS.fict - j];
                e[i][IS.Ny + IS.fict + j - 1] = e[i][IS.Ny + IS.fict - j];
            }
            else if (IS.b_up == 1) {   //free boundery
                m[i][IS.Ny + IS.fict + j - 1] = m[i][IS.Ny + IS.fict - j];
                impx[i][IS.Ny + IS.fict + j - 1] = impx[i][IS.Ny + IS.fict - j];
                impy[i][IS.Ny + IS.fict + j - 1] = impy[i][IS.Ny + IS.fict - j];
                e[i][IS.Ny + IS.fict + j - 1] = e[i][IS.Ny + IS.fict - j];
            }
        }
    }
}

//caluclate time step
double get_dt(InitialState& IS, vec1d& x, vec1d& y, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    double new_step = 1000000;
    double p, vx, vy, r, c;
    double c_step;
    double CFL = IS.CFL;
    for (int i = IS.fict; i < IS.Nx + IS.fict; ++i) {
        for (int j = IS.fict; j < IS.Ny + IS.fict; ++j) {
            cons_to_noncons(IS, p, vx, vy, r, m[i][j], impx[i][j], impy[i][j], e[i][j]);
            c = calc_c(IS.gamma, p, r);
            c_step = std::min(CFL * (x[i + 1] - x[i]) / (std::fabs(vx) + c), CFL * (y[j + 1] - y[j]) / (std::fabs(vy) + c));
            if (c_step < new_step) {
                new_step = c_step;
            }
        }
    }
    return new_step;
}

void exchange_dt(InitialState& IS, double& dt, vec1d& x, vec1d& y, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e, int myrank, int size, MPI_Comm comm){
    MPI_Status Status;
    if (myrank != 0){
        MPI_Send(&dt, 1, MPI_DOUBLE, 0, myrank, comm);
        MPI_Recv(&dt, 1, MPI_DOUBLE, 0, myrank, comm, &Status);
    }
    if(myrank == 0) {
        double temp;
        for(int i = 1; i < size; ++i){
            MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, comm, &Status);
            if(temp < dt) dt = temp;
        }
        for(int i = 1; i < size; ++i){
            MPI_Send(&dt, 1, MPI_DOUBLE, i, i, comm);
        }
    }
}


void Soda_x(InitialState& IS, vec1d& xc, vec1d& yc, vec2d& p, vec2d& vx, vec2d& vy, vec2d& r, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    double p1 = 1.0;
    double r1 = 1.0;
    double p2 = 0.1;
    double r2 = 0.125;

    for(int i = 0; i < IS.Nx + 2 * IS.fict; ++i){
        for(int j = 0; j < IS.Ny + 2 * IS.fict; ++j){
            if(xc[i] <= 0.5){
                r[i][j] = r1;
                p[i][j] = p1;
                vx[i][j] = 0.;
                vy[i][j] = 0.;
            }
            else{
                r[i][j] = r2;
                p[i][j] = p2;
                vx[i][j] = 0.;
                vy[i][j] = 0.;
            }
            noncons_to_cons(IS, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
        }
    }
}

void Soda_y(InitialState& IS, vec1d& xc, vec1d& yc, vec2d& p, vec2d& vx, vec2d& vy, vec2d& r, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    double p1 = 1.0;
    double r1 = 1.0;
    double p2 = 0.1;
    double r2 = 0.125;

    for(int i = 0; i < IS.Nx + 2 * IS.fict; ++i){
        for(int j = 0; j < IS.Ny + 2 * IS.fict; ++j){
            if(yc[j] <= 0.5){
                r[i][j] = r1;
                p[i][j] = p1;
                vx[i][j] = 0.;
                vy[i][j] = 0.;
            }
            else{
                r[i][j] = r2;
                p[i][j] = p2;
                vx[i][j] = 0.;
                vy[i][j] = 0.;
            }
            noncons_to_cons(IS, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
        }
    }
}

void Gelmgolc(InitialState& IS, vec1d& xc, vec1d& yc, vec2d& p, vec2d& vx, vec2d& vy, vec2d& r, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
    double p1 = 1.0;
    double r1 = 1.0;
    double vx1 = 0.;
    double vy1 = 0.5;
    double p2 = 1.;
    double r2 = 2.;
    double vx2 = 0.;
    double vy2 = -0.5;

    for(int i = 0; i < IS.Nx + 2 * IS.fict; ++i){
        for(int j = 0; j < IS.Ny + 2 * IS.fict; ++j){
            if(xc[i] <= 0.5){
                r[i][j] = r1;
                p[i][j] = p1;
                vx[i][j] = vx1;
                vy[i][j] = vy1;
            }
            else{
                r[i][j] = r2;
                p[i][j] = p2;
                vx[i][j] = vx2;
                vy[i][j] = vy2;
            }
            noncons_to_cons(IS, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
        }
    }
}

void Bubble(InitialState& IS, vec1d& xc, vec1d& yc, vec2d& p, vec2d& vx, vec2d& vy, vec2d& r, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
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

    for (int i = 0; i < IS.Nx + 2 * IS.fict; ++i) {
        for (int j = 0; j < IS.Ny + 2 * IS.fict; ++j) {
            if (pow((xc[i] - 0.5), 2.) + pow((yc[j] - 0.5), 2.) < pow(R, 2.)) {
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
            noncons_to_cons(IS, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
        }
    }
}

void Explosion(InitialState& IS, vec1d& xc, vec1d& yc, vec2d& p, vec2d& vx, vec2d& vy, vec2d& r, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e) {
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

    for (int i = 0; i < IS.Nx + 2 * IS.fict; ++i) {
        for (int j = 0; j < IS.Ny + 2 * IS.fict; ++j) {
            if (pow((xc[i] - 0.5), 2.) + pow((yc[j] - 0.5), 2.) < pow(R, 2.)) {
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
            noncons_to_cons(IS, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
        }
    }
}

void Init(InitialState& IS, vec1d& xc, vec1d& yc, vec2d& p, vec2d& vx, vec2d& vy, vec2d& r, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e){
    switch (IS.initial){
    case 1:
        Bubble(IS, xc, yc, p, vx, vy, r, m, impx, impy, e);
        break;
    case 2:
        Soda_x(IS, xc, yc, p, vx, vy, r, m, impx, impy, e);
        break;
    case 3:
        Explosion(IS, xc, yc, p, vx, vy, r, m, impx, impy, e);
        break;
    case 4:
        Gelmgolc(IS, xc, yc, p, vx, vy, r, m, impx, impy, e);
        break;
    default:
        Bubble(IS, xc, yc, p, vx, vy, r, m, impx, impy, e);
        break;
    }
}

//save results in the directory out_dir
void write_out(InitialState& IS, vec1d& xc, vec1d& yc, vec2d& p, vec2d& vx, vec2d& vy, vec2d& r, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e, vec2d& ei, int _n, double _time, std::string out_dir)
{
    if (!fs::exists(out_dir)) {
        fs::create_directory(out_dir);
    }
    std::string file_path = out_dir + "/" + std::to_string(_n) + ".txt";
    std::ofstream csv;
    csv.open(file_path);
    csv << _time << std::endl;
    csv << "xc;yc;vx;vy;rho;m;p;impx;impy;e;ei\n";
    for(int i = IS.fict; i < IS.Nx + IS.fict; ++i){
        for (int j = IS.fict; j < IS.Ny + IS.fict; ++j) {
            csv << xc[i] << ";" << yc[j] << ";" << vx[i][j] << ";" << vy[i][j] << ";" << r[i][j] << ";" << m[i][j] << ";" << p[i][j] << ";" << impx[i][j] << ";" << impy[i][j] << ";" << e[i][j] << ';' << ei[i][j] << "\n";
        }
    }
    csv.close();
}

void write_out_p(InitialState& IS, vec1d& xc, vec1d& yc, vec2d& p, vec2d& vx, vec2d& vy, vec2d& r, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e, vec2d& ei, int _n, double _time, \
                 std::string out_dir, int myrank, int size, MPI_Comm comm)
{
    if (myrank == 0) {
        if (!fs::exists(out_dir)) {
            fs::create_directory(out_dir);
        }
    }
    MPI_Barrier(comm);

    std::string file_path = out_dir + "/" + "Iter=" + std::to_string(_n) + ".csv";
    std :: ostringstream buffer;
    
    if (myrank == 0) {
        buffer << (std::to_string(_time) + "\n") << "xc;yc;vx;vy;rho;m;p;impx;impy;e;ei\n";
    }
    for(int i = IS.fict; i < IS.Nx + IS.fict; ++i) {
        for(int j = IS.fict; j < IS.Ny + IS.fict; ++j) {
            buffer << xc[i] << ";" << yc[j] << ";" << vx[i][j] << ";" << vy[i][j] << ";" << r[i][j] << ";" << m[i][j] << ";" << p[i][j] \
                    << ";" << impx[i][j] << ";" << impy[i][j] << ";" << e[i][j] << ';' << ei[i][j] << "\n";
        }
    }
    std::string local_data = buffer.str();
    int local_size = local_data.size();

    int offset = 0;
    MPI_Exscan(&local_size, &offset, 1, MPI_INT, MPI_SUM, comm);
    MPI_File fh;
    MPI_File_open(comm, file_path.c_str(),MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write_at_all(fh, offset, local_data.c_str(), local_size, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
}


double Newton_find_P(InitialState& IS, double pl, double vl, double rl, double cl, double pr, double vr, double rr, double cr) {
    double p_lin; //linear approximation
    double p_min, p_max; //min max pressure
    double p_ratio;
    double p_ratio_max = 2.0;
    double tl, tr;
    double g1 = (IS.gamma - 1.0) / 2. / IS.gamma;
    double g2 = (IS.gamma - 1.0) / (IS.gamma + 1.0);
    double g3 = 2.0 / (IS.gamma + 1.0);

    double eps = 1.e-6;

    p_lin = std::max(eps, 0.5 * (pl + pr) + 0.5 * (vl - vr) * 0.25 * (rl + rr) * (cl + cr)); //Toro p.299 (9.20, 9.21)
    p_min = std::min(pl, pr);
    p_max = std::max(pl, pr);
    p_ratio = p_max / p_min;

    if ((p_ratio <= p_ratio_max) && ((p_min < p_lin && p_lin < p_max))) {
        return p_lin;
    }
    else {
        if (p_lin < p_min) {
            //two-rarefaction
            return pow(((cl + cr - 0.5 * (IS.gamma - 1.0) * (vr - vl)) / (cl / pow(pl, g1) + cr / pow(pr, g1))), 1. / g1); //Toro p.301 (9.32)
        }
        else {
            // two-shock approximation
            tl = std::sqrt(g3 / rl / (g2 * pl + p_lin));
            tr = std::sqrt(g3 / rr / (g2 * pr + p_lin));
            return (tl * pl + tr * pr - (vr - vl)) / (tl + tr); //Toro p.303 (9.42)
        }
    }
}


void calc_flux(InitialState& IS, double& F, double& FD, double p_curr, double p, double v, double r, double c) {
    double g1 = (IS.gamma - 1.0) / 2.0 / IS.gamma;
    double g2 = 0.5 * (IS.gamma + 1.0) / IS.gamma;
    double g4 = 2.0 / (IS.gamma - 1.0);
    double g5 = 2.0 / (IS.gamma + 1.0);
    double g6 = (IS.gamma - 1.0) / (IS.gamma + 1.0);
    double p_ratio = p_curr / p;
    double qrt;

    if (p_curr <= p) {
        //rarefaction
        F = g4 * c * (pow(p_ratio, g1) - 1.0); //Toro p.301 (9.30)
        FD = (1.0 / r / c) * pow(p_ratio, -g2); // dF/dp
    }
    else {
        //shock wave
        qrt = sqrt((g5 / r) / (g6 * p + p_curr));
        F = (p_curr - p) * qrt; //Toro p.301 (9.30)
        FD = (1.0 - 0.5*(p_curr - p)/(g6 * p + p_curr))*qrt; // dF/dp
    }
}

void find_in_which_part(InitialState& IS, double pl, double vl, double rl, double cl, double pr, double vr, double rr, double cr, double p_cont, double v_cont, double s, double& p_res, double& v_res, double& r_res) {
    double g1 = (IS.gamma - 1.0) / 2. / IS.gamma;
    double g2 = (IS.gamma + 1.0) / 2. / IS.gamma;
    double g3 = 2.0 * IS.gamma / (IS.gamma - 1.0);
    double g4 = 2.0 / (IS.gamma - 1.0);
    double g5 = 2.0 / (IS.gamma + 1.0);
    double g6 = (IS.gamma - 1.0) / (IS.gamma + 1.0);
    double g7 = (IS.gamma - 1.0) / 2.;
    double g8 = IS.gamma - 1.0;

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
                //           
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                cml = cl * pow(p_cont / pl, g1);
                stl = v_cont - cml;
                if (s > stl) {
                    r = rl * pow(p_cont / pl, 1.0 / IS.gamma); //Toro p.134 (4.53)
                    v = v_cont;
                    p = p_cont;
                }
                else {
                    //inside left rarefaction
                    //Toro p.135 (4.56)
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
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                r = rl * (p_ratio + g6) / (p_ratio * g6 + 1.0); //Toro p.133 (4.50)
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
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                r = rr * (p_ratio + g6) / (p_ratio * g6 + 1.0); //Toro p.135 (4.57)
                v = v_cont;
                p = p_cont;
            }
        }
        else {
            //right rarefaction
            shr = vr + cr;
            if (s >= shr) {
                //            
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                cmr = cr * pow(p_cont / pr, g1);
                str = v_cont + cmr;
                if (s <= str) {
                    r = rr * pow(p_cont / pr, 1.0 / IS.gamma); //Toro p.136 (4.60)
                    v = v_cont;
                    p = p_cont;
                }
                else {
                    //inside right rarefaction
                    //Toro p.136 (4.63)
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

//Toro p.155. Subroutine STARPU
void Newton_contact(InitialState& IS, double pl, double vl, double rl, double cl, double pr, double vr, double rr, double cr, double& p_cont, double& v_cont) {
    double p_old; //previous iteration
    double fl, fr;
    double fl_diff, fr_diff; //derivative of f

    int step = 0;
    int max_iter = 500;

    double change; //convergence parameter
    double eps = 1.e-6; //accuracy

    if (2. * (cl + cr) / (IS.gamma - 1.) <= (vr - vl)) { //Toro p.127 (4.40)
        std::cerr << "Vacuum was generated" << std::endl;
    }

    //initial pressure for iterational algorithm
    p_old = Newton_find_P(IS, pl, vl, rl, cl, pr, vr, rr, cr);
    if (p_old < 0.) {
        p_old = eps;
        std::cerr << "Initial pressure < 0" << std::endl;
    }

    //itterational algorithm
    do {
        calc_flux(IS, fl, fl_diff, p_old, pl, vl, rl, cl);
        calc_flux(IS, fr, fr_diff, p_old, pr, vr, rr, cr);
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
    v_cont = 0.5 * (vl + vr + fr - fl); //Toro p.125
}


void Riemann_solver(InitialState& IS, double rl, double vl, double pl, double rr, double vr, double pr, double& p_res, double& v_res, double& r_res) {
    double p_cont, v_cont;
    double p, v, r;
    double cr, cl;

    cl = calc_c(IS.gamma, pl, rl);
    cr = calc_c(IS.gamma, pr, rr);

    if (2. * (cl + cr) / (IS.gamma - 1.) <= (vr - vl)) {
        std::cerr << "Vacuum was generated" << std::endl;
    }

    Newton_contact(IS, pl, vl, rl, cl, pr, vr, rr, cr, p_cont, v_cont); // contact pressure and velocity

    find_in_which_part(IS, pl, vl, rl, cl, pr, vr, rr, cr, p_cont, v_cont, 0.0, p, v, r); //solution selection

    p_res = p;
    v_res = v;
    r_res = r;
}


void Godunov_flux_x(InitialState& IS, double p, double vx, double vy, double r, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
    double m, impx, impy, e;
    noncons_to_cons(IS, p, vx, vy, r, m, impx, impy, e);

    Fm =  r * vx;
    Fimpx = Fm * vx + p;
    Fimpy = Fm * vy;
    Fe = (p + e) * vx;
}


void Godunov_flux_y(InitialState& IS, double p, double vx, double vy, double r, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
    double m, impx, impy, e;
    noncons_to_cons(IS, p, vx, vy, r, m, impx, impy, e);

    Fm =  r * vy;
    Fimpx = Fm * vx;
    Fimpy = Fm * vy + p;
    Fe = (p + e) * vy;
}


void cons_flux_x(InitialState& IS, double m, double impx, double impy, double e, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
    double p, vx, vy, r;
    cons_to_noncons(IS, p, vx, vy, r, m, impx, impy, e);

    Fm = r * vx;
    Fimpx = Fm * vx + p;
    Fimpy = Fm * vy;
    Fe = (p + e) * vx;
}


void cons_flux_y(InitialState& IS, double m, double impx, double impy, double e, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
    double p, vx, vy, r;
    cons_to_noncons(IS, p, vx, vy, r, m, impx, impy, e);

    Fm = r * vy;
    Fimpx = Fm * vx;
    Fimpy = Fm * vy + p;
    Fe = (p + e) * vy;
}


void Godunov_solver_x(InitialState& IS, double ml, double impxl, double impyl, double el, double mr, double impxr, double impyr, double er, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
    double p, vx, vy, r;
    double pl, vxl, vyl, rl;
    double pr, vxr, vyr, rr;

    cons_to_noncons(IS, pl, vxl, vyl, rl, ml, impxl, impyl, el);
    cons_to_noncons(IS, pr, vxr, vyr, rr, mr, impxr, impyr, er);

    Riemann_solver(IS, rl, vxl, pl, rr, vxr, pr, p, vx, r);

    if (vx >= 0)
        vy = impyl / ml;
    else
        vy = impyr / mr;
    
    Godunov_flux_x(IS, p, vx, vy, r, Fm, Fimpx, Fimpy, Fe);
}


void Godunov_solver_y(InitialState& IS, double md, double impxd, double impyd, double ed, double mu, double impxu, double impyu, double eu, double& Fm, double& Fimpx, double& Fimpy, double& Fe) {
    double p, vx, vy, r;
    double pd, vxd, vyd, rd;
    double pu, vxu, vyu, ru;

    cons_to_noncons(IS, pd, vxd, vyd, rd, md, impxd, impyd, ed);
    cons_to_noncons(IS, pu, vxu, vyu, ru, mu, impxu, impyu, eu);

    Riemann_solver(IS, rd, vyd, pd, ru, vyu, pu, p, vy, r);

    if (vy >= 0)
        vx = impxd / md;
    else
        vx = impxu / mu;
    
    Godunov_flux_y(IS, p, vx, vy, r, Fm, Fimpx, Fimpy, Fe);
}


double minmod(double a, double b, int mod_type){ // 1, 2 - Kolgan, 3 - Osher, 4 - VanLeer
    if (mod_type == 1) {
        if (abs(a) <= abs(b)) {
            return a;
        }
        else {
            return b;
        }
    }
    else if (mod_type == 2) {
        double c = (a + b) / 2.;
        if ((abs(a) <= abs(b)) & (abs(a) <= abs(c))) {
            return a;
        }
        else if ((abs(b) <= abs(a)) & (abs(b) <= abs(c))) {
            return b;
        }
        else {
            return c;
        }
    }
    else if (mod_type == 3) {
        if (a * a <= a * b) {
            return a;
        }
        else if (b * b < a * b) {
            return b;
        }
        else {
            return 0.;
        }
    }
    else if (mod_type == 4) {
        return minmod((a + b) / 2., 2. * minmod(a, b, 3), 3);
    }
	return minmod(a, b, 3);
}


//Godunov-Kolgan reconstruction on the left cell boundary
void Kolgan_left(double qi, double qi_1, double qi_2, double qi1, double xi, double xi_1, double xi_2, double xi1, double& ql, double& qr, int mod_type){
    double dif_l, dif_r, dif_c; //linear coeffs
    
    dif_l = (qi - qi_1)/(xi - xi_1);
    dif_r = (qi1 - qi)/(xi1 - xi);
    dif_c = minmod(dif_l, dif_r, mod_type);

    qr = -dif_c * (xi - xi_1) / 2. + qi;

    dif_r = dif_l;
    dif_l = (qi_1 - qi_2)/(xi_1 - xi_2);
    dif_c = minmod(dif_l, dif_r, mod_type);

    ql = dif_c * (xi - xi_1) / 2. + qi_1;    
}

//Godunov-Kolgan reconstruction on the right cell boundary
void Kolgan_right(double qi, double qi_1, double qi1, double qi2, double xi, double xi_1, double xi1, double xi2, double& ql, double& qr, int mod_type){
    double dif_l, dif_r, dif_c; //linear coeffs
    
    dif_l = (qi - qi_1)/(xi - xi_1);
    dif_r = (qi1 - qi)/(xi1 - xi);
    dif_c = minmod(dif_l, dif_r, mod_type);

    ql = dif_c * (xi1 - xi) / 2. + qi;

    dif_l = dif_r;
    dif_r = (qi2 - qi1)/(xi2 - xi1);
    dif_c = minmod(dif_l, dif_r, mod_type);

    qr = - dif_c * (xi1 - xi) / 2. + qi1;    
}

bool isPrime(int p) {
    if (p <= 3) {
        return false;
    }
    if (p % 2 == 0) {
        return false;
    }
    for (int i = 4; i <= std::sqrt(p); i += 2) {
        if (p % i == 0) {
            return false;
        }
    }
    return true;
}

void calc_Fx(InitialState& IS, vec3d& u, vec3d& F) {
    double r, vx, vy, p, e;
    for (int i = 0; i < IS.Nx + 2 * IS.fict; i++) {
        for(int j = 0; j < IS.Ny + 2 * IS.fict; j++){
            r = u[0][i][j];
            vx = u[1][i][j] / r;
            vy = u[2][i][j] / r;
            e = u[3][i][j];
            p = (IS.gamma - 1.0) * (e - 0.5 * r * (vx * vx + vy * vy));
            F[0][i][j] = vx * r;
            F[1][i][j] = vx * vx * r + p;
            F[2][i][j] = vx * vy * r;
            F[3][i][j] = vx * (e + p);
        }
    }
}

void calc_Fy(InitialState& IS, vec3d& u, vec3d& F) {
    double r, vx, vy, p, e;
    for (int i = 0; i < IS.Nx + 2 * IS.fict; i++) {
        for(int j = 0; j < IS.Ny + 2 * IS.fict; j++){
            r = u[0][i][j];
            vx = u[1][i][j] / r;
            vy = u[2][i][j] / r;
            e = u[3][i][j];
            p = (IS.gamma - 1.0) * (e - 0.5 * r * (vx * vx + vy * vy));
            F[0][i][j] = vy * r;
            F[1][i][j] = vy * vx * r;
            F[2][i][j] = vy * vy * r + p;
            F[3][i][j] = vy * (e + p);
        }
    }
}

double calc_A(InitialState& IS, vec3d& u) {
    double new_step = 0.;
    double p, vx, vy, r, e, c;
    double a_step;
    for (int i = IS.fict; i < IS.Nx + IS.fict; ++i) {
        for(int j = IS.fict; j < IS.Ny + IS.fict; ++j){
            r = u[0][i][j];
            vx = u[1][i][j] / r;
            vy = u[2][i][j] / r;
            e = u[3][i][j];
            p = (IS.gamma - 1.0) * (e - 0.5 * r * (vx * vx + vy * vy));
            c = calc_c(IS.gamma, p, r);
            a_step = c;
            if (a_step > new_step) {
                new_step = a_step;
            }
        }
    }
    return new_step;
}

double WENO_reconstruction_minus(InitialState& IS, double um2, double um1, double u, double up1, double up2) {
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
        if (b[j] < eps){
            a[j] = gammaRight[j] / pow(eps, 2.0);
        } else{
            a[j] = gammaRight[j] / pow(b[j], 2.0);
        }
        
        a_sum += a[j];
    }

    for (int j = 0; j < 3; j++) {
        w[j] = a[j] / a_sum;
    }
    return w[0] * P[0] + w[1] * P[1] + w[2] * P[2];
}

double WENO_reconstruction_plus(InitialState& IS, double um2, double um1, double u, double up1, double up2) {
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
        if (b[j] < eps){
            a[j] = gammaLeft[j] / pow(eps, 2.0);
        } else{
            a[j] = gammaLeft[j] / pow(b[j], 2.0);
        }
        a_sum += a[j];
    }

    for (int j = 0; j < 3; j++) {
        w[j] = a[j] / a_sum;
    }
    return w[0] * P[0] + w[1] * P[1] + w[2] * P[2];
}

void WENO5(InitialState& IS, vec3d& u, vec3d& L, double dx, double dy) {
    vec3d Fx(4, vec2d(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict)));
    vec3d Fy(4, vec2d(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict)));
    vec3d Fx_plus(4, vec2d(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict)));
    vec3d Fy_plus(4, vec2d(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict)));
    vec3d Fx_minus(4, vec2d(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict)));
    vec3d Fy_minus(4, vec2d(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict)));
    

    calc_Fx(IS, u, Fx);
    calc_Fy(IS, u, Fy);
    double A = calc_A(IS, u);

    for (int k = 0; k < 4; k++) {
        for (int i = 0; i < IS.Nx + 2 * IS.fict; i++) {
            for(int j = 0; j < IS.Ny + 2 * IS.fict; j++){
                Fx_plus[k][i][j] = 0.5 * (Fx[k][i][j] + A * u[k][i][j]);
                Fx_minus[k][i][j] = 0.5 * (Fx[k][i][j] - A * u[k][i][j]);
                Fy_plus[k][i][j] = 0.5 * (Fy[k][i][j] + A * u[k][i][j]);
                Fy_minus[k][i][j] = 0.5 * (Fy[k][i][j] - A * u[k][i][j]);
            }
        }
    }

    vec3d FlowL(4, vec2d(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict)));
    vec3d FlowR(4, vec2d(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict)));
    for (int k = 0; k < 4; k++) {
        for (int i = IS.fict - 1; i < IS.Nx + IS.fict + 1; i++) {
            for(int j = IS.fict - 1; j < IS.Ny + IS.fict + 1; j++){
                FlowL[k][i][j] = WENO_reconstruction_minus(IS, Fx_minus[k][i - 2][j], Fx_minus[k][i - 1][j], Fx_minus[k][i][j], Fx_minus[k][i + 1][j], Fx_minus[k][i + 2][j]);
                FlowR[k][i][j] = WENO_reconstruction_plus(IS, Fx_plus[k][i - 2][j], Fx_plus[k][i - 1][j], Fx_plus[k][i][j], Fx_plus[k][i + 1][j], Fx_plus[k][i + 2][j]);
            }
        }
    }
    // Boundary_x(IS, FlowL[0], FlowL[1], FlowL[2], FlowL[3]);
    // Boundary_x(IS, FlowR[0], FlowR[1], FlowR[2], FlowR[3]);

    vec3d FlowD(4, vec2d(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict)));
    vec3d FlowU(4, vec2d(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict)));
    for (int k = 0; k < 4; k++) {
        for (int i = IS.fict - 1; i < IS.Nx + IS.fict + 1; i++) {
            for(int j = IS.fict - 1; j < IS.Ny + IS.fict + 1; j++){
                FlowD[k][i][j] = WENO_reconstruction_minus(IS, Fy_minus[k][i][j - 2], Fy_minus[k][i][j - 1], Fy_minus[k][i][j], Fy_minus[k][i][j + 1], Fy_minus[k][i][j + 2]);
                FlowU[k][i][j] = WENO_reconstruction_plus(IS, Fy_plus[k][i][j - 2], Fy_plus[k][i][j - 1], Fy_plus[k][i][j], Fy_plus[k][i][j + 1], Fy_plus[k][i][j + 2]);
            }
        }
    }
    // Boundary_y(IS, FlowD[0], FlowD[1], FlowD[2], FlowD[3]);
    // Boundary_y(IS, FlowU[0], FlowU[1], FlowU[2], FlowU[3]);

    for (int k = 0; k < 4; k++) {
        for (int i = IS.fict; i < IS.Nx + IS.fict; i++) {
            for (int j = IS.fict; j < IS.Ny + IS.fict; j++) {
                L[k][i][j] = ((FlowR[k][i][j] + FlowL[k][i + 1][j]) - (FlowR[k][i - 1][j] + FlowL[k][i][j])) / dx + \
                 ((FlowU[k][i][j] + FlowD[k][i][j + 1]) - (FlowU[k][i][j - 1] + FlowD[k][i][j])) / dy;
            }
        }
    }
}

void split_proc(int& p, int& px, int& py){
    //if(isPrime(p)) p--;
    py = std::sqrt(p);
    while (py > 1 && p % py != 0) {
        py--;
    }

    px = p / py;
}

void split_plane(int px, int py, int myrank, int Nx, int Ny, int& Nx_proc, int& Ny_proc, int& Nx_start, int& Nx_end, int& Ny_start, int& Ny_end){
    int base_x = Nx / px;
    int base_y = Ny / py;
    int remainder_x = Nx % px;
    int remainder_y = Ny % py;

    Nx_proc = base_x + (myrank % px < remainder_x ? 1 : 0);
    Nx_start = myrank % px * base_x + std::min(myrank % px, remainder_x);
    Nx_end = Nx_start + Nx_proc;

    Ny_proc = base_y + (myrank / px < remainder_y ? 1 : 0);
    Ny_start = myrank / px * base_y + std::min(myrank / px, remainder_y);
	Ny_end = Ny_start + Ny_proc;
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

void exchange_data(InitialState& IS, int myrank, int px, int py, vec2d& m, vec2d& impx, vec2d& impy, vec2d& e, vec2d& buff_x, vec2d& buff_y, MPI_Comm comm, MPI_Status& Status){
    // X sending
    if((myrank % px) % 2 != 0){
        if((myrank % px) > 0){
            for(int i = 0; i < IS.fict; i++){
                for(int j = 0; j < IS.Ny + 2 * IS.fict; j++){
                    buff_x[i][j] = m[i + IS.fict][j];
                    buff_x[i + IS.fict][j] = impx[i + IS.fict][j];
                    buff_x[i + 2 * IS.fict][j] = impy[i + IS.fict][j];
                    buff_x[i + 3 * IS.fict][j] = e[i + IS.fict][j];
                }
            }
            send_vec2(buff_x, 4 * IS.fict, (IS.Ny + 2 * IS.fict), myrank - 1, 0, comm);
        }
        if((myrank % px) > 0){
            buff_x = recv_vec2(4 * IS.fict, IS.Ny + 2 * IS.fict, myrank - 1, 0, comm, Status);
            for(int i = 0; i < IS.fict; i ++){
                for(int j = 0; j < IS.Ny + 2 * IS.fict; j++){
                    m[i][j] = buff_x[i][j];
                    impx[i][j] = buff_x[i + IS.fict][j];
                    impy[i][j] = buff_x[i + 2 * IS.fict][j];
                    e[i][j] = buff_x[i + 3 * IS.fict][j];
                }
            }
        }
        if((myrank % px) < px - 1){
            buff_x = recv_vec2(4 * IS.fict, IS.Ny + 2 * IS.fict, myrank + 1, 0, comm, Status);
            for(int i = 0; i < IS.fict; i ++){
                for(int j = 0; j < IS.Ny + 2 * IS.fict; j++){
                    m[IS.Nx + IS.fict + i][j] = buff_x[i][j];
                    impx[IS.Nx + IS.fict + i][j] = buff_x[i + IS.fict][j];
                    impy[IS.Nx + IS.fict + i][j] = buff_x[i + 2 * IS.fict][j];
                    e[IS.Nx + IS.fict + i][j] = buff_x[i + 3 * IS.fict][j];
                }
            }
        }
        if((myrank % px) < px - 1){
            for(int i = 0; i < IS.fict; i ++){
                for(int j = 0; j < IS.Ny + 2 * IS.fict; j++){
                    buff_x[i][j] = m[IS.Nx + i][j];
                    buff_x[i + IS.fict][j] = impx[IS.Nx + i][j];
                    buff_x[i + 2 * IS.fict][j] = impy[IS.Nx + i][j];
                    buff_x[i + 3 * IS.fict][j] = e[IS.Nx + i][j];
                }
            }
            send_vec2(buff_x, 4 * IS.fict, (IS.Ny + 2 * IS.fict), myrank + 1, 0, comm);    
        }
    }
    if((myrank % px) % 2 == 0){
        if((myrank % px) < px - 1){
            buff_x = recv_vec2(4 * IS.fict, IS.Ny + 2 * IS.fict, myrank + 1, 0, comm, Status);
            for(int i = 0; i < IS.fict; i ++){
                for(int j = 0; j < IS.Ny + 2 * IS.fict; j++){
                    m[IS.Nx + IS.fict + i][j] = buff_x[i][j];
                    impx[IS.Nx + IS.fict + i][j] = buff_x[i + IS.fict][j];
                    impy[IS.Nx + IS.fict + i][j] = buff_x[i + 2 * IS.fict][j];
                    e[IS.Nx + IS.fict + i][j] = buff_x[i + 3 * IS.fict][j];
                }
            }
        }
        if((myrank % px) < px - 1){
            for(int i = 0; i < IS.fict; i ++){
                for(int j = 0; j < IS.Ny + 2 * IS.fict; j++){
                    buff_x[i][j] = m[IS.Nx + i][j];
                    buff_x[i + IS.fict][j] = impx[IS.Nx + i][j];
                    buff_x[i + 2 * IS.fict][j] = impy[IS.Nx + i][j];
                    buff_x[i + 3 * IS.fict][j] = e[IS.Nx + i][j];
                }
            }
            send_vec2(buff_x, 4 * IS.fict, (IS.Ny + 2 * IS.fict), myrank + 1, 0, comm); 
        }
        if((myrank % px) > 0){
            for(int i = 0; i < IS.fict; i ++){
                for(int j = 0; j < IS.Ny + 2 * IS.fict; j++){
                    buff_x[i][j] = m[i + IS.fict][j];
                    buff_x[i + IS.fict][j] = impx[i + IS.fict][j];
                    buff_x[i + 2 * IS.fict][j] = impy[i + IS.fict][j];
                    buff_x[i + 3 * IS.fict][j] = e[i + IS.fict][j];
                }
            }
            send_vec2(buff_x, 4 * IS.fict, (IS.Ny + 2 * IS.fict), myrank - 1, 0, comm);
        }
        if((myrank % px) > 0){
            buff_x = recv_vec2(4 * IS.fict, IS.Ny + 2 * IS.fict, myrank - 1, 0, comm, Status);
            for(int i = 0; i < IS.fict; i ++){
                for(int j = 0; j < IS.Ny + 2 * IS.fict; j++){
                    m[i][j] = buff_x[i][j];
                    impx[i][j] = buff_x[i + IS.fict][j];
                    impy[i][j] = buff_x[i + 2 * IS.fict][j];
                    e[i][j] = buff_x[i + 3 * IS.fict][j];
                }
            }
        }
    }

    // Y sending
    if((myrank / px) % 2 != 0){
        if((myrank / px) > 0){
            for(int i = 0; i < IS.fict; i++){
                for(int j = 0; j < IS.Nx + 2 * IS.fict; j++){
                    buff_y[i][j] = m[j][i + IS.fict];
                    buff_y[i + IS.fict][j] = impx[j][i + IS.fict];
                    buff_y[i + 2 * IS.fict][j] = impy[j][i + IS.fict];
                    buff_y[i + 3 * IS.fict][j] = e[j][i + IS.fict];
                }
            }
            send_vec2(buff_y, 4 * IS.fict, IS.Nx + 2 * IS.fict, myrank - px, 1, comm);
        }
        if((myrank / px) > 0){
            buff_y = recv_vec2(4 * IS.fict, IS.Nx + 2 * IS.fict, myrank - px, 1, comm, Status);
            for(int i = 0; i < IS.fict; i++){
                for(int j = 0; j < IS.Nx + 2 * IS.fict; j++){
                    m[j][i] = buff_y[i][j];
                    impx[j][i] = buff_y[i + IS.fict][j];
                    impy[j][i] = buff_y[i + 2 * IS.fict][j];
                    e[j][i] = buff_y[i + 3 * IS.fict][j];
                }
            }
        }
        if((myrank / px) < py - 1){
            buff_y = recv_vec2(4 * IS.fict, IS.Nx + 2 * IS.fict, myrank + px, 1, comm, Status);
            for(int i = 0; i < IS.fict; i++){
                for(int j = 0; j < IS.Nx + 2 * IS.fict; j++){
                    m[j][IS.Ny + IS.fict + i] = buff_y[i][j];
                    impx[j][IS.Ny + IS.fict + i] = buff_y[i + IS.fict][j];
                    impy[j][IS.Ny + IS.fict + i] = buff_y[i + 2 * IS.fict][j];
                    e[j][IS.Ny + IS.fict + i] = buff_y[i + 3 * IS.fict][j];
                }
            }
        }
        if((myrank / px) < py - 1){
            for(int i = 0; i < IS.fict; i++){
                for(int j = 0; j < IS.Nx + 2 * IS.fict; j++){
                    buff_y[i][j] = m[j][IS.Ny + i];
                    buff_y[i + IS.fict][j] = impx[j][IS.Ny + i];
                    buff_y[i + 2 * IS.fict][j] = impy[j][IS.Ny + i];
                    buff_y[i + 3 * IS.fict][j] = e[j][IS.Ny + i];
                }
            }
            send_vec2(buff_y, 4 * IS.fict, IS.Nx + 2 * IS.fict, myrank + px, 1, comm);    
        }
    }
    if((myrank / px) % 2 == 0){
        if((myrank / px) < py - 1){
            buff_y = recv_vec2(4 * IS.fict, IS.Nx + 2 * IS.fict, myrank + px, 1, comm, Status);
            for(int i = 0; i < IS.fict; i++){
                for(int j = 0; j < IS.Nx + 2 * IS.fict; j++){
                    m[j][IS.Ny + IS.fict + i] = buff_y[i][j];
                    impx[j][IS.Ny + IS.fict + i] = buff_y[i + IS.fict][j];
                    impy[j][IS.Ny + IS.fict + i] = buff_y[i + 2 * IS.fict][j];
                    e[j][IS.Ny + IS.fict + i] = buff_y[i + 3 * IS.fict][j];
                }
            }
        }
        if((myrank / px) < py - 1){
            for(int i = 0; i < IS.fict; i++){
                for(int j = 0; j < IS.Nx + 2 * IS.fict; j++){
                    buff_y[i][j] = m[j][IS.Ny + i];
                    buff_y[i + IS.fict][j] = impx[j][IS.Ny + i];
                    buff_y[i + 2 * IS.fict][j] = impy[j][IS.Ny + i];
                    buff_y[i + 3 * IS.fict][j] = e[j][IS.Ny + i];
                }
            }
            send_vec2(buff_y, 4 * IS.fict, IS.Nx + 2 * IS.fict, myrank + px, 1, comm); 
        }
        if((myrank / px) > 0){
            for(int i = 0; i < IS.fict; i++){
                for(int j = 0; j < IS.Nx + 2 * IS.fict; j++){
                    buff_y[i][j] = m[j][i + IS.fict];
                    buff_y[i + IS.fict][j] = impx[j][i + IS.fict];
                    buff_y[i + 2 * IS.fict][j] = impy[j][i + IS.fict];
                    buff_y[i + 3 * IS.fict][j] = e[j][i + IS.fict];
                }
            }
            send_vec2(buff_y, 4 * IS.fict, IS.Nx + 2 * IS.fict, myrank - px, 1, comm);
        }
        if((myrank / px) > 0){
            buff_y = recv_vec2(4 * IS.fict, IS.Nx + 2 * IS.fict, myrank - px, 1, comm, Status);
            for(int i = 0; i < IS.fict; i++){
                for(int j = 0; j < IS.Nx + 2 * IS.fict; j++){
                    m[j][i] = buff_y[i][j];
                    impx[j][i] = buff_y[i + IS.fict][j];
                    impy[j][i] = buff_y[i + 2 * IS.fict][j];
                    e[j][i] = buff_y[i + 3 * IS.fict][j];
                }
            }
        }
    }
}

void Multiproc_solve(InitialState& IS, std::string out_dir, int& myrank, int& size, MPI_Comm comm){
    MPI_Status Status;
    int px, py;
    split_proc(size, px, py);
    size = px * py;
    
    int Nx_proc, Ny_proc, Nx_start, Nx_end, Ny_start, Ny_end;
    split_plane(px, py, myrank, IS.Nx, IS.Ny, Nx_proc, Ny_proc, Nx_start, Nx_end, Ny_start, Ny_end);

    double dx = (IS.x_end - IS.x_start) / IS.Nx;
    double dy = (IS.y_end - IS.y_start) / IS.Ny;

    double xp_start, xp_end, yp_start, yp_end;
    xp_start = Nx_start * dx;
    xp_end = Nx_end * dx;
    yp_start = Ny_start * dy;
    yp_end = Ny_end * dy;
    InitialState params_proc = change_params(IS, Nx_proc, Ny_proc, xp_start, xp_end, yp_start, yp_end, myrank, px, py);
    //printf("myrank = %d, px = %d, py = %d, Nx_proc = %d, Ny_proc = %d, x = (%.5f, %.5f), y = (%.5f, %.5f)\n", myrank, px, py, Nx_proc, Ny_proc, xp_start, xp_end, yp_start, yp_end);
    //printf("myrank = %d, b_right = %d, b_left = %d, b_down = %d, b_up = %d\n", myrank, params_proc.b_right, params_proc.b_left, params_proc.b_down, params_proc.b_up);
    
    vec1d x(Nx_proc + 1 + 2 * IS.fict); 
    vec1d xc(Nx_proc + 1 + 2 * IS.fict);
    vec1d y(Ny_proc + 1 + 2 * IS.fict); 
    vec1d yc(Ny_proc + 1 + 2 * IS.fict);

    vec2d p(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d vx(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d vy(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d r(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));   
    vec2d m(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d impx(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d impy(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d e(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d ei(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));

    vec2d m_new(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d impx_new(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d impy_new(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d e_new(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));

    // for Kolgan/Rodionov
    vec2d m_kolganl(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d impx_kolganl(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d impy_kolganl(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d e_kolganl(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d m_kolganr(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d impx_kolganr(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d impy_kolganr(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d e_kolganr(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));

    vec2d m_kolgand(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d impx_kolgand(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d impy_kolgand(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d e_kolgand(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d m_kolganu(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d impx_kolganu(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d impy_kolganu(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d e_kolganu(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));

    vec2d m05(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d impx05(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d impy05(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d e05(Nx_proc + 2*IS.fict, vec1d(Ny_proc + 2*IS.fict));

    double dm_x, dimpx_x, dimpy_x, de_x;
    double dm_y, dimpx_y, dimpy_y, de_y;

    // for WENO
    std::vector<vec3d> U(4, vec3d(4, vec2d(params_proc.Nx + 2*params_proc.fict, vec1d(params_proc.Ny + 2*params_proc.fict))));
    std::vector<vec1d> rk = { {1., 0., 0., 0.}, {1., 0., 0., 1.}, {0.75, 0.25, 0., 0.25}, {1. / 3., 0., 2. / 3., 2. / 3.} };
    vec3d L(4, vec2d(params_proc.Nx + 2*params_proc.fict, vec1d(params_proc.Ny + 2*params_proc.fict)));

    vec2d buff_x(4 * IS.fict, vec1d(Ny_proc + 2*IS.fict));
    vec2d buff_y(4 * IS.fict, vec1d(Nx_proc + 2*IS.fict));

    double Fml, Fimpxl, Fimpyl, Fel, Fmr, Fimpxr, Fimpyr, Fer;  //cell x fluxes
    double Fmd, Fimpxd, Fimpyd, Fed, Fmu, Fimpxu, Fimpyu, Feu;  //cell y fluxes

    double ml, impxl, impyl, el, mr, impxr, impyr, er; // left and right cell values of cons
    double md, impxd, impyd, ed, mu, impxu, impyu, eu; // down and up cell values of cons

    int step = 0;
    double time = 0.;
    double dt;

    x_grid(params_proc, xc, x);
    y_grid(params_proc, yc, y);
    Init(params_proc, xc, yc, p, vx, vy, r, m, impx, impy, e);
    Boundary_x(params_proc, m, impx, impy, e);
    Boundary_y(params_proc, m, impx, impy, e);
    calc_ei(params_proc, p, r, ei);
    if(myrank == 0) {
        if (fs::exists(out_dir)) {
            fs::remove_all(out_dir);
        }
    }
    write_out_p(params_proc, xc, yc, p, vx, vy, r, m, impx, impy, e, ei, step, time, out_dir, myrank, size, comm);
    
    while (time < params_proc.t_end) {
        dt = get_dt(params_proc, x, y, m, impx, impy, e);
        exchange_dt(params_proc, dt, x, y, m, impx, impy, e, myrank, size, comm); // time step
        time += dt;
        step++;
        if(IS.s_type != 4){
            if (IS.s_type == 2 || IS.s_type == 3) { // Kolgan/Rodionov 
                for (int i = params_proc.fict; i < params_proc.Nx + params_proc.fict + 1; ++i) {
                    for (int j = params_proc.fict; j < params_proc.Ny + params_proc.fict + 1; ++j) {
                        Kolgan_left(m[i][j], m[i - 1][j], m[i - 2][j], m[i + 1][j], xc[i], xc[i - 1], xc[i - 2], xc[i + 1], m_kolganl[i][j], m_kolganr[i][j], params_proc.mod_type);
                        Kolgan_left(impx[i][j], impx[i - 1][j], impx[i - 2][j], impx[i + 1][j], xc[i], xc[i - 1], xc[i - 2], xc[i + 1], impx_kolganl[i][j], impx_kolganr[i][j], params_proc.mod_type);
                        Kolgan_left(impy[i][j], impy[i - 1][j], impy[i - 2][j], impy[i + 1][j], xc[i], xc[i - 1], xc[i - 2], xc[i + 1], impy_kolganl[i][j], impy_kolganr[i][j], params_proc.mod_type);
                        Kolgan_left(e[i][j], e[i - 1][j], e[i - 2][j], e[i + 1][j], xc[i], xc[i - 1], xc[i - 2], xc[i + 1], e_kolganl[i][j], e_kolganr[i][j], params_proc.mod_type);
                        
                        if(i < params_proc.Nx + params_proc.fict){
                            Kolgan_right(m[i][j], m[i - 1][j], m[i + 1][j], m[i + 2][j], xc[i], xc[i - 1], xc[i + 1], xc[i + 2], m_kolganl[i + 1][j], m_kolganr[i + 1][j], params_proc.mod_type);
                            Kolgan_right(impx[i][j], impx[i - 1][j], impx[i + 1][j], impx[i + 2][j], xc[i], xc[i - 1], xc[i + 1], xc[i + 2], impx_kolganl[i + 1][j], impx_kolganr[i + 1][j], params_proc.mod_type);
                            Kolgan_right(impy[i][j], impy[i - 1][j], impy[i + 1][j], impy[i + 2][j], xc[i], xc[i - 1], xc[i + 1], xc[i + 2], impy_kolganl[i + 1][j], impy_kolganr[i + 1][j], params_proc.mod_type);
                            Kolgan_right(e[i][j], e[i - 1][j], e[i + 1][j], e[i + 2][j], xc[i], xc[i - 1], xc[i + 1], xc[i + 2], e_kolganl[i + 1][j], e_kolganr[i + 1][j], params_proc.mod_type);
                        }

                        Kolgan_left(m[i][j], m[i][j - 1], m[i][j - 2], m[i][j + 1], yc[j], yc[j - 1], yc[j - 2], yc[j + 1], m_kolgand[i][j], m_kolganu[i][j], params_proc.mod_type);
                        Kolgan_left(impx[i][j], impx[i][j - 1], impx[i][j - 2], impx[i][j + 1], yc[j], yc[j - 1], yc[j - 2], yc[j + 1], impx_kolgand[i][j], impx_kolganu[i][j], params_proc.mod_type);
                        Kolgan_left(impy[i][j], impy[i][j - 1], impy[i][j - 2], impy[i][j + 1], yc[j], yc[j - 1], yc[j - 2], yc[j + 1], impy_kolgand[i][j], impy_kolganu[i][j], IS.mod_type);
                        Kolgan_left(e[i][j], e[i][j - 1], e[i][j - 2], e[i][j + 1], yc[j], yc[j - 1], yc[j - 2], yc[j + 1], e_kolgand[i][j], e_kolganu[i][j], params_proc.mod_type);
                        
                        if(j < params_proc.Ny + params_proc.fict){
                            Kolgan_right(m[i][j], m[i][j - 1], m[i][j + 1], m[i][j + 2], yc[j], yc[j - 1], yc[j + 1], yc[j + 2], m_kolgand[i][j + 1], m_kolganu[i][j + 1], params_proc.mod_type);
                            Kolgan_right(impx[i][j], impx[i][j - 1], impx[i][j + 1], impx[i][j + 2], yc[j], yc[j - 1], yc[j + 1], yc[j + 2], impx_kolgand[i][j + 1], impx_kolganu[i][j + 1], params_proc.mod_type);
                            Kolgan_right(impy[i][j], impy[i][j - 1], impy[i][j + 1], impy[i][j + 2], yc[j], yc[j - 1], yc[j + 1], yc[j + 2], impy_kolgand[i][j + 1], impy_kolganu[i][j + 1], params_proc.mod_type);
                            Kolgan_right(e[i][j], e[i][j - 1], e[i][j + 1], e[i][j + 2], yc[j], yc[j - 1], yc[j + 1], yc[j + 2], e_kolgand[i][j + 1], e_kolganr[i][j + 1], params_proc.mod_type);
                        }
                    }
                }
            }
            if (IS.s_type == 3) { // Rodionov
                for (int i = params_proc.fict; i < params_proc.Nx + params_proc.fict; ++i) {
                    for (int j = params_proc.fict; j < params_proc.Ny + params_proc.fict; ++j) {
                        // Predictor
                        cons_flux_x(IS, m_kolganr[i][j], impx_kolganr[i][j], impy_kolganr[i][j], e_kolganr[i][j], Fml, Fimpxl, Fimpyl, Fel);
                        cons_flux_x(IS, m_kolganl[i + 1][j], impx_kolganl[i + 1][j], impy_kolganl[i + 1][j], e_kolganl[i + 1][j], Fmr, Fimpxr, Fimpyr, Fer);
                        cons_flux_y(IS, m_kolganu[i][j], impx_kolganu[i][j], impy_kolganu[i][j], e_kolganu[i][j], Fmd, Fimpxd, Fimpyd, Fed);
                        cons_flux_y(IS, m_kolgand[i][j + 1], impx_kolgand[i][j + 1], impy_kolgand[i][j + 1], e_kolgand[i][j + 1], Fmu, Fimpxu, Fimpyu, Feu);
                        m05[i][j] = m[i][j] - dt * ((Fmr - Fml) / (xc[i] - xc[i - 1]) + (Fmu - Fmd) / (yc[j] - yc[j - 1]));
                        impx05[i][j] = impx[i][j] - dt * ((Fimpxr - Fimpxl) / (xc[i] - xc[i - 1]) + (Fimpxu - Fimpxd) / (yc[j] - yc[j - 1]));
                        impy05[i][j] = impy[i][j] - dt * ((Fimpyr - Fimpyl) / (xc[i] - xc[i - 1]) + (Fimpyu - Fimpyd) / (yc[j] - yc[j - 1]));
                        e05[i][j] = e[i][j] - dt * ((Fer - Fel) / (xc[i] - xc[i - 1]) + (Feu - Fed) / (yc[j] - yc[j - 1]));
                    }
				}
				for (int i = params_proc.fict; i < params_proc.Nx + params_proc.fict + 1; ++i) {
                    for (int j = params_proc.fict; j < params_proc.Ny + params_proc.fict; ++j) {
                        // Corrector x
                        dm_x = (m_kolganl[i + 1][j] - m_kolganr[i][j]);
                        dimpx_x = (impx_kolganl[i + 1][j] - impx_kolganr[i][j]);
                        dimpy_x = (impy_kolganl[i + 1][j] - impy_kolganr[i][j]);
                        de_x = (e_kolganl[i + 1][j] - e_kolganr[i][j]);

                        m_kolganr[i][j] = 0.5 * (m[i][j] + m05[i][j]) - 0.5 * dm_x;
                        impx_kolganr[i][j] = 0.5 * (impx[i][j] + impx05[i][j]) - 0.5 * dimpx_x;
                        impy_kolganr[i][j] = 0.5 * (impy[i][j] + impy05[i][j]) - 0.5 * dimpy_x;
                        e_kolganr[i][j] = 0.5 * (e[i][j] + e05[i][j]) - 0.5 * de_x;
                        m_kolganl[i + 1][j] = 0.5 * (m[i][j] + m05[i][j]) + 0.5 * dm_x;
                        impx_kolganl[i + 1][j] = 0.5 * (impx[i][j] + impx05[i][j]) + 0.5 * dimpx_x;
                        impy_kolganl[i + 1][j] = 0.5 * (impy[i][j] + impy05[i][j]) + 0.5 * dimpy_x;
                        e_kolganl[i + 1][j] = 0.5 * (e[i][j] + e05[i][j]) + 0.5 * de_x;

                        // Corrector y
                        dm_y = (m_kolgand[i][j + 1] - m_kolganu[i][j]);
                        dimpx_y = (impx_kolgand[i][j + 1] - impx_kolganu[i][j]);
                        dimpy_y = (impy_kolgand[i][j + 1] - impy_kolganu[i][j]);
                        de_y = (e_kolgand[i][j + 1] - e_kolganu[i][j]);

                        m_kolganu[i][j] = 0.5 * (m[i][j] + m05[i][j]) - 0.5 * dm_y;
                        impx_kolganu[i][j] = 0.5 * (impx[i][j] + impx05[i][j]) - 0.5 * dimpx_y;
                        impy_kolganu[i][j] = 0.5 * (impy[i][j] + impy05[i][j]) - 0.5 * dimpy_y;
                        e_kolganu[i][j] = 0.5 * (e[i][j] + e05[i][j]) - 0.5 * de_y;
                        m_kolgand[i][j + 1] = 0.5 * (m[i][j] + m05[i][j]) + 0.5 * dm_y;
                        impx_kolgand[i][j + 1] = 0.5 * (impx[i][j] + impx05[i][j]) + 0.5 * dimpx_y;
                        impy_kolgand[i][j + 1] = 0.5 * (impy[i][j] + impy05[i][j]) + 0.5 * dimpy_y;
                        e_kolgand[i][j + 1] = 0.5 * (e[i][j] + e05[i][j]) + 0.5 * de_y;
                    }
                }
            }
            for (int i = params_proc.fict; i < params_proc.Nx + params_proc.fict; ++i) {
                for (int j = params_proc.fict; j < params_proc.Ny + params_proc.fict; ++j) {
                    //left cell bounday flux
                    if (params_proc.s_type == 1) {  // Godunov X
                        ml = m[i - 1][j];
                        impxl = impx[i - 1][j];
                        impyl = impy[i - 1][j];
                        el = e[i - 1][j];
                        mr = m[i][j];
                        impxr = impx[i][j];
                        impyr = impy[i][j];
                        er = e[i][j];
                    }
                    else if (params_proc.s_type == 2 || params_proc.s_type == 3){ // Kolgan/Rodionov X
                        ml = m_kolganl[i][j];
                        mr = m_kolganr[i][j];
                        impxl = impx_kolganl[i][j];
                        impxr = impx_kolganr[i][j];
                        impyl = impy_kolganl[i][j];
                        impyr = impy_kolganr[i][j];
                        el = e_kolganl[i][j];
                        er = e_kolganr[i][j];
                    }
                    Godunov_solver_x(params_proc, ml, impxl, impyl, el, mr, impxr, impyr, er, Fml, Fimpxl, Fimpyl, Fel);
                    
                    //right cell boundary flux
                    if (params_proc.s_type == 1) { //Godunov X
                        ml = m[i][j];
                        impxl = impx[i][j];
                        impyl = impy[i][j];
                        el = e[i][j];
                        mr = m[i + 1][j];
                        impxr = impx[i + 1][j];
                        impyr = impy[i + 1][j];
                        er = e[i + 1][j];
                    }
                    else if (params_proc.s_type == 2 || params_proc.s_type == 3){ // Kolgan/Rodionov X
                        ml = m_kolganl[i + 1][j];
                        mr = m_kolganr[i + 1][j];
                        impxl = impx_kolganl[i + 1][j];
                        impxr = impx_kolganr[i + 1][j];
                        impyl = impy_kolganl[i + 1][j];
                        impyr = impy_kolganr[i + 1][j];
                        el = e_kolganl[i + 1][j];
                        er = e_kolganr[i + 1][j];
                    }
                    Godunov_solver_x(params_proc, ml, impxl, impyl, el, mr, impxr, impyr, er, Fmr, Fimpxr, Fimpyr, Fer);

                    //down cell bounday flux
                    if (params_proc.s_type == 1) {  //Godunov Y
                        md = m[i][j - 1];
                        impxd = impx[i][j - 1];
                        impyd = impy[i][j - 1];
                        ed = e[i][j - 1];
                        mu = m[i][j];
                        impxu = impx[i][j];
                        impyu = impy[i][j];
                        eu = e[i][j];
                    }
                    else if (params_proc.s_type == 2 || params_proc.s_type == 3){ // Kolgan/Rodionov Y
                        md = m_kolgand[i][j];
                        mu = m_kolganu[i][j];
                        impxd = impx_kolgand[i][j];
                        impxu = impx_kolganu[i][j];
                        impyd = impy_kolgand[i][j];
                        impyu = impy_kolganu[i][j];
                        ed = e_kolgand[i][j];
                        eu = e_kolganu[i][j];
                    }
                    Godunov_solver_y(params_proc, md, impxd, impyd, ed, mu, impxu, impyu, eu, Fmd, Fimpxd, Fimpyd, Fed);

                    //up cell boundary flux
                    if (params_proc.s_type == 1) { //Godunov
                        md = m[i][j];
                        impxd = impx[i][j];
                        impyd = impy[i][j];
                        ed = e[i][j];
                        mu = m[i][j + 1];
                        impxu = impx[i][j + 1];
                        impyu = impy[i][j + 1];
                        eu = e[i][j + 1];
                    }
                    else if (params_proc.s_type == 2 || params_proc.s_type == 3){ // Kolgan/Rodionov Y
                        md = m_kolgand[i][j + 1];
                        mu = m_kolganu[i][j + 1];
                        impxd = impx_kolgand[i][j + 1];
                        impxu = impx_kolganu[i][j + 1];
                        impyd = impy_kolgand[i][j + 1];
                        impyu = impy_kolganu[i][j + 1];
                        ed = e_kolgand[i][j + 1];
                        eu = e_kolganu[i][j + 1];
                    }
                    Godunov_solver_y(params_proc, md, impxd, impyd, ed, mu, impxu, impyu, eu, Fmu, Fimpxu, Fimpyu, Feu);

                    dx = xc[i] - xc[i - 1];
                    dy = yc[j] - yc[j - 1];
                    //new values of cons
                    m_new[i][j] = m[i][j] - dt * ((Fmr - Fml) / dx + (Fmu - Fmd) / dy);
                    impx_new[i][j] = impx[i][j] - dt * ((Fimpxr - Fimpxl) / dx + (Fimpxu - Fimpxd) / dy);
                    impy_new[i][j] = impy[i][j] - dt * ((Fimpyr - Fimpyl) / dx + (Fimpyu - Fimpyd) / dy);
                    e_new[i][j] = e[i][j] - dt * ((Fer - Fel) / dx + (Feu - Fed) / dy);
                }
            }
            Boundary_x(params_proc, m_new, impx_new, impy_new, e_new);
            Boundary_y(params_proc, m_new, impx_new, impy_new, e_new);
            
            exchange_data(params_proc, myrank, px, py, m_new, impx_new, impy_new, e_new, buff_x, buff_y, comm, Status);

            for (int i = 0; i < params_proc.Nx + 2 * IS.fict; i++) {
                for (int j = 0; j < params_proc.Ny + 2 * IS.fict; j++) { 
                    impx[i][j] = impx_new[i][j];
                    impy[i][j] = impy_new[i][j];
                    m[i][j] = m_new[i][j];
                    e[i][j] = e_new[i][j];
                    cons_to_noncons(params_proc, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
                }
            }
        }
        else if (IS.s_type == 4) {
            dx = xc[2] - xc[1];
            dy = yc[2] - yc[1];
            for (int i = 0; i < params_proc.Nx + 2 * params_proc.fict; i++) {
                for(int j = 0; j < params_proc.Ny + 2 * params_proc.fict; j++){
                    U[0][0][i][j] = m[i][j];
                    U[0][1][i][j] = impx[i][j];
                    U[0][2][i][j] = impy[i][j];
                    U[0][3][i][j] = e[i][j];
                }
            }

            for (int s = 1; s < 4; s++) {
                WENO5(params_proc, U[s - 1], L, dx, dy);
                for (int k = 0; k < 4; k++) {
                    for (int i = params_proc.fict; i < params_proc.Nx + params_proc.fict; i++) {
                        for (int j = params_proc.fict; j < params_proc.Ny + params_proc.fict; j++) {
                            U[s][k][i][j] = rk[s][0] * U[0][k][i][j] + rk[s][1] * U[1][k][i][j] + rk[s][2] * U[2][k][i][j] - rk[s][3] * dt * L[k][i][j];
                        }
                    }
                }
                Boundary_x(params_proc, U[s][0], U[s][1], U[s][2], U[s][3]);
                Boundary_y(params_proc, U[s][0], U[s][1], U[s][2], U[s][3]);
                exchange_data(params_proc, myrank, px, py, U[s][0], U[s][1], U[s][2], U[s][3], buff_x, buff_y, comm, Status);
            }

            for (int i = 0; i < params_proc.Nx + 2 * params_proc.fict; i++) {
                for (int j = 0; j < params_proc.Ny + 2 * params_proc.fict; j++) {
                    m[i][j] = U[3][0][i][j];
                    impx[i][j] = U[3][1][i][j];
                    impy[i][j] = U[3][2][i][j];
                    e[i][j] = U[3][3][i][j];
                    cons_to_noncons(params_proc, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
                }
            }
        }
        else{
            std::cerr << "Invalid s_type" << std::endl;
        }

        if (step % params_proc.write_interval == 0) {
            calc_ei(params_proc, p, r, ei);
            write_out_p(params_proc, xc, yc, p, vx, vy, r, m, impx, impy, e, ei, step, time, out_dir, myrank, size, comm);
        }
    }
}

void Solve(InitialState& IS, std::string out_dir){
    vec1d x(IS.Nx + 1 + 2*IS.fict); 
    vec1d xc(IS.Nx + 1 + 2*IS.fict);
    vec1d y(IS.Ny + 1 + 2*IS.fict); 
    vec1d yc(IS.Ny + 1 + 2*IS.fict);

    vec2d p(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d vx(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d vy(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d r(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));   
    vec2d m(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d impx(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d impy(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d e(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d ei(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));

    vec2d m_new(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d impx_new(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d impy_new(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d e_new(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));

    // for Kolgan/Rodionov
    vec2d m_kolganl(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d impx_kolganl(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d impy_kolganl(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d e_kolganl(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d m_kolganr(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d impx_kolganr(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d impy_kolganr(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d e_kolganr(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));

    vec2d m_kolgand(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d impx_kolgand(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d impy_kolgand(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d e_kolgand(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d m_kolganu(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d impx_kolganu(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d impy_kolganu(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d e_kolganu(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));

    vec2d m05(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d impx05(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d impy05(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));
    vec2d e05(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict));

    double dm_x, dimpx_x, dimpy_x, de_x;
    double dm_y, dimpx_y, dimpy_y, de_y;
    

    // for WENO
    std::vector<vec3d> U(4, vec3d(4, vec2d(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict))));
    std::vector<vec1d> rk = { {1., 0., 0., 0.}, {1., 0., 0., 1.}, {0.75, 0.25, 0., 0.25}, {1. / 3., 0., 2. / 3., 2. / 3.} };
    vec3d L(4, vec2d(IS.Nx + 2*IS.fict, vec1d(IS.Ny + 2*IS.fict)));
    

    double Fml, Fimpxl, Fimpyl, Fel, Fmr, Fimpxr, Fimpyr, Fer;  //cell x fluxes
    double Fmd, Fimpxd, Fimpyd, Fed, Fmu, Fimpxu, Fimpyu, Feu;  //cell y fluxes

    double ml, impxl, impyl, el, mr, impxr, impyr, er; // left and right cell values of cons
    double md, impxd, impyd, ed, mu, impxu, impyu, eu; // down and up cell values of cons

    int step = 0;
    double time = 0.;
    double dt, dx, dy;

    x_grid(IS, xc, x);
    y_grid(IS, yc, y);
    Init(IS, xc, yc, p, vx, vy, r, m, impx, impy, e);
    Boundary_x(IS, m, impx, impy, e);
    Boundary_y(IS, m, impx, impy, e);
    calc_ei(IS, p, r, ei);

    if (fs::exists(out_dir)) {
        fs::remove_all(out_dir);
    }
    write_out(IS, xc, yc, p, vx, vy, r, m, impx, impy, e, ei, step, time, out_dir);

    while (time < IS.t_end) {
        dt = get_dt(IS, x, y, m, impx, impy, e);  // time step
        time += dt;
        step++;
        if(IS.s_type != 4){
            if (IS.s_type == 2 || IS.s_type == 3) {
                for (int i = IS.fict; i < IS.Nx + IS.fict + 1; ++i) {
                    for (int j = IS.fict; j < IS.Ny + IS.fict + 1; ++j) {
                        Kolgan_left(m[i][j], m[i - 1][j], m[i - 2][j], m[i + 1][j], xc[i], xc[i - 1], xc[i - 2], xc[i + 1], m_kolganl[i][j], m_kolganr[i][j], IS.mod_type);
                        Kolgan_left(impx[i][j], impx[i - 1][j], impx[i - 2][j], impx[i + 1][j], xc[i], xc[i - 1], xc[i - 2], xc[i + 1], impx_kolganl[i][j], impx_kolganr[i][j], IS.mod_type);
                        Kolgan_left(impy[i][j], impy[i - 1][j], impy[i - 2][j], impy[i + 1][j], xc[i], xc[i - 1], xc[i - 2], xc[i + 1], impy_kolganl[i][j], impy_kolganr[i][j], IS.mod_type);
                        Kolgan_left(e[i][j], e[i - 1][j], e[i - 2][j], e[i + 1][j], xc[i], xc[i - 1], xc[i - 2], xc[i + 1], e_kolganl[i][j], e_kolganr[i][j], IS.mod_type);
                        
                        if(i < IS.Nx + IS.fict){
                            Kolgan_right(m[i][j], m[i - 1][j], m[i + 1][j], m[i + 2][j], xc[i], xc[i - 1], xc[i + 1], xc[i + 2], m_kolganl[i + 1][j], m_kolganr[i + 1][j], IS.mod_type);
                            Kolgan_right(impx[i][j], impx[i - 1][j], impx[i + 1][j], impx[i + 2][j], xc[i], xc[i - 1], xc[i + 1], xc[i + 2], impx_kolganl[i + 1][j], impx_kolganr[i + 1][j], IS.mod_type);
                            Kolgan_right(impy[i][j], impy[i - 1][j], impy[i + 1][j], impy[i + 2][j], xc[i], xc[i - 1], xc[i + 1], xc[i + 2], impy_kolganl[i + 1][j], impy_kolganr[i + 1][j], IS.mod_type);
                            Kolgan_right(e[i][j], e[i - 1][j], e[i + 1][j], e[i + 2][j], xc[i], xc[i - 1], xc[i + 1], xc[i + 2], e_kolganl[i + 1][j], e_kolganr[i + 1][j], IS.mod_type);
                        }

                        Kolgan_left(m[i][j], m[i][j - 1], m[i][j - 2], m[i][j + 1], yc[j], yc[j - 1], yc[j - 2], yc[j + 1], m_kolgand[i][j], m_kolganu[i][j], IS.mod_type);
                        Kolgan_left(impx[i][j], impx[i][j - 1], impx[i][j - 2], impx[i][j + 1], yc[j], yc[j - 1], yc[j - 2], yc[j + 1], impx_kolgand[i][j], impx_kolganu[i][j], IS.mod_type);
                        Kolgan_left(impy[i][j], impy[i][j - 1], impy[i][j - 2], impy[i][j + 1], yc[j], yc[j - 1], yc[j - 2], yc[j + 1], impy_kolgand[i][j], impy_kolganu[i][j], IS.mod_type);
                        Kolgan_left(e[i][j], e[i][j - 1], e[i][j - 2], e[i][j + 1], yc[j], yc[j - 1], yc[j - 2], yc[j + 1], e_kolgand[i][j], e_kolganu[i][j], IS.mod_type);
                        
                        if(j < IS.Ny + IS.fict){
                            Kolgan_right(m[i][j], m[i][j - 1], m[i][j + 1], m[i][j + 2], yc[j], yc[j - 1], yc[j + 1], yc[j + 2], m_kolgand[i][j + 1], m_kolganu[i][j + 1], IS.mod_type);
                            Kolgan_right(impx[i][j], impx[i][j - 1], impx[i][j + 1], impx[i][j + 2], yc[j], yc[j - 1], yc[j + 1], yc[j + 2], impx_kolgand[i][j + 1], impx_kolganu[i][j + 1], IS.mod_type);
                            Kolgan_right(impy[i][j], impy[i][j - 1], impy[i][j + 1], impy[i][j + 2], yc[j], yc[j - 1], yc[j + 1], yc[j + 2], impy_kolgand[i][j + 1], impy_kolganu[i][j + 1], IS.mod_type);
                            Kolgan_right(e[i][j], e[i][j - 1], e[i][j + 1], e[i][j + 2], yc[j], yc[j - 1], yc[j + 1], yc[j + 2], e_kolgand[i][j + 1], e_kolganr[i][j + 1], IS.mod_type);
                        }
                    }
                }
            }
            if (IS.s_type == 3) {
                for (int i = IS.fict; i < IS.Nx + IS.fict; ++i) {
                    for (int j = IS.fict; j < IS.Ny + IS.fict; ++j) {
                        // Predictor
                        cons_flux_x(IS, m_kolganr[i][j], impx_kolganr[i][j], impy_kolganr[i][j], e_kolganr[i][j], Fml, Fimpxl, Fimpyl, Fel);
                        cons_flux_x(IS, m_kolganl[i + 1][j], impx_kolganl[i + 1][j], impy_kolganl[i + 1][j], e_kolganl[i + 1][j], Fmr, Fimpxr, Fimpyr, Fer);
                        cons_flux_y(IS, m_kolganu[i][j], impx_kolganu[i][j], impy_kolganu[i][j], e_kolganu[i][j], Fmd, Fimpxd, Fimpyd, Fed);
                        cons_flux_y(IS, m_kolgand[i][j + 1], impx_kolgand[i][j + 1], impy_kolgand[i][j + 1], e_kolgand[i][j + 1], Fmu, Fimpxu, Fimpyu, Feu);
                        m05[i][j] = m[i][j] - dt * ((Fmr - Fml) / (xc[i] - xc[i - 1]) + (Fmu - Fmd) / (yc[j] - yc[j - 1]));
                        impx05[i][j] = impx[i][j] - dt * ((Fimpxr - Fimpxl) / (xc[i] - xc[i - 1]) + (Fimpxu - Fimpxd) / (yc[j] - yc[j - 1]));
                        impy05[i][j] = impy[i][j] - dt * ((Fimpyr - Fimpyl) / (xc[i] - xc[i - 1]) + (Fimpyu - Fimpyd) / (yc[j] - yc[j - 1]));
                        e05[i][j] = e[i][j] - dt * ((Fer - Fel) / (xc[i] - xc[i - 1]) + (Feu - Fed) / (yc[j] - yc[j - 1]));
                    }
				}
				for (int i = IS.fict; i < IS.Nx + IS.fict + 1; ++i) {
                    for (int j = IS.fict; j < IS.Ny + IS.fict; ++j) {
                        // Corrector x
                        dm_x = (m_kolganl[i + 1][j] - m_kolganr[i][j]);
                        dimpx_x = (impx_kolganl[i + 1][j] - impx_kolganr[i][j]);
                        dimpy_x = (impy_kolganl[i + 1][j] - impy_kolganr[i][j]);
                        de_x = (e_kolganl[i + 1][j] - e_kolganr[i][j]);

                        m_kolganr[i][j] = 0.5 * (m[i][j] + m05[i][j]) - 0.5 * dm_x;
                        impx_kolganr[i][j] = 0.5 * (impx[i][j] + impx05[i][j]) - 0.5 * dimpx_x;
                        impy_kolganr[i][j] = 0.5 * (impy[i][j] + impy05[i][j]) - 0.5 * dimpy_x;
                        e_kolganr[i][j] = 0.5 * (e[i][j] + e05[i][j]) - 0.5 * de_x;
                        m_kolganl[i + 1][j] = 0.5 * (m[i][j] + m05[i][j]) + 0.5 * dm_x;
                        impx_kolganl[i + 1][j] = 0.5 * (impx[i][j] + impx05[i][j]) + 0.5 * dimpx_x;
                        impy_kolganl[i + 1][j] = 0.5 * (impy[i][j] + impy05[i][j]) + 0.5 * dimpy_x;
                        e_kolganl[i + 1][j] = 0.5 * (e[i][j] + e05[i][j]) + 0.5 * de_x;

                        // Corrector y
                        dm_y = (m_kolgand[i][j + 1] - m_kolganu[i][j]);
                        dimpx_y = (impx_kolgand[i][j + 1] - impx_kolganu[i][j]);
                        dimpy_y = (impy_kolgand[i][j + 1] - impy_kolganu[i][j]);
                        de_y = (e_kolgand[i][j + 1] - e_kolganu[i][j]);

                        m_kolganu[i][j] = 0.5 * (m[i][j] + m05[i][j]) - 0.5 * dm_y;
                        impx_kolganu[i][j] = 0.5 * (impx[i][j] + impx05[i][j]) - 0.5 * dimpx_y;
                        impy_kolganu[i][j] = 0.5 * (impy[i][j] + impy05[i][j]) - 0.5 * dimpy_y;
                        e_kolganu[i][j] = 0.5 * (e[i][j] + e05[i][j]) - 0.5 * de_y;
                        m_kolgand[i][j + 1] = 0.5 * (m[i][j] + m05[i][j]) + 0.5 * dm_y;
                        impx_kolgand[i][j + 1] = 0.5 * (impx[i][j] + impx05[i][j]) + 0.5 * dimpx_y;
                        impy_kolgand[i][j + 1] = 0.5 * (impy[i][j] + impy05[i][j]) + 0.5 * dimpy_y;
                        e_kolgand[i][j + 1] = 0.5 * (e[i][j] + e05[i][j]) + 0.5 * de_y;
                    }
                }
            }
            for (int i = IS.fict; i < IS.Nx + IS.fict; ++i) {
                for (int j = IS.fict; j < IS.Ny + IS.fict; ++j) {
                    //left cell bounday flux
                    if (IS.s_type == 1) {  //Godunov X
                        ml = m[i - 1][j];
                        impxl = impx[i - 1][j];
                        impyl = impy[i - 1][j];
                        el = e[i - 1][j];
                        mr = m[i][j];
                        impxr = impx[i][j];
                        impyr = impy[i][j];
                        er = e[i][j];
                    }
                    else if (IS.s_type == 2 || IS.s_type == 3){
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
                    
                    //right cell boundary flux
                    if (IS.s_type == 1) { //Godunov X
                        ml = m[i][j];
                        impxl = impx[i][j];
                        impyl = impy[i][j];
                        el = e[i][j];
                        mr = m[i + 1][j];
                        impxr = impx[i + 1][j];
                        impyr = impy[i + 1][j];
                        er = e[i + 1][j];
                    }
                    else if (IS.s_type == 2 || IS.s_type == 3){
                        ml = m_kolganl[i + 1][j];
                        mr = m_kolganr[i + 1][j];
                        impxl = impx_kolganl[i + 1][j];
                        impxr = impx_kolganr[i + 1][j];
                        impyl = impy_kolganl[i + 1][j];
                        impyr = impy_kolganr[i + 1][j];
                        el = e_kolganl[i + 1][j];
                        er = e_kolganr[i + 1][j];
                    }
                    Godunov_solver_x(IS, ml, impxl, impyl, el, mr, impxr, impyr, er, Fmr, Fimpxr, Fimpyr, Fer);

                    //down cell bounday flux
                    if (IS.s_type == 1) {  //Godunov Y
                        md = m[i][j - 1];
                        impxd = impx[i][j - 1];
                        impyd = impy[i][j - 1];
                        ed = e[i][j - 1];
                        mu = m[i][j];
                        impxu = impx[i][j];
                        impyu = impy[i][j];
                        eu = e[i][j];
                    }
                    else if (IS.s_type == 2 || IS.s_type == 3){
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

                    //up cell boundary flux
                    if (IS.s_type == 1) { //Godunov
                        md = m[i][j];
                        impxd = impx[i][j];
                        impyd = impy[i][j];
                        ed = e[i][j];
                        mu = m[i][j + 1];
                        impxu = impx[i][j + 1];
                        impyu = impy[i][j + 1];
                        eu = e[i][j + 1];
                    }
                    else if (IS.s_type == 2 || IS.s_type == 3){
                        md = m_kolgand[i][j + 1];
                        mu = m_kolganu[i][j + 1];
                        impxd = impx_kolgand[i][j + 1];
                        impxu = impx_kolganu[i][j + 1];
                        impyd = impy_kolgand[i][j + 1];
                        impyu = impy_kolganu[i][j + 1];
                        ed = e_kolgand[i][j + 1];
                        eu = e_kolganu[i][j + 1];
                    }
                    Godunov_solver_y(IS, md, impxd, impyd, ed, mu, impxu, impyu, eu, Fmu, Fimpxu, Fimpyu, Feu);

                    if (IS.s_type != 0) {
                        dx = xc[i] - xc[i - 1];
                        dy = yc[j] - yc[j - 1];
                        //new values of cons
                        m_new[i][j] = m[i][j] - dt * ((Fmr - Fml) / dx + (Fmu - Fmd) / dy);
                        impx_new[i][j] = impx[i][j] - dt * ((Fimpxr - Fimpxl) / dx + (Fimpxu - Fimpxd) / dy);
                        impy_new[i][j] = impy[i][j] - dt * ((Fimpyr - Fimpyl) / dx + (Fimpyu - Fimpyd) / dy);
                        e_new[i][j] = e[i][j] - dt * ((Fer - Fel) / dx + (Feu - Fed) / dy);
                    }
                }
            }
            Boundary_x(IS, m_new, impx_new, impy_new, e_new);
            Boundary_y(IS, m_new, impx_new, impy_new, e_new);


            for (int i = 0; i < IS.Nx + 2 * IS.fict; i++) {
                for (int j = 0; j < IS.Ny + 2 * IS.fict; j++) { 
                    impx[i][j] = impx_new[i][j];
                    impy[i][j] = impy_new[i][j];
                    m[i][j] = m_new[i][j];
                    e[i][j] = e_new[i][j];
                    cons_to_noncons(IS, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
                }
            }
        }
        else if (IS.s_type == 4) {
            dx = xc[2] - xc[1];
            dy = yc[2] - yc[1];
            for (int i = 0; i < IS.Nx + 2 * IS.fict; i++) {
                for(int j = 0; j < IS.Ny + 2 * IS.fict; j++){
                    U[0][0][i][j] = m[i][j];
                    U[0][1][i][j] = impx[i][j];
                    U[0][2][i][j] = impy[i][j];
                    U[0][3][i][j] = e[i][j];
                }
            }

            for (int s = 1; s < 4; s++) {
                WENO5(IS, U[s - 1], L, dx, dy);
                for (int k = 0; k < 4; k++) {
                    for (int i = IS.fict; i < IS.Nx + IS.fict; i++) {
                        for (int j = IS.fict; j < IS.Ny + IS.fict; j++) {
                            U[s][k][i][j] = rk[s][0] * U[0][k][i][j] + rk[s][1] * U[1][k][i][j] + rk[s][2] * U[2][k][i][j] - rk[s][3] * dt * L[k][i][j];
                        }
                    }
                }
                Boundary_x(IS, U[s][0], U[s][1], U[s][2], U[s][3]);
                Boundary_y(IS, U[s][0], U[s][1], U[s][2], U[s][3]);
            }

            for (int i = 0; i < IS.Nx + 2 * IS.fict; i++) {
                for (int j = 0; j < IS.Ny + 2 * IS.fict; j++) {
                    m[i][j] = U[3][0][i][j];
                    impx[i][j] = U[3][1][i][j];
                    impy[i][j] = U[3][2][i][j];
                    e[i][j] = U[3][3][i][j];
                    cons_to_noncons(IS, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
                }
            }
        }
        else{
            std::cerr << "Invalid s_type" << std::endl;
        }

        if (step % IS.write_interval == 0) {
            calc_ei(IS, p, r, ei);
            write_out(IS, xc, yc, p, vx, vy, r, m, impx, impy, e, ei, step, time, out_dir);
        }
    }
}

int main(int argc, char* argv[]){
    int myrank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    double begin, end;
    
    InitialState IS;

    read_params(IS, "Init.txt");
	IS.s_type = 1;
	IS.fict = 1;

    MPI_Barrier(MPI_COMM_WORLD);
    begin = MPI_Wtime();
    Multiproc_solve(IS, "Godunov_2D_p", myrank, size, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    if(myrank == 0) printf("Godunov_2D_p time: %.5f\n\n", end - begin);

	IS.s_type = 2;
	IS.fict = 2;
    MPI_Barrier(MPI_COMM_WORLD);
    begin = MPI_Wtime();
    Multiproc_solve(IS, "GK_2D_p", myrank, size, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    if(myrank == 0) printf("GK_2D_p time: %.5f\n\n", end - begin);
    
	IS.s_type = 3;
	IS.fict = 2;
    MPI_Barrier(MPI_COMM_WORLD);
    begin = MPI_Wtime();
    Multiproc_solve(IS, "GKR_2D_p", myrank, size, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    if(myrank == 0) printf("GKR_2D_p time: %.5f\n\n", end - begin);

	IS.s_type = 4;
	IS.fict = 3;
    MPI_Barrier(MPI_COMM_WORLD);
    begin = MPI_Wtime();
    Multiproc_solve(IS, "WENO_2D_p", myrank, size, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    if(myrank == 0) printf("WENO_2D_p time: %.5f\n\n", end - begin);
    
    /*
    if(myrank == 0){
        begin = MPI_Wtime();
        read_params(IS, "W.txt");
        Solve(IS, "resultW");
        end = MPI_Wtime();
        printf("time sequently: %.5f\n\n", end - begin);
    }
    */
    
    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    return 0;
}