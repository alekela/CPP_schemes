#include <iostream>
#include<fstream>
#include <string>
#include <vector>
#include <map>
#include <functional>


struct Point {
    double x, y, z;
};


class PrePost {
public:
    std::map<std::string, std::vector<double>> input;
    std::vector<std::vector<double>> grid;
    std::map<std::string, std::vector<std::vector<std::vector<double>>>> res;


    PrePost() {};

    std::vector<std::string> split(std::string s, std::string delimiter) {
        size_t pos_start = 0, pos_end, delim_len = delimiter.length();
        std::string token;
        std::vector<std::string> res;

        while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
            token = s.substr(pos_start, pos_end - pos_start);
            pos_start = pos_end + delim_len;
            res.push_back(token);
        }

        res.push_back(s.substr(pos_start));
        return res;
    }

    void ReadFile(std::string name) {
        std::ifstream file;
        std::string line;
        file.open(name);
        while (std::getline(file, line))
        {
            if (line[0] != '#') {
                std::vector<std::string> v = split(line, ":");
                std::vector<double> tmp;
                for (auto i : split(v[1], ",")) {
                    tmp.push_back(stod(i));
                }
                input[v[0]] = tmp;
            }
        }
    };

    std::map<std::string, std::vector<double>> GetInputDict() {
        return input;
    }

    std::vector<std::vector<double>> GetInitGrid() {
        for (int i = 0; i < input["Dimension"][0]; i++) {
            std::vector<double> tmp;
            if (i == 0) {
                for (double x = input["x_start"][0]; x <= input["x_end"][0]; x += input["dx"][0]) {
                    tmp.push_back(x);
                }
            }
            else if (i == 1) {
                for (double y = input["y_start"][0]; y <= input["y_end"][0]; y += input["dy"][0]) {
                    tmp.push_back(y);
                }
            }
            else if (i == 2) {
                for (double z = input["z_start"][0]; z <= input["z_end"][0]; z += input["dz"][0]) {
                    tmp.push_back(z);
                }
            }
            grid.push_back(tmp);
        }
        return grid;
    }

    std::map<std::string, std::vector<std::vector<std::vector<double>>>> GetInitRes(std::vector<std::string> params, std::vector<std::function<double(Point)>> funcs, std::vector<std::vector<double>> grid) {
        for (int num = 0; num < params.size(); num++) {
            Point point;
            std::vector<std::vector<std::vector<double>>> tmp;
            for (int i = 0; i < grid[0].size(); i++) {
                std::vector<std::vector<double>> tmp2;
                for (int j = 0; j < grid[1].size(); j++) {
                    std::vector<double> tmp3;
                    for (int k = 0; k < grid[2].size(); k++) {
                        point.x = grid[0][i];
                        point.y = grid[1][j];
                        point.z = grid[2][k];
                        tmp3.push_back(funcs[num](point));
                    }
                    tmp2.push_back(tmp3);
                }
                tmp.push_back(tmp2);
            }
            res[params[num]] = tmp;
        }
        return res;
    }

};