#include <iostream>
#include<fstream>
#include <string>
#include <vector>
#include <map>
#include <functional>

class PostProc {
public:
    PostProc() {};
    void WriteCSV(std::string name, std::map<std::string, std::vector<std::vector<std::vector<double>>>> res, std::vector<std::vector<double>> grid) {
        std::ofstream outfile(name);
        std::string tmp = "X,Y,Z,";
        for (auto outputs : res) {
            tmp += outputs.first + ',';
        }
        tmp.pop_back();
        outfile << tmp << '\n';
        for (int i = 0; i < grid[0].size(); i++) {
            for (int j = 0; j < grid[1].size(); j++) {
                for (int k = 0; k < grid[2].size(); k++) {
                    tmp = "";
                    tmp += std::to_string(grid[0][i]) + ',';
                    tmp += std::to_string(grid[1][j]) + ',';
                    tmp += std::to_string(grid[2][k]) + ',';
                    for (auto outputs : res) {
                        tmp += std::to_string(outputs.second[i][j][k]) + ',';
                    }
                    tmp.pop_back();
                    outfile << tmp << '\n';
                }
            }
        }
    }
};
