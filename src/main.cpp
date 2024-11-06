#include "ChannelInputData.h"
#include "ChannelOutputData.h"
#include "ChannelSolver.h"
#include "Grid.h"

#include <fstream>
#include <iostream>

int main() {
    Grid grid;
    grid.spatial_step = 0.11;
    grid.left_bound = 0;
    grid.right_bound = 1;
    grid.time_step = 0.0001;
    grid.total_time = 10;
    MediumParameters left_parameters;
    left_parameters.density = 1;
    left_parameters.velosity = 1.0237498;
    left_parameters.pressure = 8;
    MediumParameters right_parameters;
    right_parameters.density = 0.193388;
    right_parameters.velosity = 5.2937598;
    right_parameters.pressure = 0.8018469;
    ChannelInputData input = ChannelInputData{grid, left_parameters, right_parameters};
    size_t size = static_cast<size_t>(grid.right_bound/grid.spatial_step) * (grid.total_time/grid.time_step);
    ChannelOutputData output = ChannelOutputData{size};
    ChannelSolver solver = ChannelSolver{input, output};

    ChannelOutputData output2 = solver.solve();

    std::vector<double> density;
    std::vector<double> velosity;
    std::vector<double> pressure;

    for(size_t i = 100; i < 200; i++) {
        //std::cout << output2.getDensityVector()[i] << std::endl;
        density.push_back(output2.getDensityVector()[i]);
        //std::cout << density[i] << std::endl;
        velosity.push_back(output2.getVelosityVector()[i]);
        pressure.push_back(output2.getPressureVector()[i]);
    }
    std::ofstream fout1("density.raw", std::ios_base::binary);
    fout1.write((char*)density.data(), sizeof(double) * density.size());
    fout1.close();
    std::ofstream fout2("velosity.raw", std::ios_base::binary);
    fout2.write((char*)velosity.data(), sizeof(double) * velosity.size());
    fout2.close();
    std::ofstream fout3("pressure.raw", std::ios_base::binary);
    fout3.write((char*)pressure.data(), sizeof(double) * pressure.size());
    fout3.close();
    return 0;
}