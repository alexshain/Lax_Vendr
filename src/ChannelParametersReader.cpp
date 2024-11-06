#include "ChannelParametersReader.h"

/*void ChannelParametersReader::read(std::filesystem::path medium_file, std::filesystem::path grid_file) {
    input_data_ = ChannelInputData{readGrid(grid_file), readMedium(medium_file)};
}*/

MediumParameters ChannelParametersReader::readMedium(std::filesystem::path &medium_file) {
    std::ifstream file(medium_file);
    if (!file.is_open()) {
        throw std::runtime_error("Error: File " + medium_file.string() + " does not exist or cannot be opened.");
    }
    MediumParameters parameters;
    file >> parameters.density
         >> parameters.velosity
         >> parameters.pressure;
    file.close();
    return parameters;
}

Grid ChannelParametersReader::readGrid(std::filesystem::path &grid_file) {
    std::ifstream file(grid_file);
    if (!file.is_open()) {
        throw std::runtime_error("Error: File " + grid_file.string() + " does not exist or cannot be opened.");
    }
    Grid grid;
    file >> grid.spatial_step
         >> grid.left_bound
         >> grid.right_bound
         >> grid.time_step
         >> grid.total_time;
    file.close();
    return grid;
}

ChannelInputData ChannelParametersReader::getInputData() {
    return input_data_;
}