#pragma once

#include "ChannelInputData.h"

class ChannelParametersReader {
    ChannelInputData input_data_;

    public:
    void read(std::filesystem::path medium_file, std::filesystem::path grid_file);
    ChannelInputData getInputData();

    private:
    MediumParameters readMedium(std::filesystem::path &medium_file);
    Grid readGrid(std::filesystem::path &grid_file);
};