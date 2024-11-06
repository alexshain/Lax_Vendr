#pragma once

#include <filesystem>
#include <vector>
#include <fstream>

#include "PhysicalValues.h"
#include "Grid.h"

using namespace physical_values;

class ChannelInputData {
    Grid grid_;
    MediumParameters left_parameters_;
    MediumParameters right_parameters_;

    public:
    ChannelInputData(Grid grid, MediumParameters left_parameters, MediumParameters right_parameters);
    Grid getGrid() const;
    MediumParameters getLeftParameters() const;
    MediumParameters getRightParameters() const;
    double crossSection(double x) const;
};

