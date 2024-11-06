#pragma once

#include <vector>

#include "ChannelInputData.h"

class ChannelOutputData {
    std::vector<double> density_;
    std::vector<double> velosity_;
    std::vector<double> pressure_;

    public:
    ChannelOutputData(size_t data_size);

    double getDensity(size_t ind) const;
    double getVelosity(size_t ind) const;
    double getPressure(size_t ind) const;
    std::vector<double> getDensityVector() const;
    std::vector<double> getVelosityVector() const;
    std::vector<double> getPressureVector() const;
    MediumParameters getParams(size_t ind) const;

    void setFlowParameters(physical_values::MediumParameters parameters, size_t ind);

    private:
    void setNextDensity(double next_density, size_t ind);
    void setNextVelosity(double next_mach_number, size_t ind);
    void setNextPressure(double next_pressure, size_t ind);
};