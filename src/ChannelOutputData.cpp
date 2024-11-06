#include <iostream>

#include "ChannelOutputData.h"

ChannelOutputData::ChannelOutputData(size_t data_size) : 
    density_(data_size), 
    velosity_(data_size), 
    pressure_(data_size) {}

double ChannelOutputData::getDensity(size_t ind) const {
    return density_[ind];
}

double ChannelOutputData::getVelosity(size_t ind) const {
    return velosity_[ind];
}

double ChannelOutputData::getPressure(size_t ind) const {
   return pressure_[ind];
}

std::vector<double> ChannelOutputData::getDensityVector() const {
    return density_;
}

std::vector<double> ChannelOutputData::getVelosityVector() const {
    return velosity_;
}

std::vector<double> ChannelOutputData::getPressureVector() const {
    return pressure_;
}

MediumParameters ChannelOutputData::getParams(size_t ind) const {
    MediumParameters params;
    params.density = density_[ind];
    params.velosity = velosity_[ind];
    params.pressure = pressure_[ind];
    return params;
}

void ChannelOutputData::setFlowParameters(MediumParameters parameters, size_t ind) {
    setNextDensity(parameters.density, ind);
    setNextVelosity(parameters.velosity, ind);
    setNextPressure(parameters.pressure, ind);
}

void ChannelOutputData::setNextDensity(double next_density, size_t next_ind) {
    //std::cout << next_density << std::endl;
    density_[next_ind] = next_density; 
}

void ChannelOutputData::setNextVelosity(double next_mach_number, size_t next_ind) {
    velosity_[next_ind] = next_mach_number;
}

void ChannelOutputData::setNextPressure(double next_pressure, size_t next_ind) {
    pressure_[next_ind] = next_pressure;
}