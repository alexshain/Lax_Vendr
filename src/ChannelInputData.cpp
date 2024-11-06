#include "ChannelInputData.h"

ChannelInputData::ChannelInputData(Grid grid, MediumParameters left_parameters, MediumParameters right_parameters) :
    grid_(grid), 
    left_parameters_(left_parameters),
    right_parameters_(right_parameters) {}

Grid ChannelInputData::getGrid() const {
    return grid_;
}

MediumParameters ChannelInputData::getLeftParameters() const {
    return left_parameters_;
}

MediumParameters ChannelInputData::getRightParameters() const {
    return right_parameters_;
}

double ChannelInputData::crossSection(double x) const {
    //return 1;
    return 0.5 + 0.5*pow((1 - x/0.5), 2);
}