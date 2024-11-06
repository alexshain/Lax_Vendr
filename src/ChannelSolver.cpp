#include <cmath>
#include <iostream>

#include "ChannelSolver.h"

ChannelSolver::ChannelSolver(const ChannelInputData &channel_input_data, const ChannelOutputData &channel_output_data) : 
    channel_input_data_(channel_input_data),
    channel_output_data_(channel_output_data), 
    number_of_nodes_((channel_input_data.getGrid().right_bound 
                    - channel_input_data.getGrid().left_bound)
                    / channel_input_data.getGrid().spatial_step + 1),
    number_of_time_steps_(channel_input_data.getGrid().total_time 
                        / channel_input_data.getGrid().time_step + 1),
    tau_(channel_input_data_.getGrid().time_step),
    h_(channel_input_data_.getGrid().spatial_step) {}

ChannelOutputData ChannelSolver::solve() {
    std::cout << "Начало solve" << std::endl;
    size_t x_ind_min_cross = getIndMinCrossSection();
    setInitData(x_ind_min_cross);
    MediumParameters parameters;
    MediumParameters next_first_step_params;
    MediumParameters prev_first_step_params;
    std::cout << "Входим в цикл" << std::endl;
    for(size_t t_ind = 1; t_ind < 50; t_ind++) {
        boundary(t_ind * number_of_nodes_, t_ind * number_of_nodes_ + x_ind_min_cross, t_ind * number_of_nodes_ + number_of_nodes_ - 1);
        std::cout << channel_output_data_.getDensity((t_ind-1) * number_of_nodes_) << "   ";
        std::cout << channel_output_data_.getVelosity((t_ind-1) * number_of_nodes_) << "   ";
        std::cout << channel_output_data_.getPressure((t_ind-1) * number_of_nodes_) << "   " << std::endl;
        for(size_t x_ind = 1; x_ind < number_of_nodes_ - 1; x_ind++) {
            size_t prev_ind = (t_ind-1) * number_of_nodes_ + x_ind;
            size_t next_ind = t_ind * number_of_nodes_ + x_ind;
            std::cout << channel_output_data_.getDensity(prev_ind) << "   ";
            std::cout << channel_output_data_.getVelosity(prev_ind) << "   ";
            std::cout << channel_output_data_.getPressure(prev_ind) << "   " << std::endl;
            next_first_step_params = firstStepSolve(prev_ind, prev_ind + 1, x_ind);
            prev_first_step_params = firstStepSolve(prev_ind - 1, prev_ind, x_ind);
            parameters.density = channel_input_data_.crossSection(x_ind*h_)*(channel_output_data_.getDensity(prev_ind)
                               - tau_*(fMassPart(next_first_step_params) - fMassPart(prev_first_step_params))/h_);
            parameters.velosity = (channel_output_data_.getVelosity(prev_ind)*channel_output_data_.getDensity(prev_ind)  
                                - tau_*(fMomentumPart(next_first_step_params) - fMomentumPart(prev_first_step_params)
                                + channel_output_data_.getDensity(prev_ind)*0.5*(channel_input_data_.crossSection((x_ind+1)*h_)
                                - channel_input_data_.crossSection((x_ind-1)*h_)) / channel_input_data_.crossSection(x_ind*h_))/h_) 
                                / parameters.density;      
            double energy = channel_input_data_.crossSection(x_ind*h_)*((channel_output_data_.getDensity(prev_ind)*getEnergy(channel_output_data_.getParams(prev_ind))
                          - tau_*(fEnergyPart(next_first_step_params, energy) - fEnergyPart(prev_first_step_params, energy))/h_)
                          / parameters.density);
            setPressureFromEnergy(parameters, energy);
            parameters.pressure *= channel_input_data_.crossSection(x_ind*h_);
            channel_output_data_.setFlowParameters(parameters, next_ind);
        }
        parameters.velosity = getSoundVelocity(parameters);
        channel_output_data_.setFlowParameters(parameters, t_ind * number_of_nodes_ + x_ind_min_cross);
        std::cout << channel_output_data_.getDensity((t_ind-1) * number_of_nodes_ + number_of_nodes_ - 1) << "   ";
        std::cout << channel_output_data_.getVelosity((t_ind-1) * number_of_nodes_ + number_of_nodes_ - 1) << "   ";
        std::cout << channel_output_data_.getPressure((t_ind-1) * number_of_nodes_ + number_of_nodes_ - 1) << "   " << std::endl;
        std::cout << std::endl;
    }
    return channel_output_data_;
}

void ChannelSolver::setInitData(size_t x_ind_min_cross) {
    MediumParameters parameters;
    double left_density = channel_input_data_.crossSection(0)*channel_input_data_.getLeftParameters().density;
    double left_velosity = channel_input_data_.getLeftParameters().velosity;
    double left_pressure = channel_input_data_.crossSection(0)*channel_input_data_.getLeftParameters().pressure;

    double right_density = channel_input_data_.crossSection(0)*channel_input_data_.getRightParameters().density;
    double right_velosity = channel_input_data_.getRightParameters().velosity;
    double right_pressure = channel_input_data_.crossSection(0)*channel_input_data_.getRightParameters().pressure;

    for(size_t x_ind = 0; x_ind < number_of_nodes_; x_ind++) {
        parameters.density = left_density + x_ind*(right_density - left_density)/(number_of_nodes_-1);
        parameters.velosity = left_velosity + x_ind*(right_velosity - left_velosity)/(number_of_nodes_-1);
        parameters.pressure = left_pressure + x_ind*(right_pressure - left_pressure)/(number_of_nodes_-1);
        channel_output_data_.setFlowParameters(parameters, x_ind);
    }

    parameters.velosity = getSoundVelocity(parameters);
    channel_output_data_.setFlowParameters(parameters, x_ind_min_cross);
}

MediumParameters ChannelSolver::firstStepSolve(size_t prev_ind, size_t next_ind, size_t x_ind) {
    MediumParameters intermediate_parameters;
    MediumParameters prev_params;
    MediumParameters next_params;
    next_params.density = channel_input_data_.crossSection(x_ind*h_)*channel_output_data_.getDensity(next_ind);
    prev_params.density = channel_input_data_.crossSection(x_ind*h_)*channel_output_data_.getDensity(prev_ind);
    next_params.velosity = channel_output_data_.getVelosity(next_ind);
    prev_params.velosity = channel_output_data_.getVelosity(prev_ind);
    next_params.pressure = channel_input_data_.crossSection(x_ind*h_)*channel_output_data_.getPressure(next_ind);
    prev_params.pressure = channel_input_data_.crossSection(x_ind*h_)*channel_output_data_.getPressure(prev_ind);
    intermediate_parameters.density = channel_input_data_.crossSection(x_ind*h_)*((0.5*(next_params.density + prev_params.density) 
                                    - 0.5*tau_*(fMassPart(next_params) - fMassPart(prev_params))/h_));
    intermediate_parameters.velosity = (0.5*(fMassPart(next_params) + fMassPart(prev_params)) 
                                     - 0.5*tau_*(fMomentumPart(next_params) - fMomentumPart(prev_params))/h_)
                                     / intermediate_parameters.density;
    double energy = channel_input_data_.crossSection(x_ind*h_)*((0.5*(next_params.density*getEnergy(next_params) + prev_params.density*getEnergy(prev_params)) 
                  - 0.5*tau_*(fEnergyPart(next_params, getEnergy(next_params)) - fEnergyPart(prev_params, getEnergy(prev_params)))/h_)
                  / intermediate_parameters.density);
    setPressureFromEnergy(intermediate_parameters, energy);
    intermediate_parameters.pressure *= channel_input_data_.crossSection(x_ind*h_);
    return intermediate_parameters; 
}

double ChannelSolver::fMassPart(MediumParameters &parameters) {
    return parameters.density * parameters.velosity;
}

double ChannelSolver::fMomentumPart(MediumParameters &parameters) {
    return parameters.density*pow(parameters.velosity, 2) + parameters.pressure;
}

double ChannelSolver::fEnergyPart(MediumParameters &parameters, double energy) {
    return parameters.velosity * (parameters.density*energy + parameters.pressure);
}

void ChannelSolver::boundary(size_t start, size_t ind_min_cross, size_t end) {
    MediumParameters left_parameters;
    left_parameters.density = channel_input_data_.getLeftParameters().density;
    left_parameters.velosity = channel_input_data_.getLeftParameters().velosity;
    left_parameters.pressure = channel_input_data_.getLeftParameters().pressure;
    channel_output_data_.setFlowParameters(left_parameters, start);
    //left_parameters = channel_output_data_.getParams(ind_min_cross);
    //left_parameters.velosity = getSoundVelocity(left_parameters);

    //std::cout << "left_parameters.velosity: " << left_parameters.velosity << std::endl;

    //channel_output_data_.setFlowParameters(left_parameters, ind_min_cross);

    MediumParameters right_parameters;
    right_parameters.density = channel_input_data_.getRightParameters().density;
    right_parameters.velosity = channel_input_data_.getRightParameters().velosity;
    right_parameters.pressure = channel_input_data_.getRightParameters().pressure;
    channel_output_data_.setFlowParameters(right_parameters, end);
}

size_t ChannelSolver::getIndMinCrossSection() {
    double minCross = channel_input_data_.crossSection(0); 
    size_t ind_in_min_section = 0; 
    for(size_t x_ind = 0; x_ind < number_of_nodes_; x_ind++) {
        if(channel_input_data_.crossSection(x_ind*channel_input_data_.getGrid().spatial_step) < minCross) {
            minCross = channel_input_data_.crossSection(x_ind*channel_input_data_.getGrid().spatial_step);
            ind_in_min_section = x_ind;
        }
    }
    return ind_in_min_section;
}

double ChannelSolver::getSoundVelocity(MediumParameters parameters) {
    double pressure = parameters.pressure;
    double density = parameters.density;
    return sqrt(ADIABATIC_EXPONENT*pressure/density);
}