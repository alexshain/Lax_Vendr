#pragma once

#include "IChannelSolver.h"
#include "ChannelInputData.h"
#include "ChannelOutputData.h"

class ChannelSolver: public IChannelSolver {
    ChannelInputData channel_input_data_;
    ChannelOutputData channel_output_data_;

    size_t number_of_nodes_;
    size_t number_of_time_steps_;

    double tau_;
    double h_;

    const int MACH_IN_CRITICAL = 1;

    public:
    ChannelSolver(const ChannelInputData &channel_input_data, const ChannelOutputData &channel_output_data);

    ChannelOutputData solve();

    private:
    size_t getIndMinCrossSection();
    double getSoundVelocity(MediumParameters parameters);
    double fMassPart(MediumParameters &parameters);
    double fMomentumPart(MediumParameters &parameters);
    double fEnergyPart(MediumParameters &parameters, double energy);
    void setInitData(size_t x_ind_min_cross);
    MediumParameters firstStepSolve(size_t prev_ind, size_t next_ind, size_t x_ind);
    void boundary(size_t start, size_t ind_min_cross, size_t end);
};