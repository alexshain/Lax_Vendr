#pragma once

namespace physical_values {
    const double ADIABATIC_EXPONENT = 1.4;

    struct MediumParameters {
        double density;
        double velosity;
        double pressure;
    };

    inline double getEnergy(MediumParameters parameters) {
        double p = parameters.pressure;
        double u = parameters.velosity;
        double rho = parameters.density;
        return p/(rho*(ADIABATIC_EXPONENT - 1)) + pow(u, 2)/2;
    }

    inline double getEnthalpy(MediumParameters parameters) {
        return (getEnergy(parameters) + parameters.pressure)/parameters.density;
    }

    inline void setPressureFromEnergy(MediumParameters &parameters, double energy) {
        double rho = parameters.density;
        double u = parameters.velosity;
        parameters.pressure = (ADIABATIC_EXPONENT - 1)*(energy - rho*pow(u, 2)/2);
    }
}