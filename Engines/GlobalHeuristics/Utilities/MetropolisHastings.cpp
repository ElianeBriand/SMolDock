//
// Created by briand on 2/20/19.
//

#include "MetropolisHastings.h"


// Highly inspired of autodock Vina Metropolis.cpp
bool
SmolDock::Heuristics::MetropolisAccept(double oldScore, double newScore, double temperature,
                                       std::mt19937 &rndGenerator) {

    if (newScore < oldScore) {
        return true; // Unconditional accept if better
    }

    const double acceptanceProbability = std::exp((oldScore - newScore) / temperature);

    std::uniform_real_distribution<double> realProbaDistribution(0.0, 1.0);
    const double value = realProbaDistribution(rndGenerator);

    return value < acceptanceProbability;
}
