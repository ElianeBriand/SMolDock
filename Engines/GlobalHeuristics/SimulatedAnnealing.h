//
// Created by briand on 2/22/19.
//

#ifndef SMOLDOCK_SIMULATEDANNEALING_H
#define SMOLDOCK_SIMULATEDANNEALING_H

#include <Engines/ScoringFunctions/ScoringFunctionInterface.h>
#include <Engines/LocalOptimizers/OptimizerInterface.h>

#include "HeuristicInterface.h"


namespace SmolDock::Heuristics {

    class SimulatedAnnealing : public GlobalHeuristic {
    public:
        struct Parameters;

        SimulatedAnnealing(Score::ScoringFunction *scorFunc_, Optimizer::Optimizer *optimizer_, unsigned int seed_,
                           SimulatedAnnealing::Parameters params_);

        bool search() final;

        arma::mat getResultMatrix() final;

        ~SimulatedAnnealing() final = default;

        struct Parameters {

        };

    private:
        Score::ScoringFunction *scorFunc;
        Optimizer::Optimizer *optimizer;

        std::mt19937 rndGenerator;

        arma::mat result;

        Parameters params;
    };

}
#endif //SMOLDOCK_SIMULATEDANNEALING_H
