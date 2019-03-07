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

        SimulatedAnnealing(Score::ScoringFunction* scorFunc_, Optimizer::Optimizer* optimizer_, unsigned int seed_,
                           SimulatedAnnealing::Parameters params_);

        bool search() final;

        arma::mat getResultMatrix() final;

        ~SimulatedAnnealing() final = default;

        struct Parameters {
            unsigned int maxIterations = 100000;
            double initTemp = 10000.0;
            unsigned int initialNoTempDropMoves = 1000;
            unsigned int moveCtrlSweep = 100;
            double tolerance = 1e-3;
            unsigned int maxToleranceSweep = 3;
            double maxMoveSize = 5.0;
            double initMoveSize = 0.3;
            double gain = 0.3;
        };

    private:
        Score::ScoringFunction* scorFunc;
        Optimizer::Optimizer* optimizer;

        std::mt19937 rndGenerator;

        arma::mat result;

        Parameters params;
    };

}
#endif //SMOLDOCK_SIMULATEDANNEALING_H
