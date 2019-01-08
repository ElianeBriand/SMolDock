//
// Created by eliane on 08/01/19.
//

#ifndef SMOLDOCK_RANDOMRESTART_H
#define SMOLDOCK_RANDOMRESTART_H

#include <Engines/ScoringFunctions/ScoringFunctionInterface.h>
#include <Engines/LocalOptimizers/OptimizerInterface.h>

#include "HeuristicInterface.h"

namespace SmolDock::Heuristics {

    class RandomRestart : GlobalHeuristic {
    public:

        RandomRestart(Score::ScoringFunction* scorFunc_, Optimizer::Optimizer* optimizer_, unsigned int seed_);

        bool search() final;

        arma::mat getResultMatrix();

    private:
        Score::ScoringFunction* scorFunc;
        Optimizer::Optimizer* optimizer;

        std::mt19937 rnd_generator;

        arma::mat result;

    };

}


#endif //SMOLDOCK_RANDOMRESTART_H
