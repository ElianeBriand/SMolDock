//
// Created by eliane on 08/01/19.
//

#ifndef SMOLDOCK_RANDOMRESTART_H
#define SMOLDOCK_RANDOMRESTART_H

#include <Engines/ScoringFunctions/ScoringFunctionInterface.h>
#include <Engines/LocalOptimizers/OptimizerInterface.h>

#include "HeuristicInterface.h"

namespace SmolDock::Heuristics {


    class RandomRestart : public GlobalHeuristic {
    public:

        struct Parameters;

        RandomRestart(Score::ScoringFunction *scorFunc_, Optimizer::Optimizer *optimizer_, unsigned int seed_,
                      RandomRestart::Parameters params_);

        bool search() final;

        arma::mat getResultMatrix() final;

        ~RandomRestart() final = default;


        struct Parameters {

            double proteinMaxRadius = 100.0;

        };

    private:
        Score::ScoringFunction* scorFunc;
        Optimizer::Optimizer* optimizer;

        std::mt19937 rndGenerator;

        arma::mat result;

        Parameters params;

    };

}


#endif //SMOLDOCK_RANDOMRESTART_H
