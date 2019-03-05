//
// Created by eliane on 09/02/19.
//

#ifndef SMOLDOCK_ONLYLOCAL_H
#define SMOLDOCK_ONLYLOCAL_H

#include <Engines/ScoringFunctions/ScoringFunctionInterface.h>
#include <Engines/LocalOptimizers/OptimizerInterface.h>

#include "HeuristicInterface.h"

namespace SmolDock::Heuristics {


    class OnlyLocal : public GlobalHeuristic {
    public:

        struct Parameters;

        OnlyLocal(Score::ScoringFunction* scorFunc_, Optimizer::Optimizer* optimizer_, unsigned int seed_,
                  OnlyLocal::Parameters params_);

        bool search() final;

        arma::mat getResultMatrix() final;

        ~OnlyLocal() final = default;


        struct Parameters {

        };
    private:
        Score::ScoringFunction* scorFunc;
        Optimizer::Optimizer* optimizer;

        std::mt19937 rnd_generator;

        arma::mat result;

        Parameters params;
    };

}

#endif //SMOLDOCK_ONLYLOCAL_H
