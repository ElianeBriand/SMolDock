//
// Created by briand on 2/20/19.
//

#ifndef SMOLDOCK_ITERATEDLOCALSEARCH_H
#define SMOLDOCK_ITERATEDLOCALSEARCH_H

#include <Engines/ScoringFunctions/ScoringFunctionInterface.h>
#include <Engines/LocalOptimizers/OptimizerInterface.h>

#include "HeuristicInterface.h"

namespace SmolDock::Heuristics {


    class IteratedLocalSearch : public GlobalHeuristic {
    public:
        struct Parameters;

        IteratedLocalSearch(Score::ScoringFunction *scorFunc_, Optimizer::Optimizer *optimizer_, unsigned int seed_,
                            IteratedLocalSearch::Parameters params_);

        bool search() final;

        arma::mat getResultMatrix() final;

        ~IteratedLocalSearch() final = default;

        struct Parameters {

            double proteinMaxRadius = 100.0;

        };

    private:
        Score::ScoringFunction *scorFunc;
        Optimizer::Optimizer *optimizer;

        std::mt19937 rndGenerator;

        arma::mat result;

        Parameters params;

    };

}

#endif //SMOLDOCK_ITERATEDLOCALSEARCH_H
