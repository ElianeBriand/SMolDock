//
// Created by eliane on 05/03/19.
//

#ifndef SMOLDOCK_DIFFERENTIALEVOLUTION_H
#define SMOLDOCK_DIFFERENTIALEVOLUTION_H

#include <Engines/ScoringFunctions/ScoringFunctionInterface.h>
#include <Engines/LocalOptimizers/OptimizerInterface.h>

#include "HeuristicInterface.h"



namespace SmolDock::Heuristics {

    class DifferentialEvolution : public GlobalHeuristic {
    public:
        struct Parameters;

        DifferentialEvolution(Score::ScoringFunction* scorFunc_, Optimizer::Optimizer* optimizer_, unsigned int seed_,
                              DifferentialEvolution::Parameters params_);

        bool search() final;

        arma::mat getResultMatrix() final;

        ~DifferentialEvolution() final = default;

        struct Parameters {

        };

    private:
        Score::ScoringFunction* scorFunc;
        Optimizer::Optimizer* optimizer;

        std::mt19937 rndGenerator;

        arma::mat result;

        Parameters params;
    };

}




#endif //SMOLDOCK_DIFFERENTIALEVOLUTION_H
