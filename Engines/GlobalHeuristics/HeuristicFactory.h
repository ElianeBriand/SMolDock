//
// Created by eliane on 09/01/19.
//

#ifndef SMOLDOCK_HEURISTICFACTORY_H
#define SMOLDOCK_HEURISTICFACTORY_H


#include <memory>

#include <Engines/ScoringFunctions/ScoringFunctionInterface.h>
#include <Engines/LocalOptimizers/OptimizerInterface.h>
#include "HeuristicInterface.h"

namespace SmolDock::Heuristics {

    enum class GlobalHeuristicType {
        RandomRestart
    };

    std::shared_ptr<GlobalHeuristic> globalHeuristicFactory(GlobalHeuristicType t,
                                                            Score::ScoringFunction* scorFunc,
                                                            Optimizer::Optimizer* localOptimizer,
                                                            unsigned int rng_seed);


}


#endif //SMOLDOCK_HEURISTICFACTORY_H
