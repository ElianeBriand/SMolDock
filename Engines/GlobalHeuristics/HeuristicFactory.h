//
// Created by eliane on 09/01/19.
//

#ifndef SMOLDOCK_HEURISTICFACTORY_H
#define SMOLDOCK_HEURISTICFACTORY_H


#include <memory>
#include <variant>

#include <Engines/ScoringFunctions/ScoringFunctionInterface.h>
#include <Engines/LocalOptimizers/OptimizerInterface.h>
#include "HeuristicInterface.h"

#include "HeuristicFactory.h"
#include "RandomRestart.h"
#include "OnlyLocal.h"
#include "IteratedLocalSearch.h"

namespace SmolDock::Heuristics {

    enum class GlobalHeuristicType {
        RandomRestart,
        OnlyLocal,
        IteratedLocalSearch
    };

    using HeuristicParameters = std::variant<OnlyLocal::Parameters, RandomRestart::Parameters, IteratedLocalSearch::Parameters>;

    extern HeuristicParameters emptyParameters;

    std::shared_ptr<GlobalHeuristic> globalHeuristicFactory(GlobalHeuristicType t,
                                                            Score::ScoringFunction *scorFunc,
                                                            Optimizer::Optimizer *localOptimizer,
                                                            unsigned int rng_seed,
                                                            HeuristicParameters parameters);

    HeuristicParameters heuristicParametersFactory(GlobalHeuristicType t);

}


#endif //SMOLDOCK_HEURISTICFACTORY_H
