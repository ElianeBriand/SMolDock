//
// Created by eliane on 09/01/19.
//

#include "HeuristicFactory.h"
#include "RandomRestart.h"

namespace SmolDock::Heuristics {

    std::shared_ptr<GlobalHeuristic> globalHeuristicFactory(GlobalHeuristicType t, Score::ScoringFunction* scorFunc,
                                                            Optimizer::Optimizer* localOptimizer,
                                                            unsigned int rng_seed) {

        if (t == GlobalHeuristicType::RandomRestart) {
            return std::make_shared<RandomRestart>(scorFunc, localOptimizer, rng_seed);
        } else {
            return nullptr;
        }
    }

}