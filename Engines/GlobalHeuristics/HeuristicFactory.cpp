//
// Created by eliane on 09/01/19.
//

#include "HeuristicFactory.h"
#include "RandomRestart.h"
#include "OnlyLocal.h"

namespace SmolDock::Heuristics {

    std::shared_ptr<GlobalHeuristic> globalHeuristicFactory(GlobalHeuristicType t, Score::ScoringFunction* scorFunc,
                                                            Optimizer::Optimizer* localOptimizer,
                                                            unsigned int rng_seed) {

        if (t == GlobalHeuristicType::RandomRestart) {
            return std::make_shared<RandomRestart>(scorFunc, localOptimizer, rng_seed);
        }else if (t == GlobalHeuristicType::OnlyLocal) {
            return std::make_shared<OnlyLocal>(scorFunc, localOptimizer, rng_seed);
        } else {
            return nullptr;
        }
    }

}