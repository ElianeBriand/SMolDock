//
// Created by eliane on 09/01/19.
//

#include "HeuristicFactory.h"
#include "RandomRestart.h"
#include "OnlyLocal.h"
#include "IteratedLocalSearch.h"
#include "SimulatedAnnealing.h"

namespace SmolDock::Heuristics {

    HeuristicParameters emptyParameters;

    std::shared_ptr<GlobalHeuristic> globalHeuristicFactory(GlobalHeuristicType t, Score::ScoringFunction* scorFunc,
                                                            Optimizer::Optimizer* localOptimizer,
                                                            unsigned int rng_seed,
                                                            HeuristicParameters parameters) {

        if (t == GlobalHeuristicType::RandomRestart) {
            return std::make_shared<RandomRestart>(scorFunc, localOptimizer, rng_seed,
                                                   std::get<RandomRestart::Parameters>(parameters));
        } else if (t == GlobalHeuristicType::OnlyLocal) {
            return std::make_shared<OnlyLocal>(scorFunc, localOptimizer, rng_seed,
                                               std::get<OnlyLocal::Parameters>(parameters));
        } else if (t == GlobalHeuristicType::IteratedLocalSearch) {
            return std::make_shared<IteratedLocalSearch>(scorFunc, localOptimizer, rng_seed,
                                                         std::get<IteratedLocalSearch::Parameters>(parameters));
        } else if (t == GlobalHeuristicType::SimulatedAnnealing) {
            return std::make_shared<SimulatedAnnealing>(scorFunc, localOptimizer, rng_seed,
                                                        std::get<SimulatedAnnealing::Parameters>(parameters));
        } else {
            return nullptr;
        }
    }

    HeuristicParameters heuristicParametersFactory(GlobalHeuristicType t) {
        if (t == GlobalHeuristicType::RandomRestart) {
            return HeuristicParameters(std::in_place_type_t<RandomRestart::Parameters>());
        } else if (t == GlobalHeuristicType::OnlyLocal) {
            return HeuristicParameters(std::in_place_type_t<OnlyLocal::Parameters>());
        } else if (t == GlobalHeuristicType::IteratedLocalSearch) {
            return HeuristicParameters(std::in_place_type_t<IteratedLocalSearch::Parameters>());
        } else if (t == GlobalHeuristicType::SimulatedAnnealing) {
            return HeuristicParameters(std::in_place_type_t<SimulatedAnnealing::Parameters>());
        } else {
            return emptyParameters;
        }
    }

}