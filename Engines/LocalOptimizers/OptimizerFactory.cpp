//
// Created by eliane on 09/01/19.
//

#include "OptimizerFactory.h"

#include "L_BFGS.h"
#include "GradientDescentLineSearch.h"

namespace SmolDock::Optimizer {

    std::shared_ptr<Optimizer>
    optimizerFactory(LocalOptimizerType t, Score::ScoringFunction* scorFunc, double differential_epsilon) {
        if (t == LocalOptimizerType::L_BFGS) {
            return std::make_shared<L_BFGS>(scorFunc, differential_epsilon);
        } else if (t == LocalOptimizerType::GradientDescentLineSearch) {
            return std::make_shared<GradientDescentLineSearch>(scorFunc, differential_epsilon);
        } else {
            return nullptr;
        }
    }

}