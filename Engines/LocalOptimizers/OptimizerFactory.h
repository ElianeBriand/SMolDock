//
// Created by eliane on 09/01/19.
//

#ifndef SMOLDOCK_OPTIMIZERFACTORY_H
#define SMOLDOCK_OPTIMIZERFACTORY_H

#include <memory>

#include <Engines/ScoringFunctions/ScoringFunctionInterface.h>
#include "OptimizerInterface.h"

namespace SmolDock::Optimizer {

    enum class LocalOptimizerType {
        L_BFGS,
        GradientDescentLineSearch
    };

    std::shared_ptr<Optimizer>
    optimizerFactory(LocalOptimizerType t, Score::ScoringFunction* scorFunc, double differential_epsilon);


}

#endif //SMOLDOCK_OPTIMIZERFACTORY_H
