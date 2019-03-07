//
// Created by eliane on 08/01/19.
//

#ifndef SMOLDOCK_SCORINGFUNCTIONFACTORY_H
#define SMOLDOCK_SCORINGFUNCTIONFACTORY_H

#include <memory>

#include <Engines/Internals/iProtein.h>

#include "ScoringFunctionInterface.h"

namespace SmolDock::Score {


    enum class ScoringFunctionType {
        VinaRigid,
        Vina,
        VinaCovalentReversible
    };

    std::shared_ptr<ScoringFunction> scoringFunctionFactory(ScoringFunctionType t, const iConformer &conformer,
                                                            const iProtein &protein, const iTransform &transform,
                                                            double differential_upsilon);


}

#endif //SMOLDOCK_SCORINGFUNCTIONFACTORY_H
