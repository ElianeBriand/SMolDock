//
// Created by eliane on 08/01/19.
//

#include "ScoringFunctionFactory.h"

#include "VinaLikeScoringFunction.h"

namespace SmolDock::Score {

    std::shared_ptr<ScoringFunction> scoringFunctionFactory(ScoringFunctionType t, const iConformer &conformer,
                                                            const iProtein &protein, const iTransform &transform,
                                                            double differential_upsilon) {
        if (t == ScoringFunctionType::VinaRigid) {
            return std::make_shared<VinaLikeRigidScoringFunction>(conformer, protein, transform, differential_upsilon);
        } else {
            return nullptr;
        }
    }

}