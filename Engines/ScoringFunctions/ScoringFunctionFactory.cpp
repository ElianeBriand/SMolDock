//
// Created by eliane on 08/01/19.
//

#include "ScoringFunctionFactory.h"

#include "VinaLikeScoringFunction.h"

namespace SmolDock::Score {

    std::unique_ptr<ScoringFunction> scoringFunctionFactory(ScoringFunctionType t, iConformer conformer,
                                                            iProtein protein,iTransform transform, double differential_upsilon) {
        if (t == ScoringFunctionType::VinaRigid) {
            return std::make_unique<VinaLikeRigidScoringFunction>(conformer,protein,transform,differential_upsilon);
        } else {
            return nullptr;
        }
    }

}