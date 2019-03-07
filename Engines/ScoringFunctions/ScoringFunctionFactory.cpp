//
// Created by eliane on 08/01/19.
//

#include "ScoringFunctionFactory.h"

#include "VinaLikeRigid.h"
#include "VinaLike.h"
#include "VinaLikeCovalentReversible.h"

namespace SmolDock::Score {

    std::shared_ptr<ScoringFunction> scoringFunctionFactory(ScoringFunctionType t, const iConformer &conformer,
                                                            const iProtein &protein, const iTransform &transform,
                                                            double differential_upsilon) {
        if (t == ScoringFunctionType::VinaRigid) {
            return std::make_shared<VinaLikeRigid>(conformer, protein, transform, differential_upsilon);
        } else if (t == ScoringFunctionType::Vina) {
            return std::make_shared<VinaLike>(conformer, protein, transform, differential_upsilon);
        } else if (t == ScoringFunctionType::VinaCovalentReversible) {
            return std::make_shared<VinaLikeCovalentReversible>(conformer, protein, transform, differential_upsilon);
        } else {
            return nullptr;
        }
    }

}