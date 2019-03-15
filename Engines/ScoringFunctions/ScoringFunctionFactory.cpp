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
                                                            double differential_upsilon,
                                                            bool useNonDefaultCoefficients) {
        if (t == ScoringFunctionType::VinaRigid) {
            return std::make_shared<VinaLikeRigid>(conformer, protein, transform, differential_upsilon);
        } else if (t == ScoringFunctionType::Vina) {
            return std::make_shared<VinaLike>(conformer, protein, transform, differential_upsilon, useNonDefaultCoefficients);
        } else if (t == ScoringFunctionType::VinaCovalentReversible) {
            return std::make_shared<VinaLikeCovalentReversible>(conformer, protein, transform, differential_upsilon, useNonDefaultCoefficients);
        } else {
            return nullptr;
        }
    }

    std::string scoringFunctionTypeToString(ScoringFunctionType t) {
        if (t == ScoringFunctionType::VinaRigid) {
            return "VinaRigid";
        } else if (t == ScoringFunctionType::Vina) {
            return "Vina";
        } else if (t == ScoringFunctionType::VinaCovalentReversible) {
            return "VinaCovalentReversible";
        } else {
            return "UNKNOWNTYPE";
        }
    }

}