//
// Created by eliane on 26/12/18.
//

#include "ReScorer.h"

#include <exception>

#include <Engines/Internals/InternalsUtilityFunctions.h>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>

namespace SmolDock {
    ReScorer::ReScorer(SmolDock::Protein &prot, SmolDock::Molecule &mol, Score::ScoringFunctionType scorFuncType) :
            protein(prot), molecule(mol),
            scoringFunctionType(scorFuncType) {


    }

    bool ReScorer::prepare() {


        this->iprotein = this->protein.getiProtein();
        this->iconformer = this->molecule.getInitialConformer();
        this->itransform = iTransformIdentityInit();
        this->scoringFunction = scoringFunctionFactory(this->scoringFunctionType, this->iconformer,
                                                       this->iprotein, this->itransform,
                                                       10e-3);

        return true;
    }

    double ReScorer::getScore() {
        if (this->prepared == false) {
            BOOST_LOG_TRIVIAL(error) << "Rescorer getScore() called before prepare() ";
            std::terminate();
        }
        arma::mat state = this->scoringFunction->getStartingConditions();
        double score = this->scoringFunction->Evaluate(state);
        return score;
    }


}
