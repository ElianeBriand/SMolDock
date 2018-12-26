//
// Created by eliane on 26/12/18.
//

#include "ReScorer.h"

namespace SmolDock {
    ReScorer::ReScorer(SmolDock::Protein &prot, SmolDock::Molecule &mol,
                       std::function<double(const iConformer &, const iTransform &, const iProtein &)>& scorFunc) :
            protein(prot), molecule(mol),
            scoringFunction(scorFunc){

    }

    bool ReScorer::prepare() {
        this->iprotein = this->protein.getiProtein();
        this->iconformer = this->molecule.getInitialConformer();

        return true;
    }

    double ReScorer::getScore() {
        return scoringFunction(this->iconformer, iTransformIdentityInit(), this->iprotein);
    }



}
