//
// Created by eliane on 09/02/19.
//

#include "OnlyLocal.h"


namespace SmolDock::Heuristics {

    OnlyLocal::OnlyLocal(Score::ScoringFunction* scorFunc_, Optimizer::Optimizer* optimizer_,
                         unsigned int seed_, OnlyLocal::Parameters params_) :
            scorFunc(scorFunc_), optimizer(optimizer_), rnd_generator(seed_), params(params_) {

    }

    bool OnlyLocal::search() {

        arma::mat startingCondition = scorFunc->getStartingConditions();

        optimizer->optimize(startingCondition);
        this->result = optimizer->getRawResultMatrix();

        return false;
    }

    arma::mat OnlyLocal::getResultMatrix() {
        return this->result;
    }

}