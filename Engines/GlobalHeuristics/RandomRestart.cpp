//
// Created by eliane on 08/01/19.
//

#include "RandomRestart.h"


#undef BOOST_LOG

#include <boost/log/trivial.hpp>


namespace SmolDock::Heuristics {


    RandomRestart::RandomRestart(Score::ScoringFunction* scorFunc_, Optimizer::Optimizer* optimizer_,
                                 unsigned int seed_, RandomRestart::Parameters params_) :
            scorFunc(scorFunc_), optimizer(optimizer_), rndGenerator(seed_), params(params_) {

    }

    bool RandomRestart::search() {

        arma::mat currentState = scorFunc->getStartingConditions();
        double score_ = scorFunc->Evaluate(currentState);
        unsigned int restart_count = 0;

        while (score_ == 0) {
            restart_count++;
            std::uniform_real_distribution<double> randomRestartDistribution(-this->params.proteinMaxRadius,
                                                                             this->params.proteinMaxRadius); // TODO replace with protein max radius

            for (unsigned int i = 0; i < currentState.n_rows; i++) {
                currentState[i] = randomRestartDistribution(this->rndGenerator);
            }

            score_ = scorFunc->Evaluate(currentState);
        }

        optimizer->optimize(currentState);
        score_ = optimizer->getScore();
        this->result = optimizer->getRawResultMatrix();

        BOOST_LOG_TRIVIAL(debug) << "RandomRestart: " << restart_count << " restarts";


        return true;
    }

    arma::mat RandomRestart::getResultMatrix() {
        return this->result;
    }
}