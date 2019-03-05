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
        arma::mat statingState = currentState;
        double score_ = 0.0; //scorFunc->Evaluate(currentState);
        unsigned int restart_count = 0;

        while (score_ == 0 || score_ > 10.0) {
            restart_count++;
            std::uniform_real_distribution<double> randomRestartTranslationDistribution(-this->params.proteinMaxRadius,
                                                                                        this->params.proteinMaxRadius);
            std::uniform_real_distribution<double> randomRestartQuatDistribution(-3, +3);
            std::uniform_real_distribution<double> randomRestartRotDistribution(-180.0 * (3.14 / 180),
                                                                                +180.0 * (3.14 / 180));


            for (unsigned int i = 0; i < 3; i++) {
                double perturbation = randomRestartTranslationDistribution(this->rndGenerator);
                currentState[i] = statingState[i] + perturbation;
            }

            for (unsigned int i = 3; i < 8; i++) {
                double perturbation = randomRestartQuatDistribution(this->rndGenerator);
                currentState[i] = statingState[i] + perturbation;
            }

            for (unsigned int i = 8; i < currentState.n_rows; i++) {
                double perturbation = randomRestartRotDistribution(this->rndGenerator);
                currentState[i] = statingState[i] + perturbation;
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