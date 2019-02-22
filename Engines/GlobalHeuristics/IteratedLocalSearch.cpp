//
// Created by briand on 2/20/19.
//

#include "IteratedLocalSearch.h"


#undef BOOST_LOG

#include <boost/log/trivial.hpp>
#include <Engines/GlobalHeuristics/Utilities/MetropolisHastings.h>


namespace SmolDock::Heuristics {


    IteratedLocalSearch::IteratedLocalSearch(Score::ScoringFunction *scorFunc_, Optimizer::Optimizer *optimizer_,
                                             unsigned int seed_, IteratedLocalSearch::Parameters params_) :
            scorFunc(scorFunc_), optimizer(optimizer_), rndGenerator(seed_), params(params_) {

    }

    bool IteratedLocalSearch::search() {

        arma::mat currentState = scorFunc->getStartingConditions();
        arma::mat bestState = currentState;

        double score_ = scorFunc->Evaluate(currentState);
        double oldScore = std::numeric_limits<double>::max();
        double bestScore = std::numeric_limits<double>::max();

        double temperature_constant = 50.0;
        double temperature = temperature_constant;

        unsigned int restart_count = 0;
        unsigned int metropolis_count = 0;

        unsigned int allowedMetropolisRestartCount = 5;

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
        currentState = optimizer->getRawResultMatrix();

        assert(score_ < oldScore);

        while (MetropolisAccept(oldScore, score_, temperature, this->rndGenerator)) {

            if (score_ < bestScore) {
                bestScore = score_;
                bestState = currentState;
            }

            metropolis_count++;
            oldScore = score_;
            temperature = temperature_constant * std::pow(0.9, 1 + metropolis_count);
            BOOST_LOG_TRIVIAL(debug) << "IteratedLocalSearch: Temperature = " << temperature;

            std::uniform_real_distribution<double> perturbationDistribution(-2.0, +2.0);
            for (unsigned int i = 0; i < currentState.n_rows; i++) {
                currentState[i] += perturbationDistribution(this->rndGenerator);
            }

            optimizer->optimize(currentState);
            score_ = optimizer->getScore();
            currentState = optimizer->getRawResultMatrix();

            if (score_ == 0 && oldScore == 0) {
                if (bestScore == 0 || allowedMetropolisRestartCount <= 0)
                    break;
                currentState = bestState;
                optimizer->optimize(currentState);
                score_ = optimizer->getScore();
                currentState = optimizer->getRawResultMatrix();
                allowedMetropolisRestartCount--;
                continue; // restart from good result
            }


        }

        this->result = bestState;

        BOOST_LOG_TRIVIAL(debug) << "IteratedLocalSearch: Best score = " << bestScore;
        BOOST_LOG_TRIVIAL(debug) << "IteratedLocalSearch: " << restart_count << " restarts";
        BOOST_LOG_TRIVIAL(debug) << "IteratedLocalSearch: " << metropolis_count << " Metropolis loop";

        return true;
    }

    arma::mat IteratedLocalSearch::getResultMatrix() {
        return this->result;
    }
}