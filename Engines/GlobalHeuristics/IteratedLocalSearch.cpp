//
// Created by briand on 2/20/19.
//

#include "IteratedLocalSearch.h"


#undef BOOST_LOG

#include <boost/log/trivial.hpp>
#include <Engines/GlobalHeuristics/Utilities/MetropolisHastings.h>


namespace SmolDock::Heuristics {


    IteratedLocalSearch::IteratedLocalSearch(Score::ScoringFunction* scorFunc_, Optimizer::Optimizer* optimizer_,
                                             unsigned int seed_, IteratedLocalSearch::Parameters params_) :
            scorFunc(scorFunc_), optimizer(optimizer_), rndGenerator(seed_), params(params_) {

    }

    bool IteratedLocalSearch::search() {

        arma::mat currentState = scorFunc->getStartingConditions();
        arma::mat statingState = currentState;
        arma::mat bestState = currentState;

        double score_ = 0.0; //scorFunc->Evaluate(currentState);
        double oldScore = std::numeric_limits<double>::max();
        double bestScore = std::numeric_limits<double>::max();

        double temperature_constant = 10.0;
        double temperature = temperature_constant;
        double best_temperature = 10.0;

        unsigned int restart_count = 0;
        unsigned int metropolis_count = 0;

        unsigned int allowedMetropolisRestartCount = 3;

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
        currentState = optimizer->getRawResultMatrix();

        assert(score_ < oldScore);

        while (true) {

            if (score_ < bestScore) {
                bestScore = score_;
                bestState = currentState;
                best_temperature = temperature;
            }

            metropolis_count++;
            oldScore = score_;
            temperature = temperature_constant * std::pow(0.9, 1 + metropolis_count);
            BOOST_LOG_TRIVIAL(debug) << "IteratedLocalSearch: Temperature    = " << temperature;
            BOOST_LOG_TRIVIAL(debug) << "IteratedLocalSearch: RestartAllowed = " << allowedMetropolisRestartCount;
            BOOST_LOG_TRIVIAL(debug) << "IteratedLocalSearch: BestScore      = " << bestScore;
            BOOST_LOG_TRIVIAL(debug) << "IteratedLocalSearch: CurrentScore   = " << score_;

            std::uniform_real_distribution<double> perturbationTranslationDistribution(-4.0, +4.0);
            std::uniform_real_distribution<double> perturbationQuatDistribution(-2.0, +2.0);
            std::uniform_real_distribution<double> perturbationRotationDistribution(-20.0 * (3.14 / 180),
                                                                                    +20.0 * (3.14 / 180));
            for (unsigned int i = 0; i < 3; i++) {
                double perturbation = perturbationTranslationDistribution(this->rndGenerator);
                currentState[i] += perturbation;
            }

            for (unsigned int i = 3; i < 8; i++) {
                double perturbation = perturbationQuatDistribution(this->rndGenerator);
                currentState[i] += perturbation;
            }

            for (unsigned int i = 8; i < currentState.n_rows; i++) {
                double perturbation = perturbationRotationDistribution(this->rndGenerator);
                currentState[i] += perturbation;
            }


            optimizer->optimize(currentState);
            score_ = optimizer->getScore();
            currentState = optimizer->getRawResultMatrix();

            bool metropolisAccept = MetropolisAccept(oldScore, score_, temperature, this->rndGenerator);

            if ((score_ == 0 && oldScore == 0) || !metropolisAccept) {
                if (bestScore == 0 || allowedMetropolisRestartCount <= 0)
                    break;
                currentState = bestState;
                allowedMetropolisRestartCount--;
                metropolis_count = 0;
                temperature = best_temperature;
                continue; // restart from good result
            }

        }

        this->result = bestState;

        BOOST_LOG_TRIVIAL(debug) << "IteratedLocalSearch: Best score = " << bestScore;
        BOOST_LOG_TRIVIAL(debug) << "IteratedLocalSearch: " << restart_count << " restarts";
        BOOST_LOG_TRIVIAL(debug) << "======================= END ILS ================================";

        return true;
    }

    arma::mat IteratedLocalSearch::getResultMatrix() {
        return this->result;
    }
}