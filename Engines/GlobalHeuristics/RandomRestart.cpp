//
// Created by eliane on 08/01/19.
//

#include "RandomRestart.h"


#undef BOOST_LOG

#include <boost/log/trivial.hpp>


namespace SmolDock::Heuristics {


    RandomRestart::RandomRestart(Score::ScoringFunction* scorFunc_, Optimizer::Optimizer* optimizer_,
                                 unsigned int seed_) :
            scorFunc(scorFunc_), optimizer(optimizer_), rnd_generator(seed_) {
        BOOST_LOG_TRIVIAL(debug) << "RandomRestart:seed " << seed_;
        std::uniform_real_distribution<double> dis_real_position(-100.0, 100.0);
        BOOST_LOG_TRIVIAL(debug) << "RandomRestart:GenNbr :  " << dis_real_position(this->rnd_generator);


    }

    bool RandomRestart::search() {

        arma::mat startingCondition = scorFunc->getStartingConditions();


        double score_;
        unsigned int iteration_count = 0;
        while (true) {
            score_ = scorFunc->Evaluate(startingCondition);
            if (score_ == 0) {
                std::uniform_real_distribution<double> dis_real_position(-100.0, 100.0);
                for (unsigned int i = 0; i < startingCondition.n_rows; i++) {
                    startingCondition[i] = dis_real_position(this->rnd_generator);
                }
                iteration_count++;
                continue;
            } else {
                optimizer->optimize(startingCondition);
                this->result = optimizer->getRawResultMatrix();
                break;
            }
        }

        BOOST_LOG_TRIVIAL(debug) << "RandomRestart: " << iteration_count << " restarts";


        return false;
    }

    arma::mat RandomRestart::getResultMatrix() {
        return this->result;
    }
}