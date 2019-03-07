//
// Created by eliane on 05/03/19.
//

#include "Evolution.h"



#include <ensmallen.hpp>

#undef BOOST_LOG

#include <boost/log/trivial.hpp>
#include <Engines/GlobalHeuristics/Utilities/MetropolisHastings.h>


namespace SmolDock::Heuristics {

    Evolution::Evolution(Score::ScoringFunction* scorFunc_, Optimizer::Optimizer* optimizer_,
                                                 unsigned int seed_, Evolution::Parameters params_) :
            scorFunc(scorFunc_), optimizer(optimizer_), rndGenerator(seed_), params(params_) {

    }

    bool Evolution::search() {

        arma::mat currentState = scorFunc->getStartingConditions();


        unsigned int populationSize = 500;
        unsigned int maxGenerations = 5000;
        double mutationProb = 0.1;
        double mutationSize = 0.5;
        double selectPercent = 0.2;
        double tolerance = -1;
        double objectiveChange = 1e-3;

        ens::CNE CNEvolutionHeur(populationSize, maxGenerations, mutationProb, mutationSize, selectPercent, tolerance, objectiveChange);

        CNEvolutionHeur.Optimize(*this->scorFunc, currentState);

        double scoreBeforeLocalOptim = this->scorFunc->Evaluate(currentState);

        this->optimizer->optimize(currentState);
        double scoreAfterLocalOptim = this->optimizer->getScore();

        this->result = this->optimizer->getRawResultMatrix();


        BOOST_LOG_TRIVIAL(debug) << "Evolution: Score after Evolution   = " << scoreBeforeLocalOptim;
        BOOST_LOG_TRIVIAL(debug) << "Evolution: Score after local optim = " << scoreAfterLocalOptim;

        return true;
    }

    arma::mat Evolution::getResultMatrix() {
        return this->result;
    }

}