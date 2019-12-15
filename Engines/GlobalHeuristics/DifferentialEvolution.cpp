//
// Created by eliane on 05/03/19.
//

#include "DifferentialEvolution.h"


#include <ensmallen.hpp>

#undef BOOST_LOG

#include <boost/log/trivial.hpp>
#include <Engines/GlobalHeuristics/Utilities/MetropolisHastings.h>


namespace SmolDock::Heuristics {

    DifferentialEvolution::DifferentialEvolution(Score::ScoringFunction* scorFunc_, Optimizer::Optimizer* optimizer_,
                                           unsigned int seed_, DifferentialEvolution::Parameters params_) :
            scorFunc(scorFunc_), optimizer(optimizer_), rndGenerator(seed_), params(params_) {
                BOOST_LOG_TRIVIAL(error) << "DifferentialEvolution is not implemented. ";


    }

    bool DifferentialEvolution::search() {

        arma::mat currentState = scorFunc->getStartingConditions();


        unsigned int populationSize = 100;
        unsigned int maxGenerations = 25000;
        double crossoverRate = 0.6;
        double differentialWeight = 0.02;

        ens::DE DiffentialEvolutionHeur(populationSize, maxGenerations, crossoverRate, differentialWeight);

        DiffentialEvolutionHeur.Optimize(*this->scorFunc, currentState);

        double scoreBeforeLocalOptim = this->scorFunc->Evaluate(currentState);

        this->optimizer->optimize(currentState);
        double scoreAfterLocalOptim = this->optimizer->getScore();

        this->result = this->optimizer->getRawResultMatrix();


        BOOST_LOG_TRIVIAL(debug) << "DifferentialEvolution: Score after DE          = " << scoreBeforeLocalOptim;
        BOOST_LOG_TRIVIAL(debug) << "DifferentialEvolution: Score after local optim = " << scoreAfterLocalOptim;

        return true;
    }

    arma::mat DifferentialEvolution::getResultMatrix() {
        return this->result;
    }

}