//
// Created by briand on 2/22/19.
//

#include "SimulatedAnnealing.h"

#include <ensmallen.hpp>

#undef BOOST_LOG

#include <boost/log/trivial.hpp>
#include <Engines/GlobalHeuristics/Utilities/MetropolisHastings.h>


namespace SmolDock::Heuristics {

    SimulatedAnnealing::SimulatedAnnealing(Score::ScoringFunction* scorFunc_, Optimizer::Optimizer* optimizer_,
                                           unsigned int seed_, SimulatedAnnealing::Parameters params_) :
            scorFunc(scorFunc_), optimizer(optimizer_), rndGenerator(seed_), params(params_) {

    }

    bool SimulatedAnnealing::search() {

        arma::mat currentState = scorFunc->getStartingConditions();

        auto coolingType = ens::ExponentialSchedule();

        unsigned int maxIterations = this->params.maxIterations;
        double initTemp = this->params.initTemp;
        unsigned int initialNoTempDropMoves = this->params.initialNoTempDropMoves;
        unsigned int moveCtrlSweep = this->params.moveCtrlSweep;
        double tolerance = this->params.tolerance;
        unsigned int maxToleranceSweep = this->params.maxToleranceSweep;
        double maxMoveSize = this->params.maxMoveSize;
        double initMoveSize = this->params.initMoveSize;
        double gain = this->params.gain;

        ens::SA<> SimulatedAnnealer(coolingType,
                                    maxIterations,
                                    initTemp ,
                                    initialNoTempDropMoves,
                                    moveCtrlSweep,
                                    tolerance,
                                    maxToleranceSweep,
                                    maxMoveSize,
                                    initMoveSize,
                                    gain);

        SimulatedAnnealer.Optimize(*this->scorFunc, currentState);

        //double scoreBeforeLocalOptim = this->scorFunc->Evaluate(currentState);

        this->optimizer->optimize(currentState);
        //double scoreAfterLocalOptim = this->optimizer->getScore();

        this->result = this->optimizer->getRawResultMatrix();


//        BOOST_LOG_TRIVIAL(debug) << "SimulatedAnnealing: Score after SA          = " << scoreBeforeLocalOptim;
//        BOOST_LOG_TRIVIAL(debug) << "SimulatedAnnealing: Score after local optim = " << scoreAfterLocalOptim;

        return true;
    }

    arma::mat SimulatedAnnealing::getResultMatrix() {
        return this->result;
    }

}