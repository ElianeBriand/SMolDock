/*
 * Copyright (c) 2018 Eliane Briand
 *
 * This file is part of SmolDock.
 *
 * SmolDock is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SmolDock is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SmolDock.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "L_BFGS.h"


#include <ensmallen.hpp>

#undef BOOST_LOG

#include <boost/log/trivial.hpp>


namespace SmolDock::Optimizer {


    L_BFGS::L_BFGS(Score::ScoringFunction* scoringFunc_, unsigned int maxIteration_) :
            scoringFunction(scoringFunc_), maxIteration(maxIteration_),
            result(scoringFunc_->getParamVectorDimension(), 1, arma::fill::zeros),
            score(0.0) {
    }

    bool L_BFGS::optimize(arma::mat startingPoint) {


        ens::L_BFGS lbfgs;
        lbfgs.MaxIterations() = this->maxIteration;
        lbfgs.MinGradientNorm() = 1e-5;
        lbfgs.Factr() = 1e-6;
        lbfgs.MinStep() = 1e-7;

        lbfgs.Optimize(*this->scoringFunction, startingPoint);


        this->result = startingPoint;
        this->score = this->scoringFunction->Evaluate(startingPoint);

        //*
        BOOST_LOG_TRIVIAL(debug) << "LBFGS: x =  " << startingPoint.t();
        BOOST_LOG_TRIVIAL(debug) << "LBFGS: f(x) =  " << this->score;
        //*/

        return true;

    }

    arma::mat L_BFGS::getRawResultMatrix() {
        return this->result;
    }

    double L_BFGS::getScore() {
        return this->score;
    }


};