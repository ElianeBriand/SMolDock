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
#ifndef SMOLDOCK_L_BFGS_H
#define SMOLDOCK_L_BFGS_H

#include <Engines/ScoringFunctions/ScoringFunctionInterface.h>
#include "OptimizerInterface.h"

namespace SmolDock::Optimizer {


    class L_BFGS : public Optimizer {

    public:
        explicit L_BFGS(Score::ScoringFunction* scoringFunc_, unsigned int maxIteration_ = 100);

        bool optimize(arma::mat startingPoint) final;

        arma::mat getRawResultMatrix() final;

        double getScore() final;

        ~L_BFGS() final = default;


    private:
        Score::ScoringFunction* scoringFunction;
        unsigned int maxIteration;


        arma::mat result;
        double score;
    };


}


#endif //SMOLDOCK_L_BFGS_H
