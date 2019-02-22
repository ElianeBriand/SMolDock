//
// Created by eliane on 23/12/18.
//

#ifndef SMOLDOCK_GRADIENTDESCENTLINESEARCH_H
#define SMOLDOCK_GRADIENTDESCENTLINESEARCH_H

#include <functional>

#include <armadillo>

#include "OptimizerInterface.h"
#include <Engines/ScoringFunctions/ScoringFunctionInterface.h>

#include "Engines/Internals/iConformer.h"
#include "Engines/Internals/iTransform.h"
#include "Engines/Internals/iProtein.h"


namespace SmolDock::Optimizer {

    class GradientDescentLineSearch : public Optimizer {
    public:
        explicit GradientDescentLineSearch(Score::ScoringFunction *scoringFunc_, double differentialUpsilon = 1e-3);

        bool optimize(arma::mat startingPoint) final;

        arma::mat getRawResultMatrix() final;

        double getScore() final;


        unsigned int getIterationNumber();

        ~GradientDescentLineSearch() final = default;

    private:
        Score::ScoringFunction *scoringFunction;
        double differential_epsilon;
        unsigned int paramVectorDimension;
        arma::mat result;

        unsigned int iterationNumber;
        double score;


    };

}


#endif //SMOLDOCK_GRADIENTDESCENTLINESEARCH_H
