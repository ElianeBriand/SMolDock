//
// Created by eliane on 31/12/18.
//

#ifndef SMOLDOCK_SCORINGFUNCTIONINTERFACE_H
#define SMOLDOCK_SCORINGFUNCTIONINTERFACE_H

#include <armadillo>

#include <Engines/Internals/iConformer.h>

namespace SmolDock::Score {




    //! The inferface specification for ScoringFunction object
    class ScoringFunction {
    public:

        virtual double Evaluate(const arma::mat &x) = 0;

        virtual double EvaluateWithGradient(const arma::mat &x, arma::mat &gradient) = 0;


        virtual double getDifferentialEpsilon() const = 0;

        virtual arma::mat getStartingConditions() const = 0;

        virtual unsigned int getParamVectorDimension() const = 0;

        // TODO : Something like virtual checkMatDimension(unsigned int dim)

        virtual SmolDock::iConformer getConformerForParamMatrix(const arma::mat &x) = 0;


        virtual ~ScoringFunction() = default;
    private:


    };







}

#endif //SMOLDOCK_SCORINGFUNCTIONINTERFACE_H
