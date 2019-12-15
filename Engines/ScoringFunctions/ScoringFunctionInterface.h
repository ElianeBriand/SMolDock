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

        // For use with the optimizer
        virtual double Evaluate(const arma::mat &x) = 0;
        virtual double EvaluateWithGradient(const arma::mat &x, arma::mat &gradient) = 0;

        // For use by the engine, as output of the optimization process
        virtual std::vector<std::tuple<std::string,double>> EvaluateSubcomponents(const arma::mat &x) = 0;
        virtual double EvaluateOnlyIntermolecular(const arma::mat &x) = 0;

        virtual double getDifferentialEpsilon() const = 0;
        virtual arma::mat getStartingConditions() const = 0;
        virtual unsigned int getParamVectorDimension() const = 0;
        virtual SmolDock::iConformer getConformerForParamMatrix(const arma::mat &x) = 0;

        virtual unsigned int getCoefficientsVectorWidth() = 0;
        virtual std::vector<std::string> getCoefficientsNames() = 0;
        virtual std::vector<double> getCurrentCoefficients() = 0;
        virtual bool setNonDefaultCoefficients(std::vector<double> coeffs) = 0;

        virtual ~ScoringFunction() = default;

    private:

    };


}

#endif //SMOLDOCK_SCORINGFUNCTIONINTERFACE_H
