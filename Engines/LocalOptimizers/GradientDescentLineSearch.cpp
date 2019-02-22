//
// Created by eliane on 23/12/18.
//

#include "GradientDescentLineSearch.h"


#include <Engines/ScoringFunctions/VinaLikeScoringFunction.h>

#include <Engines/Internals/InternalsUtilityFunctions.h>

#include "LocalOptimizerUtilityFunctions.h"

#undef BOOST_LOG

#include <boost/log/trivial.hpp>


namespace SmolDock::Optimizer {


    GradientDescentLineSearch::GradientDescentLineSearch(Score::ScoringFunction *scoringFunc_,
                                                         double differentialUpsilon) :
            scoringFunction(scoringFunc_),
            differential_epsilon(differentialUpsilon),
            paramVectorDimension(scoringFunc_->getParamVectorDimension()),
            result(scoringFunc_->getParamVectorDimension(), 1, arma::fill::zeros),
            iterationNumber(0) {

    }

    bool GradientDescentLineSearch::optimize(arma::mat startingPoint) {

        arma::mat &paramsToOptimize = startingPoint;
        arma::mat gradient(paramVectorDimension, 1, arma::fill::zeros);

        double score_ = this->scoringFunction->EvaluateWithGradient(paramsToOptimize, gradient);

        BOOST_LOG_TRIVIAL(debug) << "Initial score : " << score_;


        double oldscore;
        double score_progress = 100.0;


        while (score_progress > 0.001) {
            if (!hasNonZeroComponent(gradient)) {
                // Score and all gradient components are zero : cannot decide on a direction,
                // --> Search failure
                BOOST_LOG_TRIVIAL(debug)
                    << "Gradient components are zero : cannot use gradient descent for this starting position.";
                return false;
            }


            // We do a linesearch
            // Aka : we want to know how far in the gradient direction we will go

            double step = this->differential_epsilon / 10;
            double previous_step = 0.0;
            double previous_tentativeNewScore = std::numeric_limits<double>::max();

            // We go further and further in the direction of the gradient, until we are not bettering the score
            while (true) {

                arma::mat tentativeNewParams = paramsToOptimize + (step * gradient);

                double tentativeNewScore = this->scoringFunction->Evaluate(tentativeNewParams);

                double initial_tentative_scoreDifference = score_ - tentativeNewScore;
                double previous_tentative_scoreDifference = previous_tentativeNewScore - tentativeNewScore;

                // Prepare for next iteration
                previous_tentativeNewScore = tentativeNewScore;

                if (initial_tentative_scoreDifference > 0
                    && previous_tentative_scoreDifference > 0) {
                    // We are both better than the initial score, and the previous try for step
                    // We try a somewhat larger step
                    previous_step = step;
                    step = step + (this->differential_epsilon / 10);
                } else {
                    // We are not better on both count
                    break; // We stop the linesearch
                }
            }

            // Analysing the result of the linesearch
            if (previous_step == 0.0) {
                // We have not been able to do any step which improves the score
                // Most probably, it means we are at a (local) minimum
                break; // Stopping the gradient search
            }

            // Else :
            // The previous step was an improvement, but this one is not
            // We use previous_step

            arma::mat newParams = paramsToOptimize + (previous_step * gradient);
            paramsToOptimize = newParams;

            oldscore = score_;
            score_ = this->scoringFunction->EvaluateWithGradient(paramsToOptimize, gradient);
            this->iterationNumber++;
            score_progress = oldscore - score_;
            BOOST_LOG_TRIVIAL(debug) << "New score: " << score_ << " (progress: " << score_progress << " )";
        }

        BOOST_LOG_TRIVIAL(debug) << "Finished gradient descent";
        this->score = score_;
        result = paramsToOptimize;
        return true;
    }


    unsigned int GradientDescentLineSearch::getIterationNumber() {
        return this->iterationNumber;
    }

    arma::mat GradientDescentLineSearch::getRawResultMatrix() {
        return this->result;
    }

    double GradientDescentLineSearch::getScore() {
        return this->score;
    }

}