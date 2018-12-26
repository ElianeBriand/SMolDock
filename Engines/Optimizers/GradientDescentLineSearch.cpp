//
// Created by eliane on 23/12/18.
//

#include "GradientDescentLineSearch.h"


#include <Utilities/IntermediateConformerCollector.h>

extern SmolDock::IntermediateConformerCollector *conformerCollector;

#undef BOOST_LOG

#include <boost/log/trivial.hpp>


namespace SmolDock {


    GradientDescentLineSearch::GradientDescentLineSearch(
            std::function<double(const iConformer &, const iTransform &, const iProtein &)> scorFunc,
            double differentialUpsilon) :
            scoringFunction(scorFunc),
            differential_epsilon(differentialUpsilon) {

    }

    bool GradientDescentLineSearch::optimize() {
        this->currentTransform = iTransformIdentityInit();
        iGradient translationGradient;
        double score = this->scoringFunction(this->currentConformer, this->currentTransform, *(this->protein));
        double oldscore;
        double score_progress = 100.0;

        while (score_progress > 0.001) {
            oldscore = score;

            BOOST_LOG_TRIVIAL(info) << "\n";
            BOOST_LOG_TRIVIAL(debug) << "score: " << score;

            // Compute approx gradient to find in which direction to go

            // Translation
            iTransform transform_dx = this->currentTransform;
            transform_dx.transl.x += this->differential_epsilon;
            translationGradient.dx =
                    oldscore - this->scoringFunction(this->currentConformer, transform_dx, *(this->protein));


            iTransform transform_dy = this->currentTransform;
            transform_dy.transl.y += this->differential_epsilon;
            translationGradient.dy =
                    oldscore - this->scoringFunction(this->currentConformer, transform_dy, *(this->protein));

            iTransform transform_dz = this->currentTransform;
            transform_dz.transl.z += this->differential_epsilon;
            translationGradient.dz =
                    oldscore - this->scoringFunction(this->currentConformer, transform_dz, *(this->protein));

            // Rotation
            iTransform transform_ds = this->currentTransform;
            transform_ds.rota.s += this->differential_epsilon;
            normalizeQuaternionInPlace(transform_ds.rota);
            translationGradient.ds =
                    oldscore - this->scoringFunction(this->currentConformer, transform_ds, *(this->protein));

            iTransform transform_du = this->currentTransform;
            transform_du.rota.u += this->differential_epsilon;
            normalizeQuaternionInPlace(transform_du.rota);
            translationGradient.du =
                    oldscore - this->scoringFunction(this->currentConformer, transform_du, *(this->protein));

            iTransform transform_dv = this->currentTransform;
            transform_dv.rota.v += this->differential_epsilon;
            normalizeQuaternionInPlace(transform_dv.rota);
            translationGradient.dv =
                    oldscore - this->scoringFunction(this->currentConformer, transform_dv, *(this->protein));

            iTransform transform_dt = this->currentTransform;
            transform_dt.rota.t += this->differential_epsilon;
            normalizeQuaternionInPlace(transform_dt.rota);
            translationGradient.dt =
                    oldscore - this->scoringFunction(this->currentConformer, transform_dt, *(this->protein));



            BOOST_LOG_TRIVIAL(debug) << "Gradient";
            BOOST_LOG_TRIVIAL(debug) << "     ds: " << translationGradient.ds;
            BOOST_LOG_TRIVIAL(debug) << "     du: " << translationGradient.du << "   dx: " << translationGradient.dx;
            BOOST_LOG_TRIVIAL(debug) << "     dv: " << translationGradient.dv << "   dy: " << translationGradient.dy;
            BOOST_LOG_TRIVIAL(debug) << "     dt: " << translationGradient.dt << "   dx: " << translationGradient.dz;
            // Now we want to know how far in that direction
            // We do a line-search

            iTransform transform_linesearch;

            double step = this->differential_epsilon / 10;
            double previous_step = 0.0;
            double temp_score;
            double previous_temp_score = std::numeric_limits<double>::max();

            // We go further and further in the direction of the gradient, until we are not bettering the score
            while (true) {
                transform_linesearch = this->currentTransform;


                transform_linesearch.transl.x += step * translationGradient.dx;
                transform_linesearch.transl.y += step * translationGradient.dy;
                transform_linesearch.transl.z += step * translationGradient.dz;

                transform_linesearch.rota.s += step * translationGradient.ds;
                transform_linesearch.rota.u += step * translationGradient.du;
                transform_linesearch.rota.v += step * translationGradient.dv;
                transform_linesearch.rota.t += step * translationGradient.dt;
                normalizeQuaternionInPlace(transform_linesearch.rota);

                temp_score = this->scoringFunction(this->currentConformer, transform_linesearch, *(this->protein));


                double diff_oldscore_tempscore = oldscore - temp_score;
                double diff_previous_tempscore = previous_temp_score - temp_score;

                /*
                BOOST_LOG_TRIVIAL(info) << "";
                BOOST_LOG_TRIVIAL(debug) << "oldscore: " << oldscore;
                BOOST_LOG_TRIVIAL(debug) << "temp_score: " << temp_score;
                BOOST_LOG_TRIVIAL(debug) << "previous_temp_score: " << previous_temp_score;
                BOOST_LOG_TRIVIAL(debug) << "diff_oldscore_tempscore: " << diff_oldscore_tempscore;
                BOOST_LOG_TRIVIAL(debug) << "diff_previous_tempscore: " << diff_previous_tempscore;
*/

                previous_temp_score = temp_score;

                if (diff_oldscore_tempscore > 0  // temp_score is lower than the old score thus better
                     && diff_previous_tempscore > 0 ) // It is also lower than in the previous iteration
                {

                    // so we jump further to see if it can be even better
                    previous_step = step;
                    step = step * 10;
                    BOOST_LOG_TRIVIAL(debug) << "Going further: " << step;
                    continue;
                } else {
                    BOOST_LOG_TRIVIAL(debug) << "Regression, stepping back to " << previous_step;

                    // We hit a step that is not better than either oldscore or the previous iterator
                    break;
                }
            }


            // We then go back to previous_step for updating the score

            if (previous_step == 0.0) {
                // This means we could not do any jump forward
                // So this is the best we can do (to the precision of the starting step 0.0001 (angstrom))
                // So we stop the search (this is one of the stop condition, the other is score_progress)
                break;
            }

            // Otherwise we update the transform and go for another round


            this->currentTransform.transl.x += previous_step * translationGradient.dx;
            this->currentTransform.transl.y += previous_step * translationGradient.dy;
            this->currentTransform.transl.z += previous_step * translationGradient.dz;

            this->currentTransform.rota.s += previous_step * translationGradient.ds;
            this->currentTransform.rota.u += previous_step * translationGradient.du;
            this->currentTransform.rota.v += previous_step * translationGradient.dv;
            this->currentTransform.rota.t += previous_step * translationGradient.dt;
            normalizeQuaternionInPlace(this->currentTransform.rota);


            BOOST_LOG_TRIVIAL(debug) << "== Accepted transform after iteration ==";
            BOOST_LOG_TRIVIAL(debug) << "     s: " << this->currentTransform.rota.s;
            BOOST_LOG_TRIVIAL(debug) << "     u: " << this->currentTransform.rota.u << "   x: " << this->currentTransform.transl.x;
            BOOST_LOG_TRIVIAL(debug) << "     v: " << this->currentTransform.rota.v << "   y: " << this->currentTransform.transl.y;
            BOOST_LOG_TRIVIAL(debug) << "     t: " << this->currentTransform.rota.t << "   x: " << this->currentTransform.transl.z;

            score = this->scoringFunction(this->currentConformer, this->currentTransform, *(this->protein));
            score_progress = oldscore - score;

            this->iterationNbr++;

            conformerCollector->addiConformer(this->getCurrentConformerWithTransform());

            //BOOST_LOG_TRIVIAL(debug) << "     Score: " << score << " (progress: "<< score_progress <<  " )";

        }

        this->score = score;

        /*
        scores.push_back(score);
        final_transform.push_back(transform);

        record_timings(end_docking_this_conformer);

        acc_iteration(iterationNbr);
        acc_score(score);
        acc_duration(
                static_cast< std::chrono::duration<double> >(end_docking_this_conformer -
                                                             begin_docking_this_conformer).count()
        );
         */


        return true;
    }

    void GradientDescentLineSearch::setProtein(const iProtein *prot) {
        this->protein = prot;
    }

    void GradientDescentLineSearch::setStartingConformer(const iConformer *conformer) {
        this->startingConformer = conformer;
        conformerCollector->addiConformer(*this->startingConformer);
        this->currentConformer = *conformer;
    }

    iConformer GradientDescentLineSearch::getFinalConformer() {
        applyTransformInPlace(this->currentConformer, this->currentTransform);
        return this->currentConformer;
    }

    double GradientDescentLineSearch::getScore() {
        return this->score;
    }

    iConformer GradientDescentLineSearch::getCurrentConformerWithTransform() {
        iConformer current = this->currentConformer;
        applyTransformInPlace(current, this->currentTransform);
        return current;
    }

    unsigned int GradientDescentLineSearch::getIterationNumber() {
        return this->iterationNbr;
    }

}