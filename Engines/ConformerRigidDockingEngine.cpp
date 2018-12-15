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
 * Foobar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <memory>
#include <chrono>
#include <thread>



#include "ConformerRigidDockingEngine.h"
#include "Internals/iConformer.h"
#include <Engines/ScoringFunctions/BasicScoringFunction.h>

#include "Utilities/TimingsLog.h"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionRunner.h>
#include <GraphMol/ChemReactions/PreprocessRxn.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <DataStructs/ExplicitBitVect.h>

#undef BOOST_LOG
#include <boost/log/trivial.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>



namespace SmolDock {
    namespace Engine {


        ConformerRigidDockingEngine::ConformerRigidDockingEngine(unsigned int conformer_num) {
            this->conformer_num = conformer_num;
            assert(conformer_num != 0);
        }

        bool ConformerRigidDockingEngine::setProtein(Protein *p) {
            this->orig_protein = p;
            return true;
        }

        bool ConformerRigidDockingEngine::setLigand(Molecule *m) {
            this->orig_ligand = m;

            // Getting a handle on the RWMol
            this->rwmol = this->orig_ligand->getInternalRWMol();


            return true;
        }

        bool ConformerRigidDockingEngine::setupDockingEngine() {

            record_timings(begin_setup);

            std::uniform_int_distribution<> dis(0, std::numeric_limits<int>::max());

            record_timings(begin_conformersgen);

            auto seed = dis(this->rnd_generator);
            BOOST_LOG_TRIVIAL(debug) << "Seed : " << seed;
            this->orig_ligand->generateConformers(this->viConformers, this->conformer_num, seed);

            record_timings(end_conformersgen);

            this->protein = this->orig_protein->getiProtein();



            record_timings(end_iprotgen);



            record_timings(end_setup);


#ifdef SMOLDOCK_VERBOSE_DEBUG
            auto total_setup_duration = static_cast< std::chrono::duration<double> >(end_setup - begin_setup).count();
            auto duration_conformergen = static_cast< std::chrono::duration<double> >(end_conformersgen - begin_conformersgen).count();
            auto duration_iprot = static_cast< std::chrono::duration<double> >(end_iprotgen - end_conformersgen).count();

            BOOST_LOG_TRIVIAL(debug) << "Timings ConformerRigidDockingEngine setup"
                                <<"\n      TOTAL: "
                                  << total_setup_duration << "s"
                                  << "\n      Conformer: "
                                  << duration_conformergen
                                  <<"s [n=" << this->conformer_num << "->" << this->viConformers.size() << "] "
                                  << duration_conformergen/this->conformer_num << "s each"
                                  << "\n      iProt: " << duration_iprot << "s";
#endif
            BOOST_LOG_TRIVIAL(info) << "Conformer docking engine: ready";
            return true;
        }

        void ConformerRigidDockingEngine::runDockingEngine() {
            using namespace boost::accumulators;
            accumulator_set<double, stats<tag::mean, tag::moment<2> > > acc_score;
            accumulator_set<double, stats<tag::mean, tag::moment<2> > > acc_duration;
            std::vector<double> scores;
            scores.reserve(this->conformer_num);

            record_timings(begin_docking);

            for(auto& conformer : this->viConformers)
            {
                record_timings(begin_docking_this_conformer);

                iTransform transform = iTransformZeroInit();
                iGradient translationGradient;
                double score = Score::basic_scoring_func(conformer,transform, this->protein);
                double oldscore;
                double score_progress = 100.0;
                double differential_epsilon = 0.01;
                unsigned int iterationNbr = 0;

                while(score_progress > 0.001 )
                {
                    oldscore = score;

                    // Compute approx gradient to find in which direction to go
                    iTransform transform_dx = transform;
                    transform_dx.transl.x += differential_epsilon;
                    translationGradient.dx = oldscore - Score::basic_scoring_func(conformer,transform_dx, this->protein);

                    iTransform transform_dy = transform;
                    transform_dy.transl.y += differential_epsilon;
                    translationGradient.dy = oldscore -  Score::basic_scoring_func(conformer,transform_dy, this->protein);

                    iTransform transform_dz = transform;
                    transform_dz.transl.z += differential_epsilon;
                    translationGradient.dz = oldscore - Score::basic_scoring_func(conformer,transform_dz, this->protein);


                    // Now we want to know how far in that direction
                    // We do a line-search

                    iTransform transform_linesearch;

                    double step = 0.0001;
                    double previous_step = 0.0;
                    double temp_score;
                    double previous_temp_score = std::numeric_limits<double>::max();

                    // We go further and further in the direction of the gradient, until we are not bettering the score
                    while(true) {
                        transform_linesearch = transform;
                        transform_linesearch.transl.x  += step * translationGradient.dx;
                        transform_linesearch.transl.y  += step * translationGradient.dy;
                        transform_linesearch.transl.z  += step * translationGradient.dz;

                        temp_score = Score::basic_scoring_func(conformer, transform_linesearch, this->protein);

                        double diff_oldscore_tempscore = oldscore - temp_score;
                        double diff_previous_tempscore = previous_temp_score - temp_score;

                        previous_temp_score = temp_score;

                        if (diff_oldscore_tempscore > 0 && // temp_score is lower than the old score thus better
                            diff_previous_tempscore > 0 ) // It is also lower than in the previous iteration
                        {
                            // so we jump further to see if it can be even better
                            previous_step = step;
                            step = step * 10;
                            continue;
                        } else {
                            // We hit a step that is not better than either oldscore or the previous iterator
                            break;
                        }
                    }


                    // We then go back to previous_step for updating the score

                    if(previous_step == 0.0)
                    {
                        // This means we could not do any jump forward
                        // So this is the best we can do (to the precision of the starting step 0.0001 (angstrom))
                        // So we stop the search (this is one of the stop condition, the other is score_progress)
                        break;
                    }

                    // Otherwise we update the transform and go for another round

                    transform.transl.x  += previous_step * translationGradient.dx;
                    transform.transl.y  += previous_step * translationGradient.dy;
                    transform.transl.z  += previous_step * translationGradient.dz;

                    score = Score::basic_scoring_func(conformer,transform, this->protein);
                    score_progress = oldscore - score;

                    iterationNbr++;

                    //BOOST_LOG_TRIVIAL(debug) << "     Score: " << score << " (progress: "<< score_progress <<  " )";

                }

                scores.push_back(score);

                record_timings(end_docking_this_conformer);

                acc_score(score);
                acc_duration(
                        static_cast< std::chrono::duration<double> >(end_docking_this_conformer - begin_docking_this_conformer).count()
                        );

            }

            record_timings(end_docking);

#ifdef SMOLDOCK_VERBOSE_DEBUG
            auto total_docking_duration = static_cast< std::chrono::duration<double> >(end_docking - begin_docking).count();

            BOOST_LOG_TRIVIAL(debug) << "Timings ConformerRigidDockingEngine docking runs"
                                     << "\n      TOTAL: " << total_docking_duration << "s"
                                     << "\n      Mean per conformer: " << mean(acc_duration) << "s"
                                     << "\n      StdDev per conformer: " << std::sqrt(moment<2>(acc_duration)) << "s";
#endif

            BOOST_LOG_TRIVIAL(info) << "Results:";
            BOOST_LOG_TRIVIAL(info) << "   Score Mean: " << mean(acc_score);
            BOOST_LOG_TRIVIAL(info) << "   Score StdDev: " << std::sqrt(moment<2>(acc_score));
            BOOST_LOG_TRIVIAL(info) << "   Scores: ";
            for(auto& score: scores)
            {
                BOOST_LOG_TRIVIAL(info) << "      " << score;
            }

        }

        std::shared_ptr<DockingResult> ConformerRigidDockingEngine::getDockingResult() {
            return std::make_shared<DockingResult>();
        }


        bool ConformerRigidDockingEngine::setDockingBox(AbstractDockingEngine::DockingBoxSetting setting) {
            if (setting != DockingBoxSetting::everything) {
                std::cout << "[!] DockingBoxSetting (that is not DockingBoxSetting::everything) is not yet implemented."
                          << std::endl;
                std::cout << "[ ] Running as if DockingBoxSetting::everything was passed" << std::endl;
                return false;
            }
            return true;
        }


        void ConformerRigidDockingEngine::setRandomSeed(int seed) {
            this->random_seed = seed;
            this->rnd_generator.seed(this->random_seed);
        }


    }

}