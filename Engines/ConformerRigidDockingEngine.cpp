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
            assert(conformer_num > 0);

            scores.reserve(this->conformer_num);
            final_iConformer.reserve(this->conformer_num);
        }

        bool ConformerRigidDockingEngine::setProtein(Protein *p) {
            this->orig_protein = p;
            return true;
        }

        bool ConformerRigidDockingEngine::setLigand(Molecule *m) {
            this->orig_ligand = m;

            // Getting a handle on the RWMol
            // this->rwmol = this->orig_ligand->getInternalRWMol();


            return true;
        }

        bool ConformerRigidDockingEngine::setupDockingEngine() {

            record_timings(begin_setup);

            std::uniform_int_distribution<> dis_int(0, std::numeric_limits<int>::max());
            std::uniform_real_distribution<double> dis_real_position(-20.0, 20.0);

            record_timings(begin_conformersgen);

            int seed = dis_int(this->rnd_generator);
            BOOST_LOG_TRIVIAL(debug) << "Seed : " << seed;
            this->orig_ligand->generateConformers(this->viConformers, this->conformer_num, seed);

            record_timings(end_conformersgen);

            this->protein = this->orig_protein->getiProtein();

            /*

            for(iConformer& conformer: this->viConformers)
            {
                iTransform starting_pos_tr = iTransformIdentityInit();
                starting_pos_tr.transl.x = this->protein.center_x + dis_real_position(this->rnd_generator);
                starting_pos_tr.transl.y = this->protein.center_y + dis_real_position(this->rnd_generator);
                starting_pos_tr.transl.z = this->protein.center_z + dis_real_position(this->rnd_generator);
                applyTransformInPlace(conformer,starting_pos_tr);
            }
            */

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
            accumulator_set<double, stats<tag::mean, tag::moment<2> > > acc_iteration;




            record_timings(begin_docking);

            for(auto& conformer : this->viConformers)
            {

                BOOST_LOG_TRIVIAL(debug) << "Initial score : " << Score::vina_like_rigid_inter_scoring_func(conformer,iTransformIdentityInit(),this->protein);
                continue;

                record_timings(begin_docking_this_conformer);

                GradientDescentLineSearch gradOptimizer(Score::vina_like_rigid_inter_scoring_func);
                gradOptimizer.setProtein(&(this->protein));
                gradOptimizer.setStartingConformer(&conformer);

                gradOptimizer.optimize();

                double score = gradOptimizer.getScore();
                iConformer result = gradOptimizer.getFinalConformer();


                this->scores.push_back(score);
                this->final_iConformer.push_back(result);

                record_timings(end_docking_this_conformer);

                acc_iteration(gradOptimizer.getIterationNumber());
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
            BOOST_LOG_TRIVIAL(info) << "   Iteration Mean: " << mean(acc_iteration);
            BOOST_LOG_TRIVIAL(info) << "   Iteration StdDev: " << std::sqrt(moment<2>(acc_iteration));
            BOOST_LOG_TRIVIAL(info) << "   Score Mean: " << mean(acc_score);
            BOOST_LOG_TRIVIAL(info) << "   Score StdDev: " << std::sqrt(moment<2>(acc_score));
            BOOST_LOG_TRIVIAL(info) << "   Scores: ";
            for(auto& score: scores)
            {
                BOOST_LOG_TRIVIAL(info) << "      " << score;
            }

        }

        std::shared_ptr<DockingResult> ConformerRigidDockingEngine::getDockingResult() {
            auto ret = std::make_shared<DockingResult>();
            for(const auto& confr : this->final_iConformer)
            {
                Molecule finalLigand = this->orig_ligand->deepcopy();
                finalLigand.updateAtomPositionsFromiConformer(confr);
                ret->ligandPoses.emplace_back(finalLigand);
            }
            return ret;
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