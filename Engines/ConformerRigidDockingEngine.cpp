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

#include <memory>
#include <chrono>
#include <thread>


#include "ConformerRigidDockingEngine.h"
#include "Internals/iConformer.h"
#include "Internals/InternalsUtilityFunctions.h"
#include <Engines/ScoringFunctions/VinaLikeScoringFunction.h>
#include <Engines/LocalOptimizers/L_BFGS.h>
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
#include <Engines/GlobalHeuristics/RandomRestart.h>


namespace SmolDock {
    namespace Engine {

        ConformerRigidDockingEngine::ConformerRigidDockingEngine(unsigned int conformer_num_,
                                                                 Protein* protein,
                                                                 Molecule* ligand,
                                                                 Score::ScoringFunctionType scFuncType,
                                                                 Heuristics::GlobalHeuristicType heurType,
                                                                 Optimizer::LocalOptimizerType localOptimizerType_,
                                                                 unsigned int seed) :
                conformer_num(conformer_num_),
                orig_protein(protein),
                orig_ligand(ligand),
                scoringFuncType(scFuncType),
                heuristicType(heurType),
                localOptimizerType(localOptimizerType_),
                rnd_generator(seed) {


            scores.reserve(this->conformer_num);
            final_iConformer.reserve(this->conformer_num);
        }


        bool ConformerRigidDockingEngine::setupDockingEngine() {

            record_timings(begin_setup);

            std::uniform_int_distribution<> dis_int(0, std::numeric_limits<int>::max());
            std::uniform_real_distribution<double> dis_real_position(-100.0, 100.0);

            record_timings(begin_conformersgen);

            this->orig_ligand->generateConformers(this->viConformers, this->conformer_num, true,
                                                  dis_int(this->rnd_generator));

            record_timings(end_conformersgen);

            this->protein = this->orig_protein->getiProtein();


/*
            for(iConformer& conformer: this->viConformers)
            {

            }
*/

            record_timings(end_iprotgen);


            record_timings(end_setup);


#ifdef SMOLDOCK_VERBOSE_DEBUG
            auto total_setup_duration = static_cast< std::chrono::duration<double> >(end_setup - begin_setup).count();
            auto duration_conformergen = static_cast< std::chrono::duration<double> >(end_conformersgen -
                                                                                      begin_conformersgen).count();
            auto duration_iprot = static_cast< std::chrono::duration<double> >(end_iprotgen -
                                                                               end_conformersgen).count();

            BOOST_LOG_TRIVIAL(debug) << "Timings ConformerRigidDockingEngine setup"
                                     << "\n      TOTAL: "
                                     << total_setup_duration << "s"
                                     << "\n      Conformer: "
                                     << duration_conformergen
                                     << "s [n=" << this->conformer_num << "->" << this->viConformers.size() << "] "
                                     << duration_conformergen / this->conformer_num << "s each"
                                     << "\n      iProt: " << duration_iprot << "s";
#endif
            BOOST_LOG_TRIVIAL(info) << "Conformer docking engine: ready";
            return true;
        }

        void ConformerRigidDockingEngine::runDockingEngine() {
            using namespace boost::accumulators;
            accumulator_set<double, stats<tag::mean, tag::moment<2> > > acc_score;
            accumulator_set<double, stats<tag::mean, tag::moment<2> > > acc_duration;

            std::uniform_int_distribution<unsigned int> dis_uint(0, std::numeric_limits<unsigned int>::max());


            record_timings(begin_docking);

            for (auto &conformer : this->viConformers) {


                record_timings(begin_docking_this_conformer);


                iTransform starting_pos_tr = iTransformIdentityInit();
                starting_pos_tr.transl = conformer.centroidNormalizingTransform;

                // Remove this for production
                // This is to test the local optimizer by shaking up the initial position
//                starting_pos_tr.transl.x += -1.0;
//                starting_pos_tr.transl.y += 1.0;
//                starting_pos_tr.rota.s += 1.0;
//                starting_pos_tr.rota.u += 1.0;
//                starting_pos_tr.rota.v += 1.0;
//                starting_pos_tr.rota.t += 1.0;
                normalizeQuaternionInPlace(starting_pos_tr.rota);


                this->scoringFunction = scoringFunctionFactory(this->scoringFuncType,
                                                               conformer,
                                                               this->protein,
                                                               starting_pos_tr,
                                                               1e-3);

                this->localOptimizer = optimizerFactory(this->localOptimizerType,
                                                        this->scoringFunction.get(),
                                                        1e-3);


                this->globalHeuristic = globalHeuristicFactory(heuristicType, this->scoringFunction.get(),
                                                               this->localOptimizer.get(),
                                                               dis_uint(this->rnd_generator));

                this->globalHeuristic->search();

                auto rawResultMatrix = this->globalHeuristic->getResultMatrix();
                double score = this->scoringFunction->Evaluate(rawResultMatrix);
                iConformer result = this->scoringFunction->getConformerForParamMatrix(rawResultMatrix);

                double direct_score = Score::vina_like_rigid_inter_scoring_func(conformer, starting_pos_tr,
                                                                                this->protein);
                BOOST_LOG_TRIVIAL(debug) << "Score after L-BFGS: " << score << " (from: " << direct_score << ")";


                this->scores.push_back(score);
                this->final_iConformer.push_back(result);

                record_timings(end_docking_this_conformer);

                acc_score(score);

#ifdef SMOLDOCK_VERBOSE_DEBUG
                acc_duration(
                        static_cast< std::chrono::duration<double> >(end_docking_this_conformer -
                                                                     begin_docking_this_conformer).count()
                );
#endif

            }

            record_timings(end_docking);

#ifdef SMOLDOCK_VERBOSE_DEBUG
            auto total_docking_duration = static_cast< std::chrono::duration<double> >(end_docking -
                                                                                       begin_docking).count();

            BOOST_LOG_TRIVIAL(debug) << "Timings ConformerRigidDockingEngine docking runs"
                                     << "\n      TOTAL: " << total_docking_duration << "s"
                                     << "\n      Mean per conformer: " << mean(acc_duration) << "s"
                                     << "\n      StdDev per conformer: " << std::sqrt(moment<2>(acc_duration)) << "s";
#endif

            BOOST_LOG_TRIVIAL(info) << "Results:";
            BOOST_LOG_TRIVIAL(info) << "   Score Mean: " << mean(acc_score);
            BOOST_LOG_TRIVIAL(info) << "   Score StdDev: " << std::sqrt(moment<2>(acc_score));
            BOOST_LOG_TRIVIAL(info) << "   Scores: ";
            for (auto &score: scores) {
                BOOST_LOG_TRIVIAL(info) << "      " << score;
            }

        }

        std::shared_ptr<DockingResult> ConformerRigidDockingEngine::getDockingResult() {
            auto ret = std::make_shared<DockingResult>();
            for (const auto &confr : this->final_iConformer) {
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


    }

}