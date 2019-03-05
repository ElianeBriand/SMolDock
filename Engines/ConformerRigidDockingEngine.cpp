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
#include <algorithm>

#include "ConformerRigidDockingEngine.h"
#include "Internals/InternalsUtilityFunctions.h"
#include <Engines/ScoringFunctions/VinaLikeRigidScoringFunction.h>
#include <Engines/LocalOptimizers/L_BFGS.h>
#include "Utilities/TimingsLog.h"

#undef BOOST_LOG

#include <boost/log/trivial.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>


namespace SmolDock {
    namespace Engine {

        ConformerRigidDockingEngine::ConformerRigidDockingEngine(unsigned int conformer_num_,
                                                                 unsigned int retryPerConformer,
                                                                 Protein* protein,
                                                                 Molecule* ligand,
                                                                 Score::ScoringFunctionType scFuncType,
                                                                 Heuristics::GlobalHeuristicType heurType,
                                                                 Optimizer::LocalOptimizerType localOptimizerType_,
                                                                 unsigned int seed) :
                conformer_num(conformer_num_),
                retryPerConformer(retryPerConformer),
                orig_protein(protein),
                orig_ligand(ligand),
                scoringFuncType(scFuncType),
                heuristicType(heurType),
                localOptimizerType(localOptimizerType_),
                rnd_generator(seed) {


            scores.reserve(this->conformer_num);
            allGeneratediConformer.reserve(this->conformer_num * this->retryPerConformer);
        }


        bool ConformerRigidDockingEngine::setupDockingEngine() {

            record_timings(begin_setup);

            std::uniform_int_distribution<> dis_int(0, std::numeric_limits<int>::max());
            std::uniform_real_distribution<double> dis_real_position(-100.0, 100.0);

            record_timings(begin_conformersgen);

            this->orig_ligand->generateConformers(this->viConformers, this->conformer_num, true,
                                                  dis_int(this->rnd_generator));

            record_timings(end_conformersgen);

            if (this->dockBoxSettings.type == DockingBoxSetting::Type::centeredAround) {
                this->protein = this->orig_protein->getPartialiProtein_sphere(this->dockBoxSettings.center,
                                                                              this->dockBoxSettings.radius, 2.0);
            } else {
                this->protein = this->orig_protein->getiProtein();
            }

            this->fullProtein = this->orig_protein->getiProtein();


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


                iTransform starting_pos_tr = iTransformIdentityInit(conformer.num_rotatable_bond);
                starting_pos_tr.transl = conformer.centroidNormalizingTransform;


                DockingBoxSetting dbsettings = this->dockBoxSettings;
                unsigned int seed = dis_uint(this->rnd_generator);

                arma::arma_rng::set_seed(seed);

                if (dbsettings.type == DockingBoxSetting::Type::centeredAround) {
                    starting_pos_tr.transl.x() += dbsettings.center[0];
                    starting_pos_tr.transl.y() += dbsettings.center[1];
                    starting_pos_tr.transl.z() += dbsettings.center[2];
                }

                starting_pos_tr.rota.normalize();

                std::shared_ptr<Score::ScoringFunction> scoringFunction = scoringFunctionFactory(this->scoringFuncType,
                                                                                                 conformer,
                                                                                                 this->protein,
                                                                                                 starting_pos_tr,
                                                                                                 1e-3);

                std::shared_ptr<Optimizer::Optimizer> localOptimizer = optimizerFactory(this->localOptimizerType,
                                                                                        scoringFunction.get(),
                                                                                        1e-3);

                Heuristics::HeuristicParameters hParams = Heuristics::heuristicParametersFactory(heuristicType);

                if (heuristicType == Heuristics::GlobalHeuristicType::IteratedLocalSearch) {

                    if (dbsettings.type == DockingBoxSetting::Type::centeredAround) {
                        std::get<Heuristics::IteratedLocalSearch::Parameters>(
                                hParams).proteinMaxRadius = dbsettings.radius;
                    } else {
                        double protRadius = this->orig_protein->getMaxRadius();
                        std::get<Heuristics::IteratedLocalSearch::Parameters>(hParams).proteinMaxRadius = protRadius;
                    }
                }

                if (heuristicType == Heuristics::GlobalHeuristicType::RandomRestart) {
                    if (dbsettings.type == DockingBoxSetting::Type::centeredAround) {
                        std::get<Heuristics::RandomRestart::Parameters>(
                                hParams).proteinMaxRadius = dbsettings.radius;
                    } else {
                        double protRadius = this->orig_protein->getMaxRadius();
                        std::get<Heuristics::RandomRestart::Parameters>(hParams).proteinMaxRadius = protRadius;
                    }
                }


                for (unsigned int retryNum = 0; retryNum < this->retryPerConformer; retryNum++) {

                    std::shared_ptr<Heuristics::GlobalHeuristic> globalHeuristic = globalHeuristicFactory(heuristicType,
                                                                                                          scoringFunction.get(),
                                                                                                          localOptimizer.get(),
                                                                                                          seed +
                                                                                                          retryNum,
                                                                                                          hParams);

                    globalHeuristic->search();

                    auto rawResultMatrix = globalHeuristic->getResultMatrix();
                    double score = scoringFunction->Evaluate(rawResultMatrix);


                    iConformer result = scoringFunction->getConformerForParamMatrix(rawResultMatrix);

                    double starting_score = Score::vina_like_rigid_inter_scoring_func(conformer, starting_pos_tr,
                                                                                      this->protein);

                    double real_score = Score::vina_like_rigid_inter_scoring_func(result, iTransformIdentityInit(),
                                                                                  this->fullProtein);

                    BOOST_LOG_TRIVIAL(debug) << "Starting score : " << starting_score;
                    BOOST_LOG_TRIVIAL(debug) << "Score after L-BFGS: " << score;
                    BOOST_LOG_TRIVIAL(debug) << "Real score (full prot): " << real_score;

                    if (score != 0.0) {
                        this->startingScores.push_back(starting_score);
                        this->localScores.push_back(score);
                        this->scores.push_back(real_score);
                        this->allGeneratediConformer.push_back(result);
                        acc_score(real_score);
                    }

                }


                record_timings(end_docking_this_conformer);


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
            for (unsigned int i = 0; i < this->scores.size(); i++) {
                BOOST_LOG_TRIVIAL(info) << "      " << this->startingScores[i] << " -> " << this->localScores[i]
                                        << " -> " << this->scores[i];
            }

            BOOST_LOG_TRIVIAL(info) << "   Best scores : ";

            std::vector<std::tuple<int, double>> scoreAndIndices;
            for (unsigned int i = 0; i < this->scores.size(); i++) {
                scoreAndIndices.push_back(std::make_tuple(i, this->scores[i]));
            }

            std::nth_element(std::begin(scoreAndIndices),
                             std::begin(scoreAndIndices) + conformer_num,
                             std::end(scoreAndIndices),
                             [](const std::tuple<int, double> &a, const std::tuple<int, double> &b) {
                                 return std::get<1>(a) < std::get<1>(b);
                             });

            for (unsigned int i = 0; i < conformer_num; i++) {
                BOOST_LOG_TRIVIAL(info) << "      " << this->startingScores[std::get<0>(scoreAndIndices[i])]
                                        << " -> " << this->localScores[std::get<0>(scoreAndIndices[i])]
                                        << " -> " << this->scores[std::get<0>(scoreAndIndices[i])];
                this->bestiConformer.push_back(this->allGeneratediConformer[std::get<0>(scoreAndIndices[i])]);
            }

        }

        std::shared_ptr<DockingResult> ConformerRigidDockingEngine::getDockingResult() {
            auto ret = std::make_shared<DockingResult>();
            for (const auto &confr : this->bestiConformer) {
                Molecule finalLigand = this->orig_ligand->deepcopy();
                finalLigand.updateAtomPositionsFromiConformer(confr);
                ret->ligandPoses.emplace_back(finalLigand);
            }
            return ret;
        }


        bool ConformerRigidDockingEngine::setDockingBox(AbstractDockingEngine::DockingBoxSetting setting) {
            this->dockBoxSettings = setting;

            if (!(setting.type == DockingBoxSetting::Type::everything ||
                  setting.type == DockingBoxSetting::Type::centeredAround)) {
                BOOST_LOG_TRIVIAL(error) << "The passed DockingBoxSetting is not yet implemented.";
                BOOST_LOG_TRIVIAL(error) << "Running as if DockingBoxSetting::everything was passed";
                this->dockBoxSettings.type = DockingBoxSetting::Type::everything;
                return false;
            }
            return true;
        }


    }

}