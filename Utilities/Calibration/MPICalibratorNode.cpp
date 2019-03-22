//
// Created by briand on 3/21/19.
//

#include "MPICalibratorNode.h"


#include <thread>

#include "MPICalibratorCommon.h"
#include <boost/log/trivial.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>

namespace SmolDock::Calibration {


    MPICalibratorNode::MPICalibratorNode(mpi::environment& env_, mpi::communicator& world_) :
            env(env_),
            world(world_) {
    }


    void MPICalibratorNode::runNode() {
        RegisterMessage rmsg;

        rmsg.numThread = std::thread::hardware_concurrency();
        this->rank = world.rank();

        this->world.send(0, MTags::RegisterNode, rmsg);

        mpi::broadcast(world, this->workStructure, 0);

        this->rndGenerator.seed(this->workStructure.seed);
        std::uniform_int_distribution<> dis_int(0, std::numeric_limits<int>::max());
        arma::arma_rng::set_seed(this->workStructure.seed);

        for (unsigned int k = 0; k < this->workStructure.numReceptor; ++k) {
            ReceptorRecord rr;

            mpi::broadcast(world, rr, 0);

            this->pdbBlockStrings.push_back(rr.PDBBlock);
            this->dbSettings.push_back(rr.dbsetting);
        }

        for (unsigned int k = 0; k < this->workStructure.numRegisteredVariant; ++k) {
            RegisteredVariant rv;
            mpi::broadcast(world, rv, 0);
            this->variantsForAllLigands.push_back(rv);
        }


        for (unsigned int j = 0; j < this->workStructure.totalNumLigand; ++j) {
            LigandRecord lr;
            mpi::broadcast(world, lr, 0);

            this->ligandSmilesRefScore[lr.receptorID].push_back(
                    std::make_tuple(lr.smiles, lr.deltaG, nullptr, std::vector<iConformer>()));
        }


        for (auto it = this->ligandSmilesRefScore.begin(); it != this->ligandSmilesRefScore.end(); ++it) {
            for (unsigned int j = 0; j < it->second.size(); ++j) {
                std::shared_ptr<Molecule> mol_sptr = std::get<std::shared_ptr<Molecule>>(it->second.at(j));
                mol_sptr = std::make_shared<Molecule>(true);
                auto smiles = std::get<std::string>(it->second.at(j));
                mol_sptr->populateFromSMILES(smiles);

                auto& conformer_vector = std::get<std::vector<iConformer>>(it->second.at(j));
                mol_sptr->generateConformers(conformer_vector, this->workStructure.conformerPerLigand, true,
                                             dis_int(this->rndGenerator));

            }
        }

        for (unsigned int l = 0; l < this->pdbBlockStrings.size(); ++l) {

            std::shared_ptr<Protein> psptr = std::make_shared<Protein>();
            psptr->populateFromPDBString(this->pdbBlockStrings[l]);
            Engine::AbstractDockingEngine::DockingBoxSetting settings = this->dbSettings[l];
            iProtein iProt;
            iProtein fulliProt;

            if (settings.type == Engine::AbstractDockingEngine::DockingBoxSetting::Type::centeredAround) {
                iProt = psptr->getPartialiProtein_sphere(settings.center, settings.radius, 2.0);
            } else {
                iProt = psptr->getiProtein();
            }

            fulliProt = psptr->getiProtein();

            referenceReceptor.push_back(std::make_tuple(
                    psptr, settings, iProt, fulliProt
            ));
        }


        int status = 0;
        this->world.send(0, MTags::ReadyForWork, status);


        while (1) {

            Task t;
            this->world.recv(0, MTags::TaskRequest, t);

            if (t.endOfTasks) {
                this->world.send(0, MTags::NodeTerminating, 0);
                break;
            }


            const unsigned int& recId = t.receptorID;
            const unsigned int& ligandIdx = t.ligandIdx;
            std::vector<double>& coeffs = t.coefficients;

            const unsigned int& conformerPerLigand = this->workStructure.conformerPerLigand;
            const unsigned int& retryPerConformer = this->workStructure.retryPerConformer;

            const auto& confvector = std::get<std::vector<iConformer>>(this->ligandSmilesRefScore.at(recId).at(ligandIdx));
            const auto& referenceScore = std::get<double>(this->ligandSmilesRefScore.at(recId).at(ligandIdx));
            const iProtein& prot = std::get<2>(this->referenceReceptor[recId]);
            const iProtein& fullprot = std::get<3>(this->referenceReceptor[recId]);

            auto resultMutex = std::make_shared<std::mutex>();
            auto local_scores = std::make_shared<std::vector<double>>();

            tbb::parallel_for(tbb::blocked_range2d<size_t>(0, conformerPerLigand, 1, 0, retryPerConformer, 1),
                              MPICalibrator2DLoopRunner(
                                      this->workStructure,
                                      confvector,
                                      prot,
                                      fullprot,
                                      resultMutex,
                                      local_scores,
                                      coeffs)
                              );

            double score = *(std::min_element(std::begin(*local_scores), std::end(*local_scores)));

            std::vector<double> gradientVector;
            for(unsigned int& idx : this->workStructure.idxOfCoeffsToCalibrate)
            {
                std::vector<double> epsilonCoeffs = coeffs;
                epsilonCoeffs[idx] += this->workStructure.differentialEpsilon;
                auto local_scores_gradient = std::make_shared<std::vector<double>>();

                tbb::parallel_for(tbb::blocked_range2d<size_t>(0, conformerPerLigand, 1, 0, retryPerConformer, 1),
                                  MPICalibrator2DLoopRunner(
                                          this->workStructure,
                                          confvector,
                                          prot,
                                          fullprot,
                                          resultMutex,
                                          local_scores_gradient,
                                          epsilonCoeffs)
                );

                double gradientScore = *(std::min_element(std::begin(*local_scores_gradient), std::end(*local_scores_gradient)));
                gradientVector.push_back(gradientScore);
            }



            Result r;
            r.loss = score - referenceScore;
            r.lossGradient = std::move(gradientVector);

            this->world.send(0, MTags::ResultMessage, r);
        }

    }


    MPICalibrator2DLoopRunner::MPICalibrator2DLoopRunner(const WorkStructure& ws,
            const std::vector<iConformer>& conformerVector_,
                                                              const iProtein& protein_,
                                                              const iProtein& fullProtein_,
                                                              std::shared_ptr<std::mutex> resultMutex_,
                                                              std::shared_ptr<std::vector<double>> local_scores_,
                                                              std::vector<double> currentCoeffs_):
            workStructure(ws),
            conformerVector(conformerVector_),
            protein(protein_),
            fullProtein(fullProtein_),
            resultMutex(resultMutex_),
            local_scores(local_scores_),
            currentCoeffs(currentCoeffs_) {
    }

    void MPICalibrator2DLoopRunner::operator()(const tbb::blocked_range2d<size_t>& r) const {
        for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
            const iConformer& conformer = conformerVector[i];

            iTransform tr;
            tr.transl.x() += conformer.centroidNormalizingTransform.x();
            tr.transl.y() += conformer.centroidNormalizingTransform.y();
            tr.transl.z() += conformer.centroidNormalizingTransform.z();
            tr.rota.normalize();

            std::shared_ptr<Score::ScoringFunction> scoringFunction = scoringFunctionFactory(
                    this->workStructure.scoringFunctionType,
                    conformer,
                    this->protein,
                    tr,
                    1e-3,
                    true);

            scoringFunction->setNonDefaultCoefficients(this->currentCoeffs);


            std::shared_ptr<Score::ScoringFunction> fullScoringFunction = scoringFunctionFactory(
                    this->workStructure.scoringFunctionType,
                    conformer,
                    this->fullProtein,
                    tr,
                    1e-3);

            fullScoringFunction->setNonDefaultCoefficients(this->currentCoeffs);

            std::shared_ptr<Optimizer::Optimizer> localOptimizer = optimizerFactory(
                    this->workStructure.localOptimizerType,
                    scoringFunction.get(),
                    1e-3);

            for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
                std::shared_ptr<Heuristics::GlobalHeuristic> globalHeuristic = globalHeuristicFactory(
                        this->workStructure.heuristicType,
                        scoringFunction.get(),
                        localOptimizer.get(),
                        this->workStructure.seed + j,
                        Heuristics::emptyParameters);

                globalHeuristic->search();

                auto rawResultMatrix = globalHeuristic->getResultMatrix();
                double score = scoringFunction->EvaluateOnlyIntermolecular(rawResultMatrix);
                double fullScore = fullScoringFunction->EvaluateOnlyIntermolecular(rawResultMatrix);
                double delta_full = fullScore / score;

                if (delta_full > 1.2 || delta_full < 0.80) {
                    continue;
                }

                {
                    std::lock_guard lock(*this->resultMutex);
                    this->local_scores->push_back(fullScore);
                }
            }
        }
    }


}