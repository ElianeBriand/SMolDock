//
// Created by briand on 3/21/19.
//

#include "MPICalibratorNode.h"


#include <thread>

#include "MPICalibratorCommon.h"
#include <boost/log/trivial.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>

#include <Engines/Internals/InternalsUtilityFunctions.h>
#include <Utilities/LogUtils.h>
#include <Structures/Atom.h>
#include <Structures/AminoAcid.h>


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
            this->specialResidueTypings.push_back(rr.specialResTypes);
        }

        for (unsigned int k = 0; k < this->workStructure.numRegisteredVariant; ++k) {
            RegisteredVariant rv;
            mpi::broadcast(world, rv, 0);
            this->variantsForAllLigands.push_back(rv);

        }


        for (unsigned int j = 0; j < this->workStructure.totalNumLigand; ++j) {
            LigandRecord lr;
            mpi::broadcast(world, lr, 0);

            if(lr.isAnchor == false) {
            this->ligandSmilesRefScore[lr.receptorID].push_back(
                    std::make_tuple(lr.smiles, lr.deltaG, nullptr, std::vector<iConformer>()));
            } else {
                this->anchorLigands[lr.receptorID].push_back(
                        std::make_tuple(lr.molblock, nullptr, nullptr, std::vector<iConformer>()));
            }

        }


        for (auto it = this->ligandSmilesRefScore.begin(); it != this->ligandSmilesRefScore.end(); ++it) {
            for (unsigned int j = 0; j < it->second.size(); ++j) {
                std::shared_ptr<Molecule> mol_sptr = std::get<std::shared_ptr<Molecule>>(it->second.at(j));
                mol_sptr = std::make_shared<Molecule>(true);
                auto smiles = std::get<std::string>(it->second.at(j));
                mol_sptr->populateFromSMILES(smiles);

                for(RegisteredVariant& variant : this->variantsForAllLigands)
                {
                    mol_sptr->applyAtomVariant(variant.smarts,static_cast<Atom::AtomVariant>(variant.atomVariant));
                }

                auto& conformer_vector = std::get<std::vector<iConformer>>(it->second.at(j));
                mol_sptr->generateConformers(conformer_vector, this->workStructure.conformerPerLigand, true,
                                             dis_int(this->rndGenerator));

            }
        }

        for (auto it = this->anchorLigands.begin(); it != this->anchorLigands.end(); ++it) {
            for (unsigned int j = 0; j < it->second.size(); ++j) {
                std::shared_ptr<Molecule> mol_ref_sptr = std::get<1>(it->second.at(j));
                std::shared_ptr<Molecule> mol_sptr = std::get<2>(it->second.at(j));
                mol_ref_sptr = std::make_shared<Molecule>(true);
                mol_sptr = std::make_shared<Molecule>(true);
                auto molblock = std::get<std::string>(it->second.at(j));
                bool res = mol_sptr->populateFromMolBlock(molblock);
                bool res2 = mol_ref_sptr->populateFromMolBlock(molblock);
                if(!res || !res2) {
                    BOOST_LOG_TRIVIAL(error) << "Cannot parse transmitted mol block to Molecule !";
                }

                for(RegisteredVariant& variant : this->variantsForAllLigands)
                {
                    mol_sptr->applyAtomVariant(variant.smarts,static_cast<Atom::AtomVariant>(variant.atomVariant));
                }

                auto& conformer_vector = std::get<std::vector<iConformer>>(it->second.at(j));
                mol_sptr->generateConformers(conformer_vector, this->workStructure.conformerPerLigand, true,
                                             dis_int(this->rndGenerator));
            }
        }


        for (unsigned int l = 0; l < this->pdbBlockStrings.size(); ++l) {

            std::shared_ptr<Protein> psptr = std::make_shared<Protein>();
            psptr->populateFromPDBString(this->pdbBlockStrings[l]);

            for(const MPISpecialResidueTyping& srt: this->specialResidueTypings[l]) {
                psptr->applySpecialResidueTyping(
                        static_cast<AminoAcid::AAType>(srt.aaType),
                        srt.serialNumber,
                        static_cast<SpecialResidueTyping>(srt.specialTyping)
                        );

            }


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


        BOOST_LOG_TRIVIAL(debug) << "Node ready for work.";

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
            bool isAnchorTask = t.anchorLigandTask;

            const unsigned int& conformerPerLigand = this->workStructure.conformerPerLigand;
            const unsigned int& retryPerConformer = this->workStructure.retryPerConformer;

            BOOST_LOG_TRIVIAL(debug) << "Received one task";
            BOOST_LOG_TRIVIAL(debug) << " | Conformer          : " << conformerPerLigand;
            BOOST_LOG_TRIVIAL(debug) << " | Retry   nb         : " << retryPerConformer;
            BOOST_LOG_TRIVIAL(debug) << " | Coefficients       : " << vectorToString(coeffs);
            BOOST_LOG_TRIVIAL(debug) << " | Coeffs to change   : " << vectorToString(this->workStructure.idxOfCoeffsToCalibrate);
            BOOST_LOG_TRIVIAL(debug) << " | Receptor idx       : " << recId;
            BOOST_LOG_TRIVIAL(debug) << " | Ligand idx         : " << ligandIdx;
            BOOST_LOG_TRIVIAL(debug) << " | Is anchor ?        : " << isAnchorTask;

            const auto& confvector =  isAnchorTask ? std::get<std::vector<iConformer>>(this->anchorLigands.at(recId).at(ligandIdx))
                    : std::get<std::vector<iConformer>>(this->ligandSmilesRefScore.at(recId).at(ligandIdx));
            const iProtein& prot = std::get<2>(this->referenceReceptor[recId]);
            const iProtein& fullprot = std::get<3>(this->referenceReceptor[recId]);

            if(isAnchorTask) {
                // Anchor task : compute RMSD and loss
                // TODO ANCHOR


                std::vector<double> gradientVector;
                for(unsigned int& idx : this->workStructure.idxOfCoeffsToCalibrate)
                {
                     gradientVector.push_back(0.0);
                }

                Result r;
                r.loss =  0.0;
                r.lossGradient = gradientVector;

                this->world.send(0, MTags::ResultMessage, r);
            } else {
                //Non Anchor task : compute score and loss



            auto resultMutex = std::make_shared<std::mutex>();
            auto local_scores = std::make_shared<std::vector<double>>();


            tbb::parallel_for(tbb::blocked_range2d<size_t>(0, conformerPerLigand,0, retryPerConformer),
                              MPICalibrator2DLoopRunner(
                                      this->workStructure,
                                      confvector,
                                      prot,
                                      fullprot,
                                      resultMutex,
                                      local_scores,
                                      coeffs)
                              );

            // We select the top nth (n = conformerPerLigand) for averaging
            std::sort(local_scores->begin(), local_scores->end(), std::less<double>());
            double score = std::accumulate(local_scores->begin(), local_scores->begin() + conformerPerLigand, 0.0)/conformerPerLigand;

            BOOST_LOG_TRIVIAL(debug) << "Evaluation done : "
            << vectorToString(std::vector<double>(local_scores->begin(), local_scores->begin() + conformerPerLigand)) << " -> " << score;

            std::vector<double> gradientVector;
            for(unsigned int& idx : this->workStructure.idxOfCoeffsToCalibrate)
            {
                std::vector<double> epsilonCoeffs = coeffs;
                epsilonCoeffs[idx] += this->workStructure.differentialEpsilon;
                auto local_scores_gradient = std::make_shared<std::vector<double>>();

                tbb::parallel_for(tbb::blocked_range2d<size_t>(0, conformerPerLigand, 0, retryPerConformer),
                                  MPICalibrator2DLoopRunner(
                                          this->workStructure,
                                          confvector,
                                          prot,
                                          fullprot,
                                          resultMutex,
                                          local_scores_gradient,
                                          epsilonCoeffs)
                );

                // We select the top nth (n = conformerPerLigand) for averaging
                std::sort(local_scores_gradient->begin(), local_scores_gradient->end(), std::less<double>());
                double gradientScore = std::accumulate(local_scores_gradient->begin(), local_scores_gradient->begin() + conformerPerLigand, 0.0)/conformerPerLigand;
                gradientVector.push_back(gradientScore - score);
            }


            const auto& referenceScore = std::get<double>(this->ligandSmilesRefScore.at(recId).at(ligandIdx));


            BOOST_LOG_TRIVIAL(debug) << "Completed one task";
            BOOST_LOG_TRIVIAL(debug) << " | Conformer          : " << conformerPerLigand;
            BOOST_LOG_TRIVIAL(debug) << " | Retry              : " << retryPerConformer;
            BOOST_LOG_TRIVIAL(debug) << " | Reference score    : " << referenceScore;
            BOOST_LOG_TRIVIAL(debug) << " | Coefficients       : " << vectorToString(coeffs);
            BOOST_LOG_TRIVIAL(debug) << " | Score              : " << score;
            BOOST_LOG_TRIVIAL(debug) << " | Gradient           : " << vectorToString(gradientVector);

            Result r;
            r.loss = score - referenceScore;
            r.lossGradient = gradientVector;

            this->world.send(0, MTags::ResultMessage, r);

            }
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
            for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
                // BOOST_LOG_TRIVIAL(debug) << "Docking conformer " << i << " retry " << j;

                const iConformer& conformer = conformerVector[i];

                iTransform tr = iTransformIdentityInit(conformer.num_rotatable_bond);
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
                        1e-3,
                        true);

                fullScoringFunction->setNonDefaultCoefficients(this->currentCoeffs);

                std::shared_ptr<Optimizer::Optimizer> localOptimizer = optimizerFactory(
                        this->workStructure.localOptimizerType,
                        scoringFunction.get(),
                        1e-3);




                std::shared_ptr<Heuristics::GlobalHeuristic> globalHeuristic = globalHeuristicFactory(
                        this->workStructure.heuristicType,
                        scoringFunction.get(),
                        localOptimizer.get(),
                        this->workStructure.seed + j,
                        heuristicParametersFactory(this->workStructure.heuristicType));

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