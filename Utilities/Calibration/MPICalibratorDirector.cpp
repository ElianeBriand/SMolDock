//
// Created by briand on 3/21/19.
//

#include "MPICalibratorDirector.h"
#include "MPICalibratorCommon.h"


#include <cmath>
#include <thread>
#include <sstream>

#include <boost/log/trivial.hpp>

#define USE_BOOST_KARMA

#include <bprinter/table_printer.h>

#include <armadillo>
#include <ensmallen.hpp>

#include <Engines/Internals/InternalsUtilityFunctions.h>

#include <boost/serialization/string.hpp>

namespace SmolDock::Calibration {

    MPICalibratorDirector::MPICalibratorDirector(mpi::environment& env_,
                                                 mpi::communicator& world_,
                                                 Score::ScoringFunctionType scFuncType,
                                                 Heuristics::GlobalHeuristicType heurType,
                                                 Optimizer::LocalOptimizerType localOptimizerType_,
                                                 unsigned int maxLearningSteps,
                                                 double initialLearningRate_,
                                                 unsigned int rngSeed,
                                                 unsigned int conformerNumber,
                                                 unsigned int retryNumber,
                                                 unsigned int batchSize_,
                                                 Heuristics::HeuristicParameters hParams)
            :
            Calibrator(scFuncType,
                       heurType,
                       localOptimizerType_,
                       maxLearningSteps,
                       initialLearningRate_,
                       rngSeed, conformerNumber,
                       retryNumber,
                       batchSize_,
                       hParams),
            env(env_),
            world(world_) {
        this->numProcess = this->world.size();

    }

    bool MPICalibratorDirector::setupCalibration() {

        for (unsigned int j = 1; j < this->numProcess; ++j) {
            RegisterMessage rmsg;
            this->world.recv(j, MTags::RegisterNode, rmsg);
            this->nodeInfo.push_back(rmsg);
            this->numTotalCPU += rmsg.numThread;

            BOOST_LOG_TRIVIAL(info) << "Registed node " << j << " with " << rmsg.numThread << " threads.";
        }

        std::uniform_int_distribution<int> dis_int(0, std::numeric_limits<int>::max());


        this->workStructure.numReceptor = this->pdbBlockStrings.size();
        this->workStructure.conformerPerLigand = this->conformerNumber;
        this->workStructure.retryPerConformer = this->retryNumber;
        this->workStructure.numRegisteredVariant = this->variantsForAllLigands.size();
        this->workStructure.totalNumLigand = this->numLigand;
        this->workStructure.seed = dis_int(this->rndGenerator);
        this->workStructure.numCoeffToCalibrate = this->idxOfCoeffsToCalibrate.size();
        this->workStructure.idxOfCoeffsToCalibrate = this->idxOfCoeffsToCalibrate;
        this->workStructure.scoringFunctionType =this->scoringFunctionType;
        this->workStructure.heuristicType =this->heuristicType;
        this->workStructure.localOptimizerType =this->localOptimizerType;
        this->workStructure.differentialEpsilon = 1e-3;

        mpi::broadcast(world, this->workStructure, 0);


        for (unsigned int k = 0; k < this->workStructure.numReceptor; ++k) {
            ReceptorRecord rr;
            rr.PDBBlock = this->pdbBlockStrings[k];
            rr.dbsetting = this->dbSettings[k];
            mpi::broadcast(world, rr, 0);
        }

        for (unsigned int k = 0; k < this->workStructure.numRegisteredVariant; ++k) {
            mpi::broadcast(world, this->variantsForAllLigands[k], 0);
        }

        for (unsigned int k = 0; k < this->current_max_ReceptorID; ++k) {
            for (unsigned int j = 0; j < this->ligandSmilesRefScore[k].size(); ++j) {
                LigandRecord lr;
                lr.receptorID = k;
                lr.smiles = std::get<std::string>(this->ligandSmilesRefScore[k][j]);
                lr.deltaG = std::get<double>(this->ligandSmilesRefScore[k][j]);
                mpi::broadcast(world, lr, 0);
            }
            this->numLigandPerReceptor.push_back(this->ligandSmilesRefScore[k].size());
        }


        for (unsigned int l = 0; l < this->numLigand; ++l) {
            this->indexShuffler.push_back(l);
        }



        for (unsigned int j = 1; j < this->numProcess; ++j) {
            int status;
            this->world.recv(j, MTags::ReadyForWork, status);
            BOOST_LOG_TRIVIAL(info) << "Node " << j << " ready. (status = " << status <<")";
        }




        return true;
    }

    bool MPICalibratorDirector::runCalibration() {



        arma::mat coeffs_internalRepr = arma::mat(this->idxOfCoeffsToCalibrate.size(), 1);
        for (unsigned int k = 0; k < coeffs_internalRepr.n_rows; ++k) {
            coeffs_internalRepr[k] = this->currentCoeffs[this->idxOfCoeffsToCalibrate[k]];
        }

        ens::Adam optimizer(this->initialLearningRate, this->batchSize, 0.9, 0.999, 1e-8, this->maxLearningSteps, 1e-4, true);
        optimizer.Optimize(*this, coeffs_internalRepr);

        this->optResultMat = coeffs_internalRepr;

        std::vector<double> updateCoeffsVector = this->currentCoeffs;

        for (unsigned int j = 0; j < this->idxOfCoeffsToCalibrate.size(); ++j) {
            unsigned int idxCoeff = this->idxOfCoeffsToCalibrate[j];
            std::string coeffName = this->nameOfAllCoeffs[idxCoeff];
            double nonUpdatedCoeff = this->currentCoeffs[idxCoeff];

            updateCoeffsVector[idxCoeff] = coeffs_internalRepr[j];

            BOOST_LOG_TRIVIAL(info) << "COEFFICIENT " << coeffName << " ---- ";
            BOOST_LOG_TRIVIAL(info) << "   Old coeff        : " << nonUpdatedCoeff;
            BOOST_LOG_TRIVIAL(info) << "   New coeff        : " << coeffs_internalRepr[j];
            BOOST_LOG_TRIVIAL(info) << " ----------------------------- \n";
        }

        for (unsigned int j = 1; j < this->numProcess; ++j) {
            Task t;
            t.endOfTasks = true;
            this->world.send(j, MTags::TaskRequest, t);

            int status;
            this->world.recv(j, MTags::NodeTerminating, status);
            BOOST_LOG_TRIVIAL(info) << "Node " << j << " terminated. (status = " << status <<")";
        }


        return true;
    }


    bool MPICalibratorDirector::addReferenceLigand_SMILES_Ki(Calibrator::ReceptorID recID, const std::string& smiles,
                                                             double Ki, int seed) {

        const double R = 8.3144598;
        const double T = 310.15; // 37 degree celsius
        double deltaG = R * T * std::log(Ki) / 1000; // kcal/mol aka same as other docking software

        this->ligandSmilesRefScore[recID].emplace_back(std::make_tuple(smiles,deltaG));
        this->numLigand++;
        return true;
    }

    bool MPICalibratorDirector::addReferenceLigand_Mol_Ki(Calibrator::ReceptorID recID, const Molecule& mol, double Ki,
                                                          int seed) {

        BOOST_LOG_TRIVIAL(error) << "Cannot add ligand from Molecule class in MPICalibrator.";
        std::terminate();
        return Calibrator::addReferenceLigand_Mol_Ki(recID, mol, Ki, seed);
    }

    Calibrator::ReceptorID MPICalibratorDirector::addReceptorFromFile(const std::string& filename,
                                                                      Engine::AbstractDockingEngine::DockingBoxSetting dbsettings) {
        std::ifstream input(filename);
        std::stringstream sstr;
        sstr << input.rdbuf();

        this->pdbBlockStrings.push_back(sstr.str());
        this->dbSettings.push_back(dbsettings);

        this->current_max_ReceptorID++;
        return this->current_max_ReceptorID - 1;
    }

    Calibrator::ReceptorID MPICalibratorDirector::addReceptor(const Protein& prot,
                                                              Engine::AbstractDockingEngine::DockingBoxSetting dbsettings) {

        BOOST_LOG_TRIVIAL(error) << "Cannot add receptor from Protein class in MPICalibrator.";
        std::terminate();

        return Calibrator::addReceptor(prot, dbsettings);
    }

    bool MPICalibratorDirector::applyVariantToAllLigands(const std::string& SMARTSPattern, Atom::AtomVariant variant) {
        RegisteredVariant r;
        r.smarts = SMARTSPattern;
        r.atomVariant = static_cast<unsigned int>(variant);
        this->variantsForAllLigands.push_back(r);
        return true;
    }

    double MPICalibratorDirector::EvaluateWithGradient(const arma::mat& x,
                                                       const size_t i,
                                                       arma::mat& g,
                                                       const size_t batchSize) {
        static unsigned int batchCount = 0;
        batchCount++;

        unsigned int currentSendRank = 1;
        std::vector<mpi::request> taskRequests;

        std::vector<Task> tasks;
        for (unsigned int k = i; k < i + batchSize; ++k) {
            unsigned int realIdx = this->indexShuffler[k];
            std::tuple<unsigned int, unsigned int> rdlidx = RecLigIdxFromGlobalIdx(realIdx);

            auto& t = tasks.emplace_back((Task()));
            t.receptorID = std::get<0>(rdlidx);
            t.ligandIdx = std::get<1>(rdlidx);

            taskRequests.push_back(this->world.isend(currentSendRank, MTags::TaskRequest, t));

            currentSendRank++;
            if(currentSendRank == this->numProcess)
                currentSendRank = 1;
        }

        BOOST_LOG_TRIVIAL(debug) << "[Batch " << batchCount << "] Task request queued.";
        currentSendRank = 1;
        std::vector<Result> resultList;
        for (unsigned int k = i; k < i + batchSize; ++k) {

            auto& r = resultList.emplace_back((Result()));
            taskRequests.push_back(this->world.irecv(currentSendRank, MTags::ResultMessage, r));

            currentSendRank++;
            if(currentSendRank == this->numProcess)
                currentSendRank = 1;
        }

        BOOST_LOG_TRIVIAL(debug) << "[Batch " << batchCount << "] Result recv queued.";
        mpi::wait_all(taskRequests.begin(), taskRequests.end());
        BOOST_LOG_TRIVIAL(debug) << "[Batch " << batchCount << "] Results received.";

        const unsigned int numItems = resultList.size();
        const unsigned int numCoefficients = x.n_rows;
        double totalLoss = 0.0;
        std::vector<double> totalGradientLoss(static_cast<size_t>(numCoefficients)); // this constructor will initialize to 0.0
        for(const Result& result : resultList) {
            totalLoss += result.loss;
            for (unsigned int j = 0; j < numCoefficients; ++j) {
                totalGradientLoss.at(j) += result.lossGradient.at(j);
            }
        }

        const double meanLoss = totalLoss / numItems;
        for (unsigned int j = 0; j < numCoefficients; ++j) {
            g(j) = totalGradientLoss.at(j) / numItems;
        }

        BOOST_LOG_TRIVIAL(debug) << "[Batch " << batchCount << "] Mean loss : " << meanLoss;

        return meanLoss;
    }

    void MPICalibratorDirector::Shuffle() {
        std::shuffle(this->indexShuffler.begin(), this->indexShuffler.end(), this->rndGenerator);
    }

    size_t MPICalibratorDirector::NumFunctions() {
        return this->numLigand;
    }

    std::tuple<unsigned int, unsigned int> MPICalibratorDirector::RecLigIdxFromGlobalIdx(unsigned int idx) {
        unsigned int currRecId = 0;
        for(unsigned int& numLigForRec : this->numLigandPerReceptor)
        {
            if(idx - numLigForRec < 0)
            {
                return std::make_tuple(currRecId, idx);
            }
            idx = idx - numLigForRec;
            currRecId++;
        }
        BOOST_LOG_TRIVIAL(error) << "Indexing error : unable to translate global index " << idx;
        BOOST_LOG_TRIVIAL(error) << "   Total num ligand : " << this->numLigand;
        std::terminate();
        return std::make_tuple(currRecId, idx);
    }


}