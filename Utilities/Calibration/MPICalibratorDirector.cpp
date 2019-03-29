//
// Created by briand on 3/21/19.
//

#include "MPICalibratorDirector.h"
#include "MPICalibratorCommon.h"


#include <cmath>
#include <thread>
#include <chrono>
#include <sstream>

#include <boost/log/trivial.hpp>

#define USE_BOOST_KARMA

#include <bprinter/table_printer.h>

#include <armadillo>
#include <ensmallen.hpp>

#include <Engines/Internals/InternalsUtilityFunctions.h>
#include <Utilities/LogUtils.h>

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

        // Determining optimal allocation of workItem depending on available CPU on each node

        int remainingWorkItemToAllocate = this->batchSize;
        for (unsigned int j = 0; j < this->nodeInfo.size(); ++j) {
            double fractionOfWorkload = static_cast<double>(this->nodeInfo.at(j).numThread) / static_cast<double>(numTotalCPU);
            unsigned int numWorkItem = static_cast<unsigned int>(this->batchSize * fractionOfWorkload);
            this->numWorkItemPerNode.push_back(numWorkItem);
            remainingWorkItemToAllocate -= numWorkItem;
        }

        if(remainingWorkItemToAllocate >= 0)
        {
            // Spread out the remaining work items
            unsigned int currRank = 1;
            while(remainingWorkItemToAllocate > 0) {
                this->numWorkItemPerNode.at(currRank) += 1;
                remainingWorkItemToAllocate--;
                if(currRank == this->nodeInfo.size() - 1)
                    currRank = 1;
                else
                    currRank++;
            }
        }else {
            BOOST_LOG_TRIVIAL(error) << "Unexpected oversubsciption of work item per batch.";
            BOOST_LOG_TRIVIAL(error) << "   Batch size : " << this->batchSize;
            for (unsigned int j = 0; j < this->numWorkItemPerNode.size(); ++j) {
                BOOST_LOG_TRIVIAL(error) << "   - Node " << j+1 << " : " << this->numWorkItemPerNode.at(j) << " work items/batch";
            }
            BOOST_LOG_TRIVIAL(error) << "   Remaining work item : " << remainingWorkItemToAllocate;
            std::terminate();
        }


        BOOST_LOG_TRIVIAL(info) << "Batch mapping information :";
        BOOST_LOG_TRIVIAL(info) << "   Batch size : " << this->batchSize;
        for (unsigned int j = 0; j < this->numWorkItemPerNode.size(); ++j) {
            BOOST_LOG_TRIVIAL(info) << "   - Node " << j+1 <<
            " : " << this->nodeInfo.at(j).numThread <<
            " CPU for " << this->numWorkItemPerNode.at(j) << " work items/batch";
        }




        //this->numWorkItemPerNode

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
            rr.specialResTypes = this->specialResidueTypings[k];
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

        calibrationStillRunning = true;
        std::thread statusPrinter(&MPICalibratorDirector::updateAndPrintStatus, this);



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

        calibrationStillRunning = false;
        statusPrinter.join();



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
        this->specialResidueTypings.push_back(std::vector<MPISpecialResidueTyping>());

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


        // Basically this is an allocation mecanism to have an amount of work item per node
        // roughly proportional to the number of cpu on each node
        // (We have collected the CPU #/node when initializing.)
        // (We use intel TBB to spread the work on all CPU available per node, so this allow
        // us to have around the same runtime for each node, so we dont wait for the weakest node)
        // TODO : we need to handle slurm-type job manager, where we wont necessarily be able to just
        // use all the CPUs of the node.
        unsigned int currentSendRank = 1; // Careful : this is (rightfully) off-by-one from our array access index
        std::vector<unsigned int> remainingWorkItemCapacity = this->numWorkItemPerNode;

        std::vector<mpi::request> taskRequests;
        std::vector<Task> tasks;

        std::vector<double> currCoeffsUpdated = this->currentCoeffs;
        for (unsigned int k = 0; k < this->idxOfCoeffsToCalibrate.size(); ++k) {
            currCoeffsUpdated[this->idxOfCoeffsToCalibrate[k]] = x(k);
        }
        coefficientHistory.push_back(currCoeffsUpdated);

        for (unsigned int k = i; k < i + batchSize; ++k) {

            if(remainingWorkItemCapacity[currentSendRank - 1] == 0) {
                // The current rank has no more cpu for us to run
                if(currentSendRank == remainingWorkItemCapacity.size()) {
                    // We have consumed every available cpu for this cycle, we start a new one
                    currentSendRank = 1;
                    remainingWorkItemCapacity = this->numWorkItemPerNode;
                }else{
                    // We use the next node
                    currentSendRank++;
                }
            }
            // In any cases, we will consume one cpu for the current item
            remainingWorkItemCapacity[currentSendRank - 1]--;


            unsigned int realIdx = this->indexShuffler[k];
            std::tuple<unsigned int, unsigned int> rdlidx = RecLigIdxFromGlobalIdx(realIdx);

            auto& t = tasks.emplace_back((Task()));
            t.receptorID = std::get<0>(rdlidx);
            t.ligandIdx = std::get<1>(rdlidx);

            t.coefficients = currCoeffsUpdated;

            taskRequests.push_back(this->world.isend(currentSendRank, MTags::TaskRequest, t));

        }

        BOOST_LOG_TRIVIAL(debug) << "[Epoch " << batchCount << "] Task request queued.";

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        currentSendRank = 1;
        remainingWorkItemCapacity = this->numWorkItemPerNode;

        std::vector<Result> resultList;
        for (unsigned int k = i; k < i + batchSize; ++k) {
            resultList.push_back((Result()));
            Result& r = resultList.back();
            this->world.recv(boost::mpi::any_source, MTags::ResultMessage, r);
            BOOST_LOG_TRIVIAL(debug) << "Received " << k - i  << " of " << (batchSize - i) << " results.";
        }

        end = std::chrono::system_clock::now();
        int elapsed_minutes = std::chrono::duration_cast<std::chrono::minutes>
                (end-start).count();
        BOOST_LOG_TRIVIAL(debug) << "[Epoch " << batchCount << "] All results received.";

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
        this->gradientHistory.push_back(std::vector<double>());
        for (unsigned int j = 0; j < numCoefficients; ++j) {
            const double meanGradComponentValue = totalGradientLoss.at(j) / numItems;
            g(j) = meanGradComponentValue;
            this->gradientHistory.back().push_back(meanGradComponentValue);
        }

        this->lossHistory.push_back(meanLoss);

        BOOST_LOG_TRIVIAL(debug) << "[Epoch " << batchCount << "] Mean loss : " << meanLoss;
        BOOST_LOG_TRIVIAL(debug) << "   Duration : " << elapsed_minutes << " minutes";
        BOOST_LOG_TRIVIAL(debug) << "   History : ";
        for (unsigned int l = 0; l < this->lossHistory.size(); ++l) {
            BOOST_LOG_TRIVIAL(debug) << "     Epoch " << l << "   ->  " << this->lossHistory.at(l);
            BOOST_LOG_TRIVIAL(debug) << "       --> Coeffs " << vectorToString(this->coefficientHistory.at(l));
            BOOST_LOG_TRIVIAL(debug) << "       --> Gradient " << vectorToString(this->gradientHistory.at(l));

        }

        batchCount++;


        BOOST_LOG_TRIVIAL(info) << "\n\n    ------\n\n";

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
            if(static_cast<long int>(idx) - static_cast<long int>(numLigForRec) < 0)
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

    bool MPICalibratorDirector::applySpecialResidueTypingFromRecID(Calibrator::ReceptorID recID,
                                                                   const AminoAcid::AAType resType,
                                                                   const unsigned int serialNumber,
                                                                   const SpecialResidueTyping specialType) {
        MPISpecialResidueTyping srt;
        srt.aaType = static_cast<unsigned int>(resType);
        srt.serialNumber = serialNumber;
        srt.specialTyping = static_cast<unsigned int>(specialType);
        this->specialResidueTypings[static_cast<unsigned int>(recID)].push_back(srt);
        return false;
    }

    void MPICalibratorDirector::updateAndPrintStatus() {
        while(calibrationStillRunning == true)
        {
            using namespace std::chrono_literals;
            std::this_thread::sleep_for(10s);
        }

    }


}