//
// Created by briand on 3/21/19.
//

#include "MPICalibratorDirector.h"
#include "MPICalibratorCommon.h"


#include <cmath>
#include <thread>
#include <chrono>
#include <sstream>
#include <set>

#include <boost/log/trivial.hpp>

#define USE_BOOST_KARMA

#include <bprinter/table_printer.h>

#include <armadillo>
#include <ensmallen.hpp>

#include <Engines/Internals/InternalsUtilityFunctions.h>
#include <Utilities/LogUtils.h>


#include <bprinter/table_printer.h>

#include <boost/serialization/string.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

using boost::lexical_cast;

namespace SmolDock::Calibration {

    MPICalibratorDirector::MPICalibratorDirector(mpi::environment& env_,
                                                 mpi::communicator& world_,
                                                 Score::ScoringFunctionType scFuncType,
                                                 Heuristics::GlobalHeuristicType heurType,
                                                 Optimizer::LocalOptimizerType localOptimizerType_,
                                                 unsigned int maxLearningSteps,
                                                 double stepSize_,
                                                 unsigned int rngSeed,
                                                 unsigned int conformerNumber,
                                                 unsigned int retryNumber,
                                                 unsigned int batchSize_,
                                                 Heuristics::HeuristicParameters hParams,
                                                 const std::string& restoreArchivePrefix_)
            :
            Calibrator(scFuncType,
                       heurType,
                       localOptimizerType_,
                       maxLearningSteps,
                       stepSize_,
                       rngSeed, conformerNumber,
                       retryNumber,
                       batchSize_,
                       hParams),
            env(env_),
            world(world_),
            restoreArchivePrefix(restoreArchivePrefix_) {
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

        // Rebalancing slightly to avoid node with 0 work items
        // NB : This may not be optimal if the nodes are very imbalanced
        //      but it is expected that one would prefer to have all nodes
        //      do at least some work ? To be discussed if problem arises.

        // Find nodes with max workload
        unsigned int maxWorkItem = 0;
        std::vector<unsigned int> idxMaxWorkItem;
        for (unsigned int j = 0; j < this->nodeInfo.size(); ++j) {
            unsigned int numWorkItemForNode = this->numWorkItemPerNode.at(j);
            if(numWorkItemForNode > maxWorkItem) {
                maxWorkItem = numWorkItemForNode;
                idxMaxWorkItem.push_back(j);
            }
        }

        // Find nodes with 0 workload
        std::vector<unsigned int> idxZeroItem;
        for (unsigned int j = 0; j < this->nodeInfo.size(); ++j) {
            unsigned int numWorkItemForNode = this->numWorkItemPerNode.at(j);
            if(numWorkItemForNode == 0) {
                idxZeroItem.push_back(j);
            }
        }

        auto it_maxLoadedNodeIdx = idxMaxWorkItem.begin();
        for(auto zeroWorkNodeIdx : idxZeroItem) {
            if(this->numWorkItemPerNode.at(std::distance(idxMaxWorkItem.begin(), it_maxLoadedNodeIdx)) > 2) {
                this->numWorkItemPerNode.at(zeroWorkNodeIdx) += 1;
                this->numWorkItemPerNode.at(std::distance(idxMaxWorkItem.begin(), it_maxLoadedNodeIdx)) -= 1;
            }
            if((it_maxLoadedNodeIdx + 1) == idxMaxWorkItem.end()) {
                    // No more max loaded node to rebalance
                    // We continue to unload this one
                    continue;
            } else {
                it_maxLoadedNodeIdx++;
            }
        }

        // Final consistency check
        unsigned int sumWorkItem = 0;
        for (unsigned int j = 0; j < this->nodeInfo.size(); ++j) {
            unsigned int numWorkItemForNode = this->numWorkItemPerNode.at(j);
            sumWorkItem += numWorkItemForNode;
        }

        if(sumWorkItem != this->batchSize) {
            BOOST_LOG_TRIVIAL(error) << "Work load balancing unexpectedly failed ??";
            BOOST_LOG_TRIVIAL(error) << "   Batch size : " << this->batchSize;
            for (unsigned int j = 0; j < this->numWorkItemPerNode.size(); ++j) {
                BOOST_LOG_TRIVIAL(error) << "   - Node " << j+1 << " : " << this->numWorkItemPerNode.at(j) << " work items/batch";
            }
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
        this->workStructure.totalNumLigand = this->totalNumLigand;
        this->workStructure.numAnchorLigand = this->numAnchorLigand;
        this->workStructure.seed = dis_int(this->rndGenerator);
        this->workStructure.numCoeffToCalibrate = this->idxOfCoeffsToCalibrate.size();
        this->workStructure.idxOfCoeffsToCalibrate = this->idxOfCoeffsToCalibrate;
        this->workStructure.scoringFunctionType =this->scoringFunctionType;
        this->workStructure.heuristicType =this->heuristicType;
        this->workStructure.localOptimizerType =this->localOptimizerType;
        this->workStructure.differentialEpsilon = 0.5e-7;

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
                lr.isAnchor = false;
                lr.molblock = "";
                mpi::broadcast(world, lr, 0);
            }
            for (unsigned int j = 0; j < this->anchorLigands[k].size(); ++j) {
                LigandRecord lr;
                lr.receptorID = k;
                lr.smiles = this->anchorLigands[k][j]->writeToSMILES();
                lr.deltaG = 0.0;
                lr.isAnchor = true;
                lr.molblock = this->anchorLigands[k][j]->writeToMolBlock();
                mpi::broadcast(world, lr, 0);
            }
            this->numLigandPerReceptor.push_back(this->ligandSmilesRefScore[k].size());
            this->numAnchorLigandPerReceptor.push_back(this->anchorLigands[k].size());
        }


        for (unsigned int l = 0; l < this->totalNumLigand; ++l) {
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


        ens::Adam optimizer(this->stepSize, this->batchSize, 0.9, 0.999, 1e-8, this->maxLearningSteps, 1e-6, true);
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
                                                             double Ki) {

        const double R = 8.3144598;
        const double T = 310.15; // 37 degree celsius
        double deltaG = R * T * std::log(Ki) / 1000; // kcal/mol aka same as other docking software

        this->addReferenceLigand_SMILES_deltaG(recID,smiles,deltaG);
        return true;
    }

    bool MPICalibratorDirector::addReferenceLigand_SMILES_deltaG(Calibrator::ReceptorID recID,
                                                             const std::string& smiles,
                                                             double deltaG) {
        this->ligandSmilesRefScore[recID].emplace_back(std::make_tuple(smiles,deltaG));
        this->totalNumLigand++;
        return true;
    }

    bool MPICalibratorDirector::addReferenceLigand_Mol_Ki(Calibrator::ReceptorID recID, const Molecule& mol, double Ki) {

        BOOST_LOG_TRIVIAL(error) << "Cannot add ligand from Molecule class in MPICalibrator.";
        std::terminate();
        return Calibrator::addReferenceLigand_Mol_Ki(recID, mol, Ki);
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

        for (unsigned int k = i; k < (i + batchSize); ++k) {

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

            if(realIdx < (totalNumLigand - numAnchorLigand)) {
                // it a normal (non-anchor) ligand
                t.anchorLigandTask = false;
            } else {
                // its an anchor ligand
                t.anchorLigandTask = true;
            }

            taskRequests.push_back(this->world.isend(currentSendRank, MTags::TaskRequest, t));

        }

        BOOST_LOG_TRIVIAL(info) << "[Epoch " << this->batchCount << "] Task request queued.";

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        currentSendRank = 1;
        remainingWorkItemCapacity = this->numWorkItemPerNode;

        std::vector<Result> resultList;
        for (unsigned int k = i; k < (i + batchSize); ++k) {
            resultList.push_back((Result()));
            Result& r = resultList.back();
            this->world.recv(boost::mpi::any_source, MTags::ResultMessage, r);
            BOOST_LOG_TRIVIAL(debug) << "Received " << (k - i + 1)  << " of " << batchSize << " results.";
            BOOST_LOG_TRIVIAL(debug) << "  Result  ";
            BOOST_LOG_TRIVIAL(debug) << "     |   loss          :  " << r.loss;
            BOOST_LOG_TRIVIAL(debug) << "     |   loss gradient :  " << vectorToString(r.lossGradient);
        }

        end = std::chrono::system_clock::now();
        int elapsed_minutes = std::chrono::duration_cast<std::chrono::minutes>
                (end-start).count();
        BOOST_LOG_TRIVIAL(info) << "[Epoch " << this->batchCount << "] All results received.";
        this->durationHistory.push_back(elapsed_minutes);

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

        BOOST_LOG_TRIVIAL(info) << "[Epoch " << this->batchCount << "] Mean loss : " << meanLoss;
        BOOST_LOG_TRIVIAL(info) << "   Duration : " << elapsed_minutes << " minutes";
        BOOST_LOG_TRIVIAL(info) << "   History : ";

        bprinter::TablePrinter tp_lossDuration(&std::cout);
        tp_lossDuration.AddColumn("Epoch", 8);
        tp_lossDuration.AddColumn("Mean loss", 20);
        tp_lossDuration.AddColumn("Duration (min)", 20);

        tp_lossDuration.PrintHeader();
        for (unsigned int l = 0; l < this->lossHistory.size(); ++l) {
            tp_lossDuration << l << this->lossHistory.at(l) << this->durationHistory.at(l);

        }
        tp_lossDuration.PrintFooter();

        bprinter::TablePrinter tp_coeffsGradients(&std::cout);
        tp_coeffsGradients.AddColumn("Epoch", 8);
        tp_coeffsGradients.AddColumn("Metric", 12);
        tp_coeffsGradients.AddColumn("Value", 60);

        tp_coeffsGradients.PrintHeader();
        for (unsigned int l = 0; l < this->lossHistory.size(); ++l) {
            tp_lossDuration << l << "Loss" <<  this->lossHistory.at(l);
            tp_lossDuration << l << "Coeffs" <<  vectorToString(this->coefficientHistory.at(l));
            tp_lossDuration << l << "Gradient" <<  vectorToString(this->gradientHistory.at(l));
        }
        tp_coeffsGradients.PrintFooter();


        this->batchCount++;


        const boost::regex my_filter(this->restoreArchivePrefix + std::string("([0-9]+)\\.restore-state") );

        int max_num = -1;
        boost::filesystem::directory_iterator end_itr; // Default ctor yields past-the-end
        for( boost::filesystem::directory_iterator i("./"); i != end_itr; ++i )
        {
            // Skip if not a file
            if( !boost::filesystem::is_regular_file( i->status() ) ) continue;

            boost::smatch what;

            // Skip if no match for V2:
            //if( !boost::regex_match( i->leaf(), what, my_filter ) ) continue;
            // For V3:
            std::string path_as_string = i->path().filename().string();
            bool res = boost::regex_match( path_as_string, what, my_filter, boost::match_extra) ;
            if(!res) continue;
            int numFile = lexical_cast<int>(what[1]);
            if(numFile > max_num) {
                max_num = numFile;
            }
        }


        std::string filename;
        if(max_num != -1) {
            filename = this->restoreArchivePrefix + std::to_string(max_num + 1) + std::string(".restore-state");
        } else {
            filename = this->restoreArchivePrefix + std::to_string(1) + std::string(".restore-state");
        }
        BOOST_LOG_TRIVIAL(info) << "Saving to " << filename;

        std::ofstream ofs(filename);
        {
            boost::archive::text_oarchive oa(ofs);
            MPICalibratorDirector_ResumeObject resumeObject =  this->createResumeState();
            oa << resumeObject;
        }
        BOOST_LOG_TRIVIAL(info) << "Saved.";


        BOOST_LOG_TRIVIAL(info) << "\n\n    ------\n\n";

        return meanLoss;
    }

    void MPICalibratorDirector::Shuffle() {
        std::shuffle(this->indexShuffler.begin(), this->indexShuffler.end(), this->rndGenerator);
    }

    size_t MPICalibratorDirector::NumFunctions() {
        return this->totalNumLigand;
    }

    std::tuple<unsigned int, unsigned int> MPICalibratorDirector::RecLigIdxFromGlobalIdx(unsigned int idx) {
        /**
         * We have N receptor with M_n ligands each
         * The global index goes ligand 0 of receptor 0, ligand 1 of receptor 0 ... ligand M_1 of receptor 0
         * then ligand 0 of receptor 1, etc etc
         * So we use the numLigandPerReceptor to get back to the local index of the ligand for a given receptor
         *
         * Additionally, after the end of the normal (non-anchor) ligand, we have the same system for the anchor ligands
         *
         * NB: Code is duplicated for ease of comprehension
         * */
         if(idx < (this->totalNumLigand - this->numAnchorLigand)) {
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
         } else  {
             unsigned int anchor_idx = idx - (this->totalNumLigand - this->numAnchorLigand);
             unsigned int currRecId = 0;
             for(unsigned int& numAnchorLigForRec : this->numAnchorLigandPerReceptor)
             {
                 if(static_cast<long int>(anchor_idx) - static_cast<long int>(numAnchorLigForRec) < 0)
                 {
                     return std::make_tuple(currRecId, anchor_idx);
                 }
                 anchor_idx = anchor_idx - numAnchorLigForRec;
                 currRecId++;
             }
         }

        BOOST_LOG_TRIVIAL(error) << "Unexpected indexing error : unable to translate global index " << idx;
        BOOST_LOG_TRIVIAL(error) << "   Total num ligand : " << this->totalNumLigand;
        BOOST_LOG_TRIVIAL(error) << "   Anchor ligand : " << this->numAnchorLigand;
        BOOST_LOG_TRIVIAL(error) << "   Non-anchor (deduced) : " << (this->totalNumLigand - this->numAnchorLigand);
        BOOST_LOG_TRIVIAL(error) << "Internal state not sound, terminating.";
        std::terminate();
        return std::make_tuple(0, idx);
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

    bool MPICalibratorDirector::addAnchorLigandFromMol2File(ReceptorID recID, std::string& filename) {
        auto aLig = std::make_shared<Molecule>();
        aLig->populateFromMol2File(filename);

        this->anchorLigands[recID].push_back(aLig);
        this->numAnchorLigand++;
        this->totalNumLigand++;
        return true;
    }



    MPICalibratorDirector_ResumeObject MPICalibratorDirector::createResumeState() {
        MPICalibratorDirector_ResumeObject ro;
        ro.lossHistory = this->lossHistory;
        ro.durationHistory = this->durationHistory;
        ro.coefficientHistory = this->coefficientHistory;
        ro.gradientHistory = this->gradientHistory;


        ro.currentCoeffs = this->currentCoeffs;

        ro.coeffsToCalibrate = this->coeffsToCalibrate;
        ro.nameOfAllCoeffs = this->nameOfAllCoeffs;
        ro.idxOfCoeffsToCalibrate = this->idxOfCoeffsToCalibrate;

        ro.scoringFunctionType = this->scoringFunctionType;
        ro.heuristicType = this->heuristicType;
        ro.localOptimizerType = this->localOptimizerType;
        return ro;
    }


    bool MPICalibratorDirector::restoreResumeState(const MPICalibratorDirector_ResumeObject& state) {
        this->lossHistory = state.lossHistory;
        this->durationHistory = state.durationHistory;
        this->coefficientHistory = state.coefficientHistory;
        this->gradientHistory = state.gradientHistory;

        this->currentCoeffs = state.currentCoeffs;

        this->coeffsToCalibrate = state.coeffsToCalibrate;
        this->nameOfAllCoeffs = state.nameOfAllCoeffs;
        this->idxOfCoeffsToCalibrate = state.idxOfCoeffsToCalibrate;

        this->scoringFunctionType = state.scoringFunctionType;
        this->heuristicType = state.heuristicType;
        this->localOptimizerType = state.localOptimizerType;

        // Duplicate removal
        {
            std::set<int> s_idx;
            unsigned size = this->idxOfCoeffsToCalibrate.size();
            for( unsigned i = 0; i < size; ++i ) s_idx.insert( this->idxOfCoeffsToCalibrate[i] );
            this->idxOfCoeffsToCalibrate.assign( s_idx.begin(), s_idx.end() );
        }   

        {
            std::set<std::string> s_allCoeffsName;
            unsigned size = this->nameOfAllCoeffs.size();
            for( unsigned i = 0; i < size; ++i ) s_allCoeffsName.insert( this->nameOfAllCoeffs[i] );
            this->nameOfAllCoeffs.assign( s_allCoeffsName.begin(), s_allCoeffsName.end() );
        }   

        {
            std::set<std::string> s_coeffsToCalibrate;
            unsigned size = this->nameOfAllCoeffs.size();
            for( unsigned i = 0; i < size; ++i ) s_coeffsToCalibrate.insert( this->coeffsToCalibrate[i] );
            this->coeffsToCalibrate.assign( s_coeffsToCalibrate.begin(), s_coeffsToCalibrate.end() );
        }   

        // FIX for currentCoeffs containing the first coeffs actually
        this->currentCoeffs = this->coefficientHistory.back();


        {
            iConformer dummy_cf;
            dummy_cf.num_rotatable_bond = 0;
            iProtein dummy_prot;
            iTransform dummy_tr = iTransformIdentityInit(0);
            auto dummy_sf = scoringFunctionFactory(this->scoringFunctionType,
                                                    dummy_cf,
                                                    dummy_prot,
                                                    dummy_tr,
                                                    1e-3,
                                                    true);
            std::vector<std::string> names = dummy_sf->getCoefficientsNames();
            unsigned int num_coeffs = names.size();
            this->currentCoeffs.erase(this->currentCoeffs.begin() + num_coeffs, this->currentCoeffs.end());
        } 

        this->batchCount = this->coefficientHistory.size();

        return true;
    }


}
