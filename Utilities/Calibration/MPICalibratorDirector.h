//
// Created by briand on 3/21/19.
//

#ifndef SMOLDOCK_MPICALIBRATORDIRECTOR_H
#define SMOLDOCK_MPICALIBRATORDIRECTOR_H


#include <string>
#include <atomic>
#include <thread>

#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

namespace mpi = boost::mpi;

#include "Calibrator.h"
#include "MPICalibratorCommon.h"
#include <Structures/Atom.h>

namespace SmolDock::Calibration {


    class MPICalibratorDirector_ResumeObject {
        // NB: update the MPICalibratorDirector resume methods when changing this
    public:

        std::vector<double> lossHistory;
        std::vector<int> durationHistory;
        std::vector<std::vector<double>> coefficientHistory;
        std::vector<std::vector<double>> gradientHistory;

        std::vector<double> currentCoeffs;

        std::vector<std::string> coeffsToCalibrate;
        std::vector<std::string> nameOfAllCoeffs;
        std::vector<unsigned int> idxOfCoeffsToCalibrate;

        Score::ScoringFunctionType scoringFunctionType;
        Heuristics::GlobalHeuristicType heuristicType;
        Optimizer::LocalOptimizerType localOptimizerType;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & lossHistory;
            ar & durationHistory;
            ar & coefficientHistory;
            ar & gradientHistory;

            ar & currentCoeffs;

            ar & coeffsToCalibrate;
            ar & nameOfAllCoeffs;
            ar & idxOfCoeffsToCalibrate;


            ar & scoringFunctionType;
            ar & heuristicType;
            ar & localOptimizerType;

        }
    };


    class MPICalibratorDirector : public Calibrator {

    public:

        MPICalibratorDirector(mpi::environment &env_,
                              mpi::communicator &world_,
                              Score::ScoringFunctionType scFuncType,
                              Heuristics::GlobalHeuristicType heurType,
                              Optimizer::LocalOptimizerType localOptimizerType_,
                              unsigned int maxLearningSteps = 1000,
                              double stepSize_ = 1e-7,
                              unsigned int rngSeed = 374,
                              unsigned int conformerNumber = 4,
                              unsigned int retryNumber = 4,
                              unsigned int batchSize_ = 5,
                              Heuristics::HeuristicParameters hParams = Heuristics::emptyParameters,
                              const std::string& restoreArchivePrefix_ = "mpicalibrator_out");

        virtual bool setupCalibration();
        virtual bool runCalibration();

        ReceptorID addReceptorFromFile(const std::string& filename, Engine::AbstractDockingEngine::DockingBoxSetting dbsettings);
        virtual ReceptorID addReceptor(const Protein& prot, Engine::AbstractDockingEngine::DockingBoxSetting dbsettings);

        bool applySpecialResidueTypingFromRecID(ReceptorID recID, const AminoAcid::AAType resType,
                                                const unsigned int serialNumber,
                                                const SpecialResidueTyping specialType);

        virtual bool addReferenceLigand_SMILES_Ki(ReceptorID recID, const std::string& smiles, double Ki);
        virtual bool addReferenceLigand_SMILES_deltaG(ReceptorID recID,const std::string& smiles, double deltaG);

        virtual bool addReferenceLigand_Mol_Ki(ReceptorID recID, const Molecule& mol, double Ki);

        virtual bool addAnchorLigandFromMol2File(ReceptorID recID, std::string& filename);

        bool applyVariantToAllLigands(const std::string& SMARTSPattern, Atom::AtomVariant variant);

        double EvaluateWithGradient(const arma::mat& x,
                                            const size_t i,
                                            arma::mat& g,
                                            const size_t batchSize);

        void Shuffle();

        size_t NumFunctions();


        MPICalibratorDirector_ResumeObject createResumeState();
        bool restoreResumeState(const MPICalibratorDirector_ResumeObject& state);

    private:

        std::tuple<unsigned int, unsigned int> RecLigIdxFromGlobalIdx(unsigned int idx);

        void updateAndPrintStatus();

        unsigned int batchCount = 0;

        mpi::environment& env;
        mpi::communicator& world;

        std::vector<RegisterMessage> nodeInfo;
        std::vector<unsigned int> numWorkItemPerNode;

        unsigned int numProcess;
        unsigned int numTotalCPU;
        unsigned int totalNumLigand = 0;
        unsigned int numAnchorLigand = 0;

        std::vector<std::string> pdbBlockStrings;
        std::vector<Engine::AbstractDockingEngine::DockingBoxSetting> dbSettings;
        std::vector<std::vector<MPISpecialResidueTyping>> specialResidueTypings;

        std::map<Calibrator::ReceptorID , std::vector<std::shared_ptr<Molecule>>> anchorLigands;
        std::vector<unsigned int> numAnchorLigandPerReceptor;

        std::map<Calibrator::ReceptorID ,std::vector<std::tuple<std::string,double>>> ligandSmilesRefScore;
        std::vector<unsigned int> numLigandPerReceptor;
        std::vector<unsigned int> indexShuffler;

        std::vector<RegisteredVariant> variantsForAllLigands;

        WorkStructure workStructure;

        std::atomic<bool> calibrationStillRunning;

        std::vector<double> lossHistory;
        std::vector<int> durationHistory;
        std::vector<std::vector<double>> coefficientHistory;
        std::vector<std::vector<double>> gradientHistory;

        std::string restoreArchivePrefix;



    };

}

#endif //SMOLDOCK_MPICALIBRATORDIRECTOR_H
