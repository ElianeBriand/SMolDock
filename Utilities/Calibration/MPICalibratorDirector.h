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

namespace mpi = boost::mpi;

#include "Calibrator.h"
#include "MPICalibratorCommon.h"
#include <Structures/Atom.h>

namespace SmolDock::Calibration {


    class MPICalibratorDirector : public Calibrator {

    public:

        MPICalibratorDirector(mpi::environment &env_,
                              mpi::communicator &world_,
                              Score::ScoringFunctionType scFuncType,
                              Heuristics::GlobalHeuristicType heurType,
                              Optimizer::LocalOptimizerType localOptimizerType_,
                              unsigned int maxLearningSteps = 1000,
                              double initialLearningRate_ = 0.5,
                              unsigned int rngSeed = 374,
                              unsigned int conformerNumber = 4,
                              unsigned int retryNumber = 4,
                              unsigned int batchSize_ = 5,
                              Heuristics::HeuristicParameters hParams = Heuristics::emptyParameters);

        virtual bool setupCalibration();
        virtual bool runCalibration();

        ReceptorID addReceptorFromFile(const std::string& filename, Engine::AbstractDockingEngine::DockingBoxSetting dbsettings);
        virtual ReceptorID addReceptor(const Protein& prot, Engine::AbstractDockingEngine::DockingBoxSetting dbsettings);

        bool applySpecialResidueTypingFromRecID(ReceptorID recID, const AminoAcid::AAType resType,
                                                const unsigned int serialNumber,
                                                const SpecialResidueTyping specialType);

        virtual bool addReferenceLigand_SMILES_Ki(ReceptorID recID, const std::string& smiles, double Ki, int seed = 364);
        virtual bool addReferenceLigand_Mol_Ki(ReceptorID recID, const Molecule& mol, double Ki, int seed = 364);

        virtual bool addAnchorLigandFromMol2File(ReceptorID recID, std::string& filename);

        bool applyVariantToAllLigands(const std::string& SMARTSPattern, Atom::AtomVariant variant);

        double EvaluateWithGradient(const arma::mat& x,
                                            const size_t i,
                                            arma::mat& g,
                                            const size_t batchSize);

        void Shuffle();

        size_t NumFunctions();


    private:

        std::tuple<unsigned int, unsigned int> RecLigIdxFromGlobalIdx(unsigned int idx);

        void updateAndPrintStatus();


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
        std::vector<std::vector<double>> coefficientHistory;
        std::vector<std::vector<double>> gradientHistory;



    };

}

#endif //SMOLDOCK_MPICALIBRATORDIRECTOR_H
