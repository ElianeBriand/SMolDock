//
// Created by briand on 3/21/19.
//

#ifndef SMOLDOCK_MPICALIBRATORDIRECTOR_H
#define SMOLDOCK_MPICALIBRATORDIRECTOR_H


#include <string>
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


        virtual bool addReferenceLigand_SMILES_Ki(ReceptorID recID, const std::string& smiles, double Ki, int seed = 364);
        virtual bool addReferenceLigand_Mol_Ki(ReceptorID recID, const Molecule& mol, double Ki, int seed = 364);

        bool applyVariantToAllLigands(const std::string& SMARTSPattern, Atom::AtomVariant variant);

        double EvaluateWithGradient(const arma::mat& x,
                                            const size_t i,
                                            arma::mat& g,
                                            const size_t batchSize);

        void Shuffle();

        size_t NumFunctions();


    private:

        std::tuple<unsigned int, unsigned int> RecLigIdxFromGlobalIdx(unsigned int idx);


        mpi::environment& env;
        mpi::communicator& world;

        std::vector<RegisterMessage> nodeInfo;

        unsigned int numProcess;
        unsigned int numTotalCPU;
        unsigned int numLigand = 0;

        std::vector<std::string> pdbBlockStrings;
        std::vector<Engine::AbstractDockingEngine::DockingBoxSetting> dbSettings;

        std::map<Calibrator::ReceptorID ,std::vector<std::tuple<std::string,double>>> ligandSmilesRefScore;
        std::vector<unsigned int> numLigandPerReceptor;
        std::vector<unsigned int> indexShuffler;

        std::vector<RegisteredVariant> variantsForAllLigands;

        WorkStructure workStructure;



    };

}

#endif //SMOLDOCK_MPICALIBRATORDIRECTOR_H
