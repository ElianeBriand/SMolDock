//
// Created by briand on 3/21/19.
//

#ifndef SMOLDOCK_MPICALIBRATORNODE_H
#define SMOLDOCK_MPICALIBRATORNODE_H

#include <thread>
#include <memory>

#include <string>
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>

namespace mpi = boost::mpi;

#include "Calibrator.h"
#include "MPICalibratorCommon.h"



namespace SmolDock::Calibration {


    class MPICalibrator2DLoopRunner {
    public:

        MPICalibrator2DLoopRunner(const WorkStructure& ws,
                const std::vector<iConformer>& conformerVector_,
                const iProtein& protein_,
                const iProtein& fullProtein_,
                std::shared_ptr<std::mutex> resultMutex_,
                std::shared_ptr<std::vector<double>> local_scores_,
                std::shared_ptr<std::vector<iConformer>> local_iconformers_,
                std::vector<double> currentCoeffs_);

        void operator()(const tbb::blocked_range2d<size_t>& r) const;

    private:
        const WorkStructure& workStructure;
        const std::vector<iConformer>& conformerVector;
        const iProtein& protein;
        const iProtein& fullProtein;
        std::shared_ptr<std::mutex> resultMutex;
        std::shared_ptr<std::vector<double>> local_scores;
        std::shared_ptr<std::vector<iConformer>> local_iconformers;
        std::vector<double> currentCoeffs;

    };

    class MPICalibratorNode {

    public:
        MPICalibratorNode(mpi::environment& env_,
                          mpi::communicator& world_);

        void runNode();


    private:


        mpi::environment& env;
        mpi::communicator& world;

        std::mt19937 rndGenerator;

        int rank;

        WorkStructure workStructure;
        std::vector<std::string> pdbBlockStrings;
        std::vector<Engine::AbstractDockingEngine::DockingBoxSetting> dbSettings;
        std::vector<std::vector<MPISpecialResidueTyping>> specialResidueTypings;

        std::map<Calibrator::ReceptorID, std::vector<std::tuple<std::string, double, std::shared_ptr<Molecule>, std::vector<iConformer>>>>
                ligandSmilesRefScore;

        std::map<Calibrator::ReceptorID, std::vector<std::tuple<std::string, std::shared_ptr<Molecule>,std::shared_ptr<Molecule>, std::vector<iConformer>>>>
                anchorLigands;

        std::vector<std::tuple<std::shared_ptr<Protein>, Engine::AbstractDockingEngine::DockingBoxSetting, iProtein, iProtein> > referenceReceptor;


        std::vector<RegisteredVariant> variantsForAllLigands;

    };

}
#endif //SMOLDOCK_MPICALIBRATORNODE_H
