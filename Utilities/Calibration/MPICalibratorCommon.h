//
// Created by briand on 3/21/19.
//

#ifndef SMOLDOCK_MPICALIBRATORCOMMON_H
#define SMOLDOCK_MPICALIBRATORCOMMON_H

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

#include <Engines/AbstractDockingEngine.h>

namespace SmolDock::Calibration {

    enum MTags {
        RegisterNode,
        ReadyForWork,
        TaskRequest,
        ResultMessage,
        NodeTerminating
    };

    struct RegisterMessage {
        unsigned int rank;
        unsigned int numThread;


        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & rank;
            ar & numThread;
        }

    };

    struct WorkStructure {
        unsigned int numReceptor;
        unsigned int conformerPerLigand;
        unsigned int retryPerConformer;
        unsigned int numRegisteredVariant;
        unsigned int totalNumLigand;
        int seed;

        unsigned int numCoeffToCalibrate;
        std::vector<unsigned int> idxOfCoeffsToCalibrate;

        Score::ScoringFunctionType scoringFunctionType;
        Heuristics::GlobalHeuristicType heuristicType;
        Optimizer::LocalOptimizerType localOptimizerType;

        double differentialEpsilon;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & numReceptor;
            ar & conformerPerLigand;
            ar & retryPerConformer;
            ar & numRegisteredVariant;
            ar & totalNumLigand;
            ar & seed;

            ar & numCoeffToCalibrate;
            ar & idxOfCoeffsToCalibrate;

            ar & scoringFunctionType;
            ar & heuristicType;
            ar & localOptimizerType;

            ar & differentialEpsilon;
        }

    };

    struct RegisteredVariant {
        std::string smarts;
        unsigned int atomVariant;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & smarts;
            ar & atomVariant;
        }

    };

    struct MPISpecialResidueTyping {

        unsigned int aaType;
        unsigned int serialNumber;
        unsigned int specialTyping;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & aaType;
            ar & serialNumber;
            ar & specialTyping;
        }
    };


    struct ReceptorRecord {
        std::string PDBBlock;
        Engine::AbstractDockingEngine::DockingBoxSetting dbsetting;
        std::vector<MPISpecialResidueTyping> specialResTypes;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & PDBBlock;
            ar & dbsetting;
            ar & specialResTypes;
        }

    };

    struct LigandRecord {
        unsigned int receptorID;
        std::string smiles;
        double deltaG;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & receptorID;
            ar & smiles;
            ar & deltaG;
        }

    };



    struct Task {
        unsigned int receptorID;
        unsigned int ligandIdx;
        std::vector<double> coefficients;

        bool endOfTasks = false;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & receptorID;
            ar & ligandIdx;
            ar & coefficients;
            ar & endOfTasks;
        }


    };

    struct Result {
        double loss;
        std::vector<double> lossGradient;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & loss;
            ar & lossGradient;
        }


    };




}

namespace boost { namespace mpi {
        template <>
        struct is_mpi_datatype<SmolDock::Calibration::RegisterMessage> : mpl::true_ { };
    } }


#endif //SMOLDOCK_MPICALIBRATORCOMMON_H
