//
// Created by eliane on 13/03/19.
//

#ifndef SMOLDOCK_REFERENCELIGANDHOLDER_H
#define SMOLDOCK_REFERENCELIGANDHOLDER_H


#include <string>
#include <vector>
#include <random>

#include <tbb/tbb.h>

#include <Structures/Molecule.h>
#include <Structures/Protein.h>
#include <Engines/AbstractDockingEngine.h>

#include <Engines/ScoringFunctions/ScoringFunctionInterface.h>
#include <Engines/ScoringFunctions/ScoringFunctionFactory.h>
#include <Engines/LocalOptimizers/OptimizerInterface.h>
#include <Engines/LocalOptimizers/OptimizerFactory.h>
#include <Engines/GlobalHeuristics/HeuristicInterface.h>
#include <Engines/GlobalHeuristics/HeuristicFactory.h>



namespace SmolDock::Calibration {




    struct CalibratorWorkItem {
        Score::ScoringFunctionType scFuncType_;
        Heuristics::GlobalHeuristicType heurType_;
        Optimizer::LocalOptimizerType localOptimizerType_;
        iTransform transform_;
        iConformer* conformer_ = nullptr;
        const iProtein* prot_ = nullptr;
        const iProtein* fullProt_ = nullptr;
        unsigned int seed_;
        unsigned int retryNumber_;
        Heuristics::HeuristicParameters hParams_;
        double referenceScore_;

    };

    class CalibratorLoopRunner {
    public:

        explicit CalibratorLoopRunner(std::shared_ptr<std::vector<CalibratorWorkItem>> workItemList,
                                      std::shared_ptr<std::mutex> resultMutex_,
                                    std::shared_ptr<std::vector<double>> local_scores_,
                                    std::shared_ptr<std::vector<double>> local_referenceScores_,
                                    std::shared_ptr<std::vector<std::vector<std::tuple<std::string,double>>>> local_scoreComponents_,
                                    std::vector<double> currentCoeffs_,
                                    std::vector<unsigned int> indexShufflingArray_ = {});

        void operator()(const tbb::blocked_range<size_t>& r) const;

    private:
        std::shared_ptr<std::vector<CalibratorWorkItem>> workItemList;

        std::shared_ptr<std::mutex> resultMutex;
        std::shared_ptr<std::vector<double>> local_scores;
        std::shared_ptr<std::vector<double>> local_referenceScores;
        std::shared_ptr<std::vector<std::vector<std::tuple<std::string,double>>>> local_scoreComponents;

        std::vector<double> currentCoeffs;

        std::vector<unsigned int> indexShufflingArray;

    };

    class CalibratorEnsmallenLayer {
    public:

        CalibratorEnsmallenLayer(std::shared_ptr<std::vector<CalibratorWorkItem>> workitemVector_,
                std::vector<double> startingCoeffs_,
                std::vector<unsigned int> idxOfCoeffsToCalibrate_,
                unsigned int seed_ = 312,
                double differentialEpsilon_ = 0.5);

        double Evaluate(const arma::mat& x, size_t i, size_t batchSize);

        void Gradient(const arma::mat& x, size_t i,arma::mat& g, size_t batchSize);

        void Shuffle();

        size_t NumFunctions();

        arma::mat getInitialParamMatrix();

    private:

        double doRealEvaluate(const arma::mat& x, size_t i, size_t batchSize);

        std::shared_ptr<std::vector<CalibratorWorkItem>> workitemVector;
        std::vector<double> startingCoeffs;
        std::vector<unsigned int> idxOfCoeffsToCalibrate;

        unsigned long numWorkItem;
        unsigned int paramSize;
        std::vector<unsigned int> shuffledIndexArray;

        double differentialEpsilon;
        std::mt19937 rndGenerator;



    };



    class Calibrator {
    public:
        Calibrator(Score::ScoringFunctionType scFuncType,
                   Heuristics::GlobalHeuristicType heurType,
                   Optimizer::LocalOptimizerType localOptimizerType_,
                   unsigned int maxLearningSteps = 10,
                   double initialLearningRate_ = 0.5,
                   unsigned int rngSeed= 374,
                   unsigned int conformerNumber = 4,
                   unsigned int retryNumber = 4,
                   unsigned int batchSize_ = 5,
                   Heuristics::HeuristicParameters hParams = Heuristics::emptyParameters);

        ~Calibrator() = default;

        using ReceptorID = unsigned int;

        ReceptorID addReceptor(const Protein& prot, Engine::AbstractDockingEngine::DockingBoxSetting dbsettings);

        bool addReferenceLigand_SMILES_Ki(ReceptorID recID, const std::string& smiles, double Ki, int seed = 364);
        bool addReferenceLigand_Mol_Ki(ReceptorID recID, const Molecule& mol, double Ki, int seed = 364);

        bool coefficientsToCalibrate(std::vector<std::string> nameOfCoeffs);

        bool setupCalibration();
        bool runCalibration();

        bool runCalibration2();



    private:

        void fillWorkItemVector(std::shared_ptr<std::vector<CalibratorWorkItem>> workItemVector);

        std::vector<std::string> coeffsToCalibrate;
        std::vector<std::string> nameOfAllCoeffs;
        std::vector<unsigned int> idxOfCoeffsToCalibrate;
        std::vector<double> currentCoeffs;

        Score::ScoringFunctionType scoringFunctionType;
        Heuristics::GlobalHeuristicType heuristicType;
        Optimizer::LocalOptimizerType localOptimizerType;

        std::shared_ptr<Score::ScoringFunction> dummy_sf;

        unsigned int maxLearningSteps;
        double initialLearningRate;


        std::mt19937 rndGenerator;

        unsigned int conformerNumber;
        unsigned int retryNumber;
        unsigned int batchSize;

        Heuristics::HeuristicParameters hParams;


        ReceptorID current_max_ReceptorID = 0;
        std::vector< std::tuple<std::shared_ptr<Protein>, Engine::AbstractDockingEngine::DockingBoxSetting, iProtein, iProtein> > referenceReceptor;
        std::map<ReceptorID ,std::vector<std::tuple<std::shared_ptr<Molecule>, double, std::vector<iConformer>>>> referenceLigands;


        arma::mat optResultMat;
        std::vector<double> optResultVector;

    };
}



#endif //SMOLDOCK_REFERENCELIGANDHOLDER_H
