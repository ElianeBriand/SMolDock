//
// Created by eliane on 09/02/19.
//

#ifndef SMOLDOCK_VINAREFINER_H
#define SMOLDOCK_VINAREFINER_H

#include <Structures/Protein.h>
#include <Structures/Molecule.h>
#include <Engines/ScoringFunctions/ScoringFunctionFactory.h>
#include <Engines/GlobalHeuristics/HeuristicFactory.h>
#include <Engines/LocalOptimizers/OptimizerFactory.h>

namespace SmolDock::Engine {

    class PoseRefiner {
    public:
        PoseRefiner(Protein* protein,
                    Molecule* ligand,
                    Score::ScoringFunctionType scFuncType,
                    Optimizer::LocalOptimizerType localOptimizerType_,
                    unsigned int seed);

        bool refinePose();

        double getInitialScore();

        double getFinalScore();

        double getScoreDifference();

        bool applyToLigand();

    private:
        Protein* orig_protein;
        Molecule* orig_ligand;

        Score::ScoringFunctionType scoringFuncType;
        Optimizer::LocalOptimizerType localOptimizerType;

        std::shared_ptr<Score::ScoringFunction> scoringFunction;
        std::shared_ptr<Heuristics::GlobalHeuristic> globalHeuristic;
        std::shared_ptr<Optimizer::Optimizer> localOptimizer;

        std::mt19937 rnd_generator;

        iProtein protein;

        double initial_score = 0.0;
        double final_score = 0.0;

        iConformer result;
    };

}


#endif //SMOLDOCK_VINAREFINER_H
