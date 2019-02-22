//
// Created by eliane on 09/02/19.
//

#include <Engines/Internals/InternalsUtilityFunctions.h>
#include "PoseRefiner.h"


namespace SmolDock::Engine {


    PoseRefiner::PoseRefiner(Protein *protein, Molecule *ligand, Score::ScoringFunctionType scFuncType,
                             Optimizer::LocalOptimizerType localOptimizerType_,
                             unsigned int seed) :
            orig_protein(protein),
            orig_ligand(ligand),
            scoringFuncType(scFuncType),
            localOptimizerType(localOptimizerType_),
            rnd_generator(seed) {

    }

    bool PoseRefiner::refinePose() {

        std::uniform_int_distribution<unsigned int> dis_uint(0, std::numeric_limits<unsigned int>::max());

        iConformer starting_conformer = this->orig_ligand->getInitialConformer(true);
        iTransform starting_tr = iTransformIdentityInit();
        starting_tr.transl = starting_conformer.centroidNormalizingTransform;

        this->protein = this->orig_protein->getiProtein();

        this->scoringFunction = scoringFunctionFactory(this->scoringFuncType,
                                                       starting_conformer,
                                                       this->protein,
                                                       starting_tr,
                                                       1e-5);

        this->localOptimizer = optimizerFactory(this->localOptimizerType,
                                                this->scoringFunction.get(),
                                                1e-5);


        Heuristics::HeuristicParameters hParams = Heuristics::heuristicParametersFactory(
                Heuristics::GlobalHeuristicType::OnlyLocal);

        this->globalHeuristic = globalHeuristicFactory(Heuristics::GlobalHeuristicType::OnlyLocal,
                                                       this->scoringFunction.get(),
                                                       this->localOptimizer.get(),
                                                       dis_uint(this->rnd_generator),
                                                       hParams);


        arma::mat startingCond = this->scoringFunction->getStartingConditions();
        this->initial_score = this->scoringFunction->Evaluate(startingCond);


        this->globalHeuristic->search();

        arma::mat rawResultMatrix = this->globalHeuristic->getResultMatrix();

        this->final_score = this->scoringFunction->Evaluate(rawResultMatrix);
        this->result = this->scoringFunction->getConformerForParamMatrix(rawResultMatrix);


        return true;
    }

    double PoseRefiner::getInitialScore() {
        return this->initial_score;
    }

    double PoseRefiner::getFinalScore() {
        return this->final_score;
    }

    double PoseRefiner::getScoreDifference() {
        return (this->initial_score - this->final_score);
    }

    bool PoseRefiner::applyToLigand() {
        this->orig_ligand->updateAtomPositionsFromiConformer(this->result);
        return true;
    }
}