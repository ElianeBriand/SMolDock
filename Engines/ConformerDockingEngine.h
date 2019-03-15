/*
 * Copyright (c) 2018 Eliane Briand
 *
 * This file is part of SmolDock.
 *
 * SmolDock is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SmolDock is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SmolDock.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef SMOLDOCK_CONFORMERRIGIDDOCKINGENGINE_H
#define SMOLDOCK_CONFORMERRIGIDDOCKINGENGINE_H

#include <GraphMol/RWMol.h>
#include <random>
#include <Engines/LocalOptimizers/OptimizerFactory.h>
#include <Engines/GlobalHeuristics/HeuristicFactory.h>


#include "AbstractDockingEngine.h"
#include "../Structures/Protein.h"
#include "../Structures/Molecule.h"

#include "Engines/LocalOptimizers/GradientDescentLineSearch.h"

#include "Engines/ScoringFunctions/ScoringFunctionFactory.h"

namespace SmolDock::Engine {

    //! Generates a lot of conformer, then runs rigid docking on them.
    /*!
     * The general idea is generating lots of conformer, then docking them with constant bond angle
     * and distance ("rigid" docking). The theory being that this may be faster than real docking
     * and still produce a similar affinity ranking. (In practice, I do not know if it would work :
     * hence this to try it
     */
    class ConformerDockingEngine : public AbstractDockingEngine {

    public:

        //!
        /*!
         * \param conformer_num Number of conformer to generate
        */
        explicit ConformerDockingEngine(unsigned int conformer_num,
                                             unsigned int retryPerConformer,
                                             Protein* protein,
                                             Molecule* ligand,
                                             Score::ScoringFunctionType scFuncType,
                                             Heuristics::GlobalHeuristicType heurType,
                                             Optimizer::LocalOptimizerType localOptimizerType_,
                                             unsigned int seed,
                                             Heuristics::HeuristicParameters hParams = Heuristics::emptyParameters,
                                             bool rigidDocking = false);

        // /// Parameters /////////////
        bool setDockingBox(DockingBoxSetting setting) final;


        // /// Actions /////////////
        bool setupDockingEngine() final;

        void runDockingEngine() final;

        // /// Results /////////////
        std::shared_ptr<DockingResult> getDockingResult() final;

        std::tuple<double,double> getMeanStdDevDuration() const;
        std::tuple<double,double> getMeanStdDevScore() const;
        double getBestScore();

        virtual ~ConformerDockingEngine() = default;

    private:



        unsigned int conformer_num;
        unsigned int retryPerConformer;

        Protein* orig_protein;
        Molecule* orig_ligand;

        Score::ScoringFunctionType scoringFuncType;
        Heuristics::GlobalHeuristicType heuristicType;
        Optimizer::LocalOptimizerType localOptimizerType;

        //std::shared_ptr<Score::ScoringFunction> scoringFunction;
        //std::shared_ptr<Heuristics::GlobalHeuristic> globalHeuristic;
        //std::shared_ptr<Optimizer::Optimizer> localOptimizer;

        std::mt19937 rnd_generator;


        Heuristics::HeuristicParameters heurParams;

        bool rigidDocking = false;

        //std::shared_ptr<RDKit::RWMol> rwmol;

        std::vector<iConformer> viConformers;

        iProtein protein;

        iProtein fullProtein;

        DockingBoxSetting dockBoxSettings;


        std::vector<double> scores;
        std::vector<double> localScores;
        std::vector<double> startingScores;

        std::vector<iConformer> allGeneratediConformer;
        std::vector<iConformer> bestiConformer;

        double meanDuration = 0.0;
        double stdDevDuration = 0.0;

        double meanScore = 0.0;
        double stdDevScore = 0.0;
        double bestScore = 0.0;

    };

}


#endif //SMOLDOCK_CONFORMERRIGIDDOCKINGENGINE_H
