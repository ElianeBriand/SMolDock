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

#ifndef SMOLDOCK_VINACOMPATIBLEDOCKINGENGINE_H
#define SMOLDOCK_VINACOMPATIBLEDOCKINGENGINE_H


#include "AbstractDockingEngine.h"


#include <rdkit/GraphMol/RWMol.h>
#include <random>
#include <Engines/LocalOptimizers/OptimizerFactory.h>
#include <Engines/GlobalHeuristics/HeuristicFactory.h>


#include "AbstractDockingEngine.h"
#include "../Structures/Protein.h"
#include "../Structures/Molecule.h"

#include "Engines/LocalOptimizers/GradientDescentLineSearch.h"

#include "Engines/ScoringFunctions/ScoringFunctionFactory.h"

#include "ConformerDockingEngine.h"

namespace SmolDock::Engine {

    class VinaCompatibleDockingEngine : AbstractDockingEngine {

    public:

        VinaCompatibleDockingEngine(Protein* protein,
                                    Molecule* ligand,
                                    unsigned int seed);



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

        virtual ~VinaCompatibleDockingEngine() = default;


    private:

        std::shared_ptr<ConformerDockingEngine> internalEngine;

    };

}


#endif //SMOLDOCK_VINACOMPATIBLEDOCKINGENGINE_H
