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

#include "VinaCompatibleDockingEngine.h"

namespace SmolDock::Engine {


    VinaCompatibleDockingEngine::VinaCompatibleDockingEngine(Protein* protein,
                                                             Molecule* ligand,
                                                             unsigned int seed) {
        this->internalEngine = std::make_shared<ConformerDockingEngine>(10,2,protein,ligand,Score::ScoringFunctionType::Vina,
                                          Heuristics::GlobalHeuristicType::SimulatedAnnealing,Optimizer::LocalOptimizerType::L_BFGS, seed);
    }

    bool VinaCompatibleDockingEngine::setDockingBox(AbstractDockingEngine::DockingBoxSetting setting) {
        return this->internalEngine->setDockingBox(setting);
    }

    bool VinaCompatibleDockingEngine::setupDockingEngine() {
        return this->internalEngine->setupDockingEngine();
    }

    void VinaCompatibleDockingEngine::runDockingEngine() {
        this->internalEngine->runDockingEngine();
    }

    std::shared_ptr<DockingResult> VinaCompatibleDockingEngine::getDockingResult() {
        return this->internalEngine->getDockingResult();
    }

    std::tuple<double, double> VinaCompatibleDockingEngine::getMeanStdDevDuration() const {
        return this->internalEngine->getMeanStdDevDuration();
    }

    std::tuple<double, double> VinaCompatibleDockingEngine::getMeanStdDevScore() const {
        return this->internalEngine->getMeanStdDevScore();
    }

    double VinaCompatibleDockingEngine::getBestScore() {
        return this->internalEngine->getBestScore();
    }
}