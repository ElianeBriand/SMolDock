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
 * Foobar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "MDStyleDockingEngine.h"


namespace SmolDock::Engine {


    bool MDStyleDockingEngine::setProtein(Protein *p) {
        this->protein = p;
        return true;
    }

    bool MDStyleDockingEngine::setLigand(Molecule *m) {
        this->ligand = m;
        return true;
    }

    bool MDStyleDockingEngine::setupDockingEngine() {
        return true;
    }

    void MDStyleDockingEngine::runDockingEngine() {

    }

    std::shared_ptr<DockingResult> MDStyleDockingEngine::getDockingResult() {
        return std::make_shared<DockingResult>();
    }


    bool MDStyleDockingEngine::setDockingBox(AbstractDockingEngine::DockingBoxSetting setting) {
        if (setting != DockingBoxSetting::everything) {
            std::cout << "[!] DockingBoxSetting (that is not DockingBoxSetting::everything) is not yet implemented."
                      << std::endl;
            std::cout << "[ ] Running as if DockingBoxSetting::everything was passed" << std::endl;
            return false;
        }
        return true;
    }
}