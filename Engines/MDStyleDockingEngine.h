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

#ifndef SMOLDOCK_MDSTYLEDOCKINGENGINE_H
#define SMOLDOCK_MDSTYLEDOCKINGENGINE_H

#include "AbstractDockingEngine.h"
#include "../Structures/Protein.h"
#include "../Structures/Molecule.h"

namespace SmolDock::Engine {

    class MDStyleDockingEngine : AbstractDockingEngine {


    public:
        MDStyleDockingEngine() = default;


        bool setDockingBox(DockingBoxSetting setting) final;

        bool setupDockingEngine() final;

        void runDockingEngine() final;

        std::shared_ptr<DockingResult> getDockingResult() final;


    private:
/*
        Protein *protein = nullptr;
        Molecule *ligand = nullptr;
*/
    };

}

#endif //SMOLDOCK_MDSTYLEDOCKINGENGINE_H
