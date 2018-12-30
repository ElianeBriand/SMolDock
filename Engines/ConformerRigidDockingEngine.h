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


#include "AbstractDockingEngine.h"
#include "../Structures/Protein.h"
#include "../Structures/Molecule.h"

#include "Engines/Optimizers/GradientDescentLineSearch.h"


namespace SmolDock::Engine {

    //! Generates a lot of conformer, then runs rigid docking on them.
    /*!
     * The general idea is generating lots of conformer, then docking them with constant bond angle
     * and distance ("rigid" docking). The theory being that this may be faster than real docking
     * and still produce a similar affinity ranking. (In practice, I do not know if it would work :
     * hence this to try it
     */
    class ConformerRigidDockingEngine : public AbstractDockingEngine {

    public:

        //!
        /*!
         * \param conformer_num Number of conformer to generate
        */
        explicit ConformerRigidDockingEngine(unsigned int conformer_num);

        // /// Parameters /////////////
        bool setDockingBox(DockingBoxSetting setting) final;

        bool setProtein(Protein *p) final;

        bool setLigand(Molecule *m) final;

        void setRandomSeed(int seed) final;


        // /// Actions /////////////
        bool setupDockingEngine() final;

        void runDockingEngine() final;

        // /// Results /////////////
        std::shared_ptr<DockingResult> getDockingResult() final;


    private:

        unsigned int conformer_num;

        Protein *orig_protein;
        Molecule *orig_ligand;

        int random_seed = 1;
        std::mt19937 rnd_generator;

        //std::shared_ptr<RDKit::RWMol> rwmol;

        std::vector<iConformer> viConformers;

        iProtein protein;


        std::vector<double> scores;
        std::vector<iConformer> final_iConformer;
    };

}


#endif //SMOLDOCK_CONFORMERRIGIDDOCKINGENGINE_H
