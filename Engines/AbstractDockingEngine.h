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

#ifndef SMOLDOCK_ABSTRACTDOCKINGENGINE_H
#define SMOLDOCK_ABSTRACTDOCKINGENGINE_H

#include "../Structures/Protein.h"
#include "../Structures/Results/DockingResult.h"
#include "../Structures/Molecule.h"

/*! \namespace SmolDock Global namespace for SmolDock */
namespace SmolDock {

    /*! \namespace SmolDock::Engine Docking engines */
    namespace Engine {

        //! The expected interface for all docking engines.
        /*!
         *
         *
         */
        class AbstractDockingEngine {
        public:

            enum class DockingBoxSetting {
                everything,
                solventExposed
            };

            enum class DockingBoxShape {
                sphere,
                cube
            };

            // /// Parameters /////////////

            //! Load a protein for docking.
            /*!
             * Load a protein for docking. Note that it *can* be modified by the engine.
             *
             * \param p Pointer to the Protein
             * \return Returns whether it has been successful
             * \sa setLigand(), setDockingBox()
            */
            virtual bool setProtein(Protein *p) = 0;

            //! Load a molecule ligand for docking.
            /*!
             * Load a molecule ligand for docking. Note that it *can* be modified by the engine.
             *
             * \param m Pointer to the Molecule
             * \return Returns whether it has been successful
             * \sa setProtein(), setDockingBox()
            */
            virtual bool setLigand(Molecule *m) = 0;

            //! Set the protein domain/box to consider for docking.
            /*! A restriction of the search space can be necessary : this provides a way to specify a domain/box in the protein
             * where the docking will take place.
             *
             * \param setting General-scope DockingBox specification
             * \return Returns whether it has been successful
             * \sa setProtein(), setLigand()
            */
            virtual bool setDockingBox(DockingBoxSetting setting) = 0;
            // virtual bool setDockingBox(std::tuple<double,double,double> center, double r, DockingBoxShape shape = DockingBoxShape::sphere) = 0;

            //! Set the RNG seed for reproducibility.
            /*!
             * \param seed RNG see
             * \return Returns whether it has been successful
             * \sa setProtein(), setLigand(), setDockingBox()
            */
            virtual void setRandomSeed(int seed) = 0;

            // /// Actions /////////////

            //! Prepare the engine for docking.
            /*!
             * Prepare the engine for docking. It is expected in future versions to separate protein, ligand, and prot-ligand
             * initilization, so as to allow maximum work sharing in screening scenario (same prot = same prot prep, just
             * different ligand)
             *
             * \return Returns whether it has been successful
             * \sa
            */
            virtual bool setupDockingEngine() = 0;
            // virtual bool setupLigand() = 0;
            // virtual bool setupProtein() = 0;
            // virtual bool setupProteinLigand() = 0;
            // virtual ??? getExportableLigandPrep() =0;
            // virtual ??? getExportableProteinPrep() =0;


            //! Run the docking
            /*!
             *
             * \return Returns whether it has been successful
             * \sa
            */
            virtual void runDockingEngine() = 0;

            // /// Exports /////////////


            // /// Results /////////////

            //! Get the docking results
            /*!
             * Behaviour is undefined if called before runDockingEngine()
             *
             * \return The results of the docking.
             * \sa runDockingEngine()
            */
            virtual std::shared_ptr<DockingResult> getDockingResult() = 0;
        };

    }

}

#endif //SMOLDOCK_ABSTRACTDOCKINGENGINE_H
