//
// Created by eliane on 25/11/18.
//

#ifndef SMOLDOCK_ABSTRACTDOCKINGENGINE_H
#define SMOLDOCK_ABSTRACTDOCKINGENGINE_H

#include "../Structures/Protein.h"
#include "../Structures/Results/DockingResult.h"
#include "../Structures/Molecule.h"

namespace SmolDock {

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
