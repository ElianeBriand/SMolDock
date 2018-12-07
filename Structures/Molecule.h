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
#ifndef SMOLDOCK_MOLECULE_H
#define SMOLDOCK_MOLECULE_H

#include <vector>


#include <GraphMol/RDKitBase.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
/*
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionRunner.h>
#include <GraphMol/ChemReactions/PreprocessRxn.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>
 */
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>

#include "Structure.h"
#include "Atom.h"
#include "Bond.h"
#include "../Engines/Internals/iConformer.h"

namespace SmolDock {

    /*!
     * \brief Class representing "rich" molecule. Contains all extended attributes and functions.
     *
     */
    class Molecule : Structure {
        friend class MoleculeTraversal;

        friend class UnitTestHelper;

    public:
        Molecule();


        std::shared_ptr<RDKit::RWMol> getInternalRWMol() const;

        //! Populate atoms and bonds from SMILES.
        /*!
         * Populate atoms and bonds from SMILES. This uses the RDKit backend.
         * \param smiles SMILES string
         * \param seed RNG seed
         * \return whether the parsing and initial conformer generation was successful.
         * \sa generateConformers()
        */
        bool populateFromSMILES(const std::string &smiles, unsigned int seed = 367454);

        unsigned int numberOfAtoms();

        unsigned int numberOfBonds();

        //! Generate a conformer of this molecule
        /*!
         * \param conformer Reference to conformer object to be filled. Not modified if generation fails.
         * \param seed RNG seed
         * \return Returns true if successfully generated conformer, false otherwise
         * \sa generateConformers()
        */
        bool generateConformer(iConformer &conformer, int seed = 367454);

        //! Generate a vector of conformer of this molecule
        /*!
         * \param viConformers Pointer to vector to be filled with conformer. May already contain something.
         * \param num Number of conformers to attempt to generate.
         * \param seed RNG seed
         * \return Returns the number of conformer actually generated. May be 0.
         * \sa generateConformer()
        */
        unsigned int generateConformers(std::vector<iConformer>& viConformers, unsigned int num, int seed = 367454);


    private:
        std::vector<std::shared_ptr<Atom> > atoms;
        std::vector<std::shared_ptr<Bond> > bonds;

        std::shared_ptr<RDKit::RWMol> rwmol;
        std::string smiles;

    };

}

#endif //SMOLDOCK_MOLECULE_H
