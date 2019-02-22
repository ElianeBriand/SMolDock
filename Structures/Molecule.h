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
#ifndef SMOLDOCK_MOLECULE_H
#define SMOLDOCK_MOLECULE_H

#include <vector>


#include <GraphMol/RDKitBase.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/FileParsers/FileParsers.h>

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
#include <Structures/InputPostProcessors/InputPostProcessorInterface.h>

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

        friend class PDBWriter;


    public:
        Molecule();



        //! Populate atoms and bonds from SMILES.
        /*!
         * Populate atoms and bonds from SMILES. This uses the RDKit backend.
         * \param smiles SMILES string
         * \param seed RNG seed
         * \return whether the parsing and initial conformer generation was successful.
         * \sa populateFromPDB()
        */
        bool populateFromSMILES(const std::string &smiles, unsigned int seed = 36754,
                                std::vector<std::shared_ptr<InputPostProcessor::InputPostProcessor> > postProcessors = {});


        bool populateFromMolFile(const std::string &filename, unsigned int seed = 36754,
                                 std::vector<std::shared_ptr<InputPostProcessor::InputPostProcessor> > postProcessors = {});

        bool populateFromMol2File(const std::string &filename, unsigned int seed = 36754,
                                  std::vector<std::shared_ptr<InputPostProcessor::InputPostProcessor> > postProcessors = {});


        std::string writeToMolBlock();

        bool writeToMolFile(const std::string &filename, bool overwrite = false);

        //! Populate atoms and bonds from a PDB file. Additionally, a SMILES string can be used for bond order hint.
        /*!
         *
         * Due to limitation of the PDB file format, it is difficult to infer bond order (single/double/aromatic bond) from
         * connectivity data. The underlying implementation in RDKit does not attempts it. A SMILES string with the correct
         * aromaticity/bonds can be supplied to help obtain accurate internal representation.
         *
         * Only the first model and first chain of the PDB file is used to construct the molecule.
         *
         * Populate atoms and bonds from SMILES. This uses the RDKit backend.
         * \param filename Path to PDB file
         * \param smiles_hint SMILES for the loaded molecule. If none available, use the empty string.
         * \return whether the parsing and initial conformer generation was successful.
         * \sa populateFromSMILES()
        */
        bool
        populateFromPDBFile(const std::string &filename, const std::string &smiles_hint = "", unsigned int seed = 36754,
                            std::vector<std::shared_ptr<InputPostProcessor::InputPostProcessor> > postProcessors = {});


        //! Return the number of atom in the Molecule, Hydrogen excepted
        unsigned int numberOfAtoms();

        //! Return the number of bonds in the Molecule
        unsigned int numberOfBonds();

        //! Obtain the conformer that was either loaded from file or generated when populating the Molecule
        /*!
         *
         * If the Molecule was loaded from a file format containing coordinate (PDB, ...), the initial conformer has those
         * coordinate. Otherwise, it is a RDKit generated conformer.
         *
         * \param centroidNormalization If true, gives the coordinates with the molecule centroid as 0,0,0. This is useful for easier application of rotation
         * \return Returns iConformer for initial conformer
         * \sa generateConformer(), generateConformers()
        */
        iConformer getInitialConformer(bool centroidNormalization = false) const;

        //! Generate a conformer of this molecule
        /*!
         * \param conformer Reference to conformer object to be filled. Not modified if generation fails.
         * \param centroidNormalization If true, gives the coordinates with the molecule centroid as 0,0,0. This is useful for easier application of rotation
         * \param seed RNG seed
         * \return Returns true if successfully generated conformer, false otherwise
         * \sa generateConformers()
        */
        bool generateConformer(iConformer &conformer, bool centroidNormalization = false, int seed = 367454);

        //! Generate a vector of conformer of this molecule
        /*!
         * \param viConformers Pointer to vector to be filled with conformer. May already contain something.
         * \param num Number of conformers to attempt to generate.
         * \param centroidNormalization If true, gives the coordinates with the molecule centroid as 0,0,0. This is useful for easier application of rotation
         * \param seed RNG seed
         * \return Returns the number of conformer actually generated. May be 0.
         * \sa generateConformer()
        */
        unsigned int
        generateConformers(std::vector<iConformer> &viConformers, unsigned int num, bool centroidNormalization = false,
                           int seed = 367454);

        //! Get the 3-letter residue code used for this molecule in PDB file
        /*!
         * \return A 3-letter residue name string
         * \sa setResidueName()
        */
        std::string getResidueName() const;

        //! Set the 3-letter residue code used for this molecule in PDB file
        /*!
         * \param res_name Residue name. Characters beyond the first three are ignored.
         * \return A 3-letter residue name string
         * \sa setResidueName()
        */
        void setResidueName(const std::string &res_name);

        bool updateAtomPositionsFromiConformer(const iConformer &conformer);

        //! Make a copy of the Molecule object where the underlying Atom and Bond are new objects instead of references to
        /*!
         * \return A new Molecule object with copied values
        */
        Molecule deepcopy();

        //! Number of rotatable bonds between heavy (= non-hydrogen) atoms
        /*!
         * \return Number of rotatable bonds
        */
        unsigned int getNumRotatableBond();


        bool operator==(const Molecule &rhs) const;

        bool operator!=(const Molecule &rhs) const;

    private:

        static unsigned int LastMolID;
        unsigned int molID;

        std::vector<std::shared_ptr<Atom> > atoms;
        std::vector<std::shared_ptr<Bond> > bonds;

        // unsigned int, unsigned int : idx of end atom in RWMol, vector<unsigned int> : atom idx rotated by such bonds
        std::vector<std::tuple<std::shared_ptr<Bond>, unsigned int, unsigned int, std::vector<unsigned int> > > rotatableBonds;

    public:
        std::shared_ptr<RDKit::RWMol> rwmol;
    private:
        std::string smiles;


        std::string residue_name = "LIG";

        //! Use this if using RDKit RWMol as entry.
        bool populateInternalAtomAndBondFromRWMol(unsigned int seed,
                                                  std::vector<std::shared_ptr<InputPostProcessor::InputPostProcessor> > postProcessors);

        int initial_conformer_id = -1;

        iConformer
        generateIConformerForGivenRDKitConformerID(unsigned int id, bool centroidNormalization = false) const;


        unsigned int numberOfRotatableBonds = 0;

        unsigned int updatedConformerID = 0;
    };

}

#endif //SMOLDOCK_MOLECULE_H
