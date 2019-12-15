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

#ifndef SMOLDOCK_PROTEIN_H
#define SMOLDOCK_PROTEIN_H


#include <vector>
#include <memory>


#include <ESBTL/default.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/weighted_atom_iterator.h>
#include <ESBTL/compressed_ifstream.h>


#include <Engines/Internals/iProtein.h>
#include <Structures/InputModifiers/InputModifierInterface.h>

#include "Structure.h"
#include "AminoAcid.h"

namespace SmolDock {

    /*!
     * \brief Class representing "rich" protein. Contains all extended attributes and functions.
     *
     */




    class Protein : Structure {

    public:
        Protein();

        //! Populate Protein from a PDB file.
        /*!
         *
         * Alternative location, temperature factor and other subtleties of the PDB format are ignored.
         *
         * \param filename Path to PDB file
         * \param postProcessors A vector of pointer to post-processors object that will be applied to the protein
         * \return whether the parsing was successful.
        */
        bool populateFromPDB(const std::string &filename,
                             std::vector<std::shared_ptr<InputModifier::InputModifier> > modifiers = {});

        bool populateFromPDBString(const std::string &PDB_Block,
                             std::vector<std::shared_ptr<InputModifier::InputModifier> > modifiers = {});

        //! Return an iProtein object for use in docking engine from the complete Protein
        /*!
         *
         * \return the equivalent iProtein object
        */
        iProtein getiProtein() const;

        /*! Returns an an iProtein object for use in docking engine, comprised of only the residue in the sphere of
         *  given center and radius + margin, considering the centroid of the residue.
         *
         * The radius and center value are filled in the resulting iProtein structure. The margin value is added to the
         * radius during atom selection, but not to the iProtein. This allow the selection of more residues than strictly
         * necessary to smooth out docking at the edge of the sphere. Recommended value 2 angstrom
         *
         * \param center center of the sphere (coordinate unit : angstrom)
         * \param radius radius of the sphere (unit : angstrom)
         * \param margin added to the radius during selection, but not to the radius specified in the iProtein (unit : angstrom)
         * \return
         */
        iProtein getPartialiProtein_sphere(std::array<double, 3> center, double radius, double margin) const;

        /**
         * Returns the radius of the smallest sphere bounding the entire protein. It is also the maximum distance between an atom
         * and the protein centroid.
         *
         */
        double getMaxRadius() const;

        /**
         * Tag the specified residue atom(s) with the special type given. Allow the use of specific scoring function like CovalentReversible
         * variant of common fields.
         *
         * When using this function for CovalentReversible scoring, you must tag the relevant active site residue that would form
         * covalent bond, for example serine in a catalytic triad bases protease.
         *
         *
         * \param resType Residue type
         * \param serialNumber Serial number of the residue in the chain
         * \param specialType The type to tag
         * \param ignoreMismatchingResType If the residue at the given serial number does not match the residue type, and this is true, warn but continue anayway
         * \return true on success
         *
         * \sa SpecialResidueTyping
         */
        bool applySpecialResidueTyping(const AminoAcid::AAType resType,
                                        const unsigned int serialNumber,
                                        const SpecialResidueTyping specialType,
                                        const bool ignoreMismatchingResType = false);



    protected:

        void populateFromESBTLSystems(std::string friendlyName, std::vector<ESBTL::Default_system>& systems,
                std::vector<std::shared_ptr<InputModifier::InputModifier> > modifiers);

        std::vector<std::shared_ptr<AminoAcid> > aminoacids;
        std::vector<std::shared_ptr<Atom> > heteroatoms;

        double center_x = 0;
        double center_y = 0;
        double center_z = 0;

        double min_x = 0;
        double min_y = 0;
        double min_z = 0;

        double max_x = 0;
        double max_y = 0;
        double max_z = 0;

        double max_distance_to_center = 0;




    };

}

#endif //SMOLDOCK_PROTEIN_H
