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
#include <Structures/InputPostProcessors/InputPostProcessorInterface.h>

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
                             std::vector<std::shared_ptr<InputPostProcessor::InputPostProcessor> > postProcessors = {});

        //! Return an iProtein object for use in docking engine
        /*!
         *
         * \return the equivalent iProtein object
        */
        iProtein getiProtein() const;

        double getMaxRadius() const;

    private:
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
