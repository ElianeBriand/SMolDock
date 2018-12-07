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

#ifndef SMOLDOCK_PROTEIN_H
#define SMOLDOCK_PROTEIN_H


#include <vector>
#include <memory>



#include <ESBTL/default.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/weighted_atom_iterator.h>
#include <ESBTL/compressed_ifstream.h>


#include <Engines/Internals/iProtein.h>

#include "Structure.h"
#include "AminoAcid.h"

namespace SmolDock {


    class Protein : Structure {

    public:
        Protein() = default;

        bool populateFromPDB(const std::string &filename);

        iProtein getiProtein();

    private:
        std::vector<std::shared_ptr<AminoAcid> > aminoacids;
        std::vector<std::shared_ptr<Atom> > heteroatoms;

        double center_x = 0;
        double center_y = 0;
        double center_z = 0;

    };

}

#endif //SMOLDOCK_PROTEIN_H
