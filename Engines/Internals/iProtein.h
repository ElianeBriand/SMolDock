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

#ifndef SMOLDOCK_IPROTEIN_H
#define SMOLDOCK_IPROTEIN_H

#include <vector>
#include <memory>
#include <map>


namespace SmolDock {

    //! Internal, fast representation of Protein for efficiency
    /*!
     * \sa iAtom, iConformer
    */
    struct iProtein {
        //! Pseudo-center of protein as a mean of each atom coordinate
        double center_x, center_y, center_z;
        double radius;

        std::vector<double> x, y, z;
        std::vector<double> atomicRadius;
        std::vector<unsigned char> type;
        std::vector<unsigned int> variant;

        //! Map the AAId of residues to the corresponding position of the atoms in the vector
        std::map<unsigned int, std::tuple<unsigned long, unsigned long> > AAId_to_AtomPositionInVect;
    };

}


#endif //SMOLDOCK_IPROTEIN_H
