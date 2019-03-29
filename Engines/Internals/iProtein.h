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

#include <Vc/Vc>



namespace SmolDock {

    //! Internal, fast representation of Protein for efficiency
    /*!
     * \sa iAtom, iConformer
    */
    struct iProtein {

        iProtein() {

        }


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


    struct iProtein_vectorized {

        iProtein_vectorized() = delete;

        iProtein_vectorized(const iProtein& p) :
                center_x(p.center_x), center_y(p.center_y), center_z(p.center_z),
                radius(p.radius),
                x(p.x.size()),
                y(p.y.size()),
                z(p.z.size()),
                atomicRadius(p.atomicRadius.size()),
                type(p.type.size()),
                variant(p.variant.size()),
                AAId_to_AtomPositionInVect(p.AAId_to_AtomPositionInVect) {
            for (unsigned int i = 0; i < x.entriesCount(); ++i) {
                x.scalar(i) = p.x[i];
                y.scalar(i) = p.y[i];
                z.scalar(i) = p.z[i];

                atomicRadius.scalar(i) = p.atomicRadius[i];
                type.scalar(i) = p.type[i];
                variant.scalar(i) = p.variant[i];
            }

        }


        //! Pseudo-center of protein as a mean of each atom coordinate
        double center_x, center_y, center_z;
        double radius;

        Vc::Memory <Vc::Vector<double>> x, y, z;
        Vc::Memory <Vc::Vector<double>> atomicRadius;
        Vc::Memory <Vc::Vector<unsigned char>> type;
        Vc::Memory <Vc::Vector<unsigned int>> variant;

        //! Map the AAId of residues to the corresponding position of the atoms in the vector
        std::map<unsigned int, std::tuple<unsigned long, unsigned long> > AAId_to_AtomPositionInVect;
    };

}


#endif //SMOLDOCK_IPROTEIN_H
