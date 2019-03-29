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

#ifndef SMOLDOCK_ICONFORMER_H
#define SMOLDOCK_ICONFORMER_H

#include <vector>
#include <memory>

#include <Vc/Vc>

#include "iTransform.h"

namespace SmolDock {

    //! Internal, fast representation of a given Molecule conformer for efficiency
    /*!
     * \sa iAtom, iProtein
    */


    struct iConformer {

        unsigned int num_rotatable_bond;

        Eigen::Translation<double, 3> centroidNormalizingTransform; /*!< If the coordinate have been normalized so as to have the centroid
                                                        as [0,0,0], then this contain the translation that restore it to its original position,
                                                        aka the coordinate of the centroid (because |normalized> = |originalPos> - |centroidPos> )
                                                        Else it's just 0,0,0 */


        std::vector<double> x, y, z;
        std::vector<double> atomicRadius;
        std::vector<unsigned char> type;
        std::vector<unsigned int> variant;

        std::vector<unsigned int> bondEnds1Index;
        std::vector<unsigned int> bondEnds2Index;
        std::vector<std::vector<unsigned int> > rotatableGroups;

    };

    struct iConformer_Vectorized {

        iConformer_Vectorized() = delete;

        iConformer_Vectorized(const iConformer& iconformer) :
                num_rotatable_bond(iconformer.num_rotatable_bond),
                centroidNormalizingTransform(iconformer.centroidNormalizingTransform),
                x(iconformer.x.size()),
                y(iconformer.y.size()),
                z(iconformer.z.size()),

                atomicRadius(iconformer.atomicRadius.size()),
                type(iconformer.type.size()),
                variant(iconformer.variant.size()),

                bondEnds1Index(iconformer.bondEnds1Index),
                bondEnds2Index(iconformer.bondEnds2Index),
                rotatableGroups(iconformer.rotatableGroups) {
            for (unsigned int i = 0; i < x.entriesCount(); ++i) {
                x.scalar(i) = iconformer.x[i];
                y.scalar(i) = iconformer.y[i];
                z.scalar(i) = iconformer.z[i];

                atomicRadius.scalar(i) = iconformer.atomicRadius[i];
                type.scalar(i) = iconformer.type[i];
                variant.scalar(i) = iconformer.variant[i];
            }

            for (unsigned int i = 0; i < x.vectorsCount(); ++i) {
                for (unsigned int j = 0; j < Vc::Vector<double>::Size ; ++j) {
                    std::cout << " -> " << i << ", " << j << std::endl;
                    std::cout << "    " << x.vector(i)[j] << " <-> " << iconformer.x[Vc::Vector<double>::Size * i + j] << std::endl;
                    std::cout << "    " << y.vector(i)[j] << " <-> " << iconformer.y[Vc::Vector<double>::Size * i + j] << std::endl;
                    std::cout << "    " << z.vector(i)[j] << " <-> " << iconformer.z[Vc::Vector<double>::Size * i + j] << std::endl;
                }
            }
        }

        unsigned int num_rotatable_bond;

        Eigen::Translation<double, 3> centroidNormalizingTransform;

        Vc::Memory<Vc::Vector<double>> x, y, z;

        Vc::Memory<Vc::Vector<double>> atomicRadius;
        Vc::Memory<Vc::Vector<unsigned char>> type;
        Vc::Memory<Vc::Vector<unsigned int>> variant;

        std::vector<unsigned int> bondEnds1Index;
        std::vector<unsigned int> bondEnds2Index;
        std::vector<std::vector<unsigned int> > rotatableGroups;
    };

}


#endif //SMOLDOCK_ICONFORMER_H
