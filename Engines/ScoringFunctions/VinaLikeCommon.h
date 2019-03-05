//
// Created by eliane on 04/03/19.
//

#ifndef SMOLDOCK_VINALIKECOMMON_H
#define SMOLDOCK_VINALIKECOMMON_H

#include <boost/log/trivial.hpp>

#include <cmath>
#include <cassert>

#include <Structures/Atom.h>

namespace SmolDock::Score {


    // ////////////////// VINA CODE //////////////////////////////////////////////////
    // Until the closing "VINA CODE END" comment, the following code is adapted from Autodock Vina
    // Copyright (c) 2006-2010, The Scripps Research Institute
    // Released under the Apache 2.0 licence   http://www.apache.org/licenses/LICENSE-2.0
    // See COPYING for more details on licence information

    inline bool isHydrophobic(unsigned char atomicNumber, unsigned int variantFlags) {
        return (atomicNumber == 6 && (variantFlags & ((unsigned int) Atom::AtomVariant::apolar))) || // C && apolar
               atomicNumber == 9 || // F
               atomicNumber == 17 || // Cl
               atomicNumber == 35 || // Br
               atomicNumber == 53; // I
    }


    inline bool isHydrogenAcceptor(unsigned char atomicNumber, unsigned int atomVariantFlags) {
        return (atomicNumber == 7 || //N
                atomicNumber == 8) // O
               && (atomVariantFlags & ((unsigned int) Atom::AtomVariant::hydrogenAcceptor));
    }

    inline bool isHydrogenDonor(unsigned char atomicNumber, unsigned int atomVariantFlags) {
        return (atomicNumber == 7 ||
                atomicNumber == 8)
               && (atomVariantFlags & ((unsigned int) Atom::AtomVariant::hydrogenDonor));
    }

    inline bool hydrogenDonorAcceptorPair(unsigned char atomicNumber1, unsigned int atom1VariantFlags,
                                          unsigned char atomicNumber2, unsigned int atom2VariantFlags) {
        return isHydrogenDonor(atomicNumber1, atom1VariantFlags) &&
               isHydrogenAcceptor(atomicNumber2, atom2VariantFlags);
    }

    inline bool hydrogenBondingPossible(unsigned char atomicNumber1, unsigned int atom1VariantFlags,
                                        unsigned char atomicNumber2, unsigned int atom2VariantFlags) {
        return hydrogenDonorAcceptorPair(atomicNumber1, atom1VariantFlags, atomicNumber2, atom2VariantFlags) ||
               hydrogenDonorAcceptorPair(atomicNumber2, atom2VariantFlags, atomicNumber1, atom1VariantFlags);
    }

    // ////////////////// VINA CODE END //////////////////////////////////////////////////


    // The following scoring function was implemented without knowledge of the Vina code proper, using papers such as
    // DOI:10.1371/journal.pone.0155183 or the original Vina paper, that describe the scoring function but do not
    // show code. The functions/constants above this comment (VINA CODE) were used after testing revealed small numerical
    // differences between this software and Vina.
    // The following function is thus thought to not be derived work under copyright law. However, if it actually is,
    // it could be distributed under the GPLv 3 as the original Vina code is released under Apache 2.0.



    inline double
    scoreForAtomCouple(double distance, unsigned char atom1AtomicNumber, unsigned int atom1AtomVariant,
                       unsigned char atom2AtomicNumber, unsigned int atom2AtomVariant) {

        double score_intermol = 0.0;

        double score_gauss1 = 0;
        double score_gauss2 = 0;
        double score_repulsion = 0;
        double score_hydrophobic = 0;
        double score_hydrogen = 0;

        // exp −(d/0.5Å)^2
        double gauss1 = std::exp(-1 * std::pow(distance / 0.5, 2));
        score_intermol += -0.035579 * gauss1;
        score_gauss1 += gauss1;

        // exp −((d−3Å)/2Å)^2
        double gauss2 = std::exp(-1 * std::pow((distance - 3.0) / 2.0, 2));
        score_intermol += -0.005156 * gauss2;
        score_gauss2 += gauss2;

        if (distance < 0) {
            double repuls = std::pow(distance, 2);
            score_intermol += 0.840245 * repuls;
            score_repulsion += repuls;
        }

        if (isHydrophobic(atom1AtomicNumber, atom1AtomVariant) &&
            isHydrophobic(atom2AtomicNumber, atom2AtomVariant)) // "Hydrophobic" atoms
        {
            double hydrophobic_contrib = 0;
            if (distance >= 1.5)
                hydrophobic_contrib = 0;
            if (distance <= 0.5)
                hydrophobic_contrib = 1;
            if (0.5 < distance && distance < 1.5)
                hydrophobic_contrib = (1.5 - distance);

            score_intermol += -0.035069 * hydrophobic_contrib;
            score_hydrophobic += hydrophobic_contrib;
        }

        if (hydrogenBondingPossible(atom1AtomicNumber, atom1AtomVariant, atom2AtomicNumber,
                                    atom2AtomVariant)) // Hydrogen donor and acceptor
        {
            double hbond_contrib = 0;
            if (distance < -0.7) {
                hbond_contrib = 1;
            } else if (distance < 0) //  // ==> distance between -0.7 and 0
            {
                hbond_contrib = -distance / 0.7;
            }
            score_intermol += -0.587439 * hbond_contrib;
            score_hydrogen += hbond_contrib;

        }
        return score_intermol;
    }


}


#endif //SMOLDOCK_VINALIKECOMMON_H
