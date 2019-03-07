//
// Created by briand on 3/7/19.
//

#ifndef SMOLDOCK_VINAEXTENDEDCOMMON_H
#define SMOLDOCK_VINAEXTENDEDCOMMON_H

#include "VinaLikeCommon.h"


namespace SmolDock::Score::VinaExtended {

    constexpr const double coeff_CovalentReversible = -25.0;


    __attribute__((const)) inline double covalentReversibleComponent(const double distance, const unsigned char atomicNumberAtom1, const unsigned int variantFlagsAtom1,
                                                                     const unsigned char atomicNumberAtom2, const unsigned int variantFlagsAtom2) noexcept
    {
        if (covalentReversibleBondingPossible(atomicNumberAtom1, variantFlagsAtom1, atomicNumberAtom2,
                                              variantFlagsAtom2)) // Hydrogen donor and acceptor
        {
            // Given an ideal bond length of 1.43-1.63 Angstrom
            // and linear interpolation between 2 angstrom and 1.53, and 1.43 and -inf, respectively
            // We would have :
            //            if (distance < 2.0)
            //            {
            //                if(distance > 1.43 && distance < 1.53)
            //                {
            //                    return 1.0;
            //                }
            //                if(distance < 1.43)
            //                {
            //                    return -0.8 * distance;
            //                }
            //                return 0.47 * distance;
            //            }
            // But the distance argument is calculated by substrating VdW radii.
            // Given a C-O VdW radii of 1.9+1.7 = 3.6, we need to add that back to the distance.
            const double realDistance = distance + 3.6;
            // As that would give negative ideal distance, we need to use the RepulsionExceptForCovalentComponent instead of simple repulsion
            // TODO: take into account the atom type (this is currently only for C-O covalent reversible)
            // Thiw would mean changing the distance to add.

            if (realDistance < 2.0)
            {
                if(realDistance > 1.43 && realDistance < 1.53)
                {
                    return 1.0;
                }
                if(realDistance < 1.43)
                {
                    return -0.8 * realDistance;
                }
                return 0.47 * realDistance;
            }
        }
        return 0.0;
    }

    __attribute__((const)) inline double RepulsionExceptForCovalentComponent(const double distance, const double cutoff,
                                                                             const unsigned char atomicNumberAtom1, const unsigned int variantFlagsAtom1,
                                                                             const unsigned char atomicNumberAtom2, const unsigned int variantFlagsAtom2) noexcept {
        if (! covalentReversibleBondingPossible(atomicNumberAtom1, variantFlagsAtom1, atomicNumberAtom2,
                                              variantFlagsAtom2))
        {
            if (distance < cutoff) {
                return std::pow(distance, 2);
            }else{
                return 0.0;
            }
        }
        return 0.0;

    }
}

#endif //SMOLDOCK_VINAEXTENDEDCOMMON_H
