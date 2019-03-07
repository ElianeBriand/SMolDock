//
// Created by eliane on 04/03/19.
//

#ifndef SMOLDOCK_VINALIKECOMMON_H
#define SMOLDOCK_VINALIKECOMMON_H

#include <boost/log/trivial.hpp>

#include <cmath>
#include <cassert>

#include <Structures/Atom.h>


#define CONST_FUN_ATTR  ;;

namespace SmolDock::Score {


    // ////////////////// VINA CODE //////////////////////////////////////////////////
    // Until the closing "VINA CODE END" comment, the following code is adapted from Autodock Vina
    // Copyright (c) 2006-2010, The Scripps Research Institute
    // Released under the Apache 2.0 licence   http://www.apache.org/licenses/LICENSE-2.0
    // See COPYING for more details on licence information


    CONST_FUN_ATTR inline bool isHydrophobic(const unsigned char atomicNumber, const unsigned int variantFlags) noexcept {
        return (atomicNumber == 6 && (variantFlags & ((const unsigned int) Atom::AtomVariant::apolar))) || // C && apolar
               atomicNumber == 9 || // F
               atomicNumber == 17 || // Cl
               atomicNumber == 35 || // Br
               atomicNumber == 53; // I
    }


    CONST_FUN_ATTR inline bool isHydrogenAcceptor(const unsigned char atomicNumber,const unsigned int atomVariantFlags) noexcept {
        return (atomicNumber == 7 || //N
                atomicNumber == 8) // O
               && (atomVariantFlags & ((const unsigned int) Atom::AtomVariant::hydrogenAcceptor));
    }

    CONST_FUN_ATTR inline bool isHydrogenDonor(const unsigned char atomicNumber,const  unsigned int atomVariantFlags) noexcept {
        return (atomicNumber == 7 ||
                atomicNumber == 8)
               && (atomVariantFlags & ((const unsigned int) Atom::AtomVariant::hydrogenDonor));
    }

    CONST_FUN_ATTR inline bool hydrogenDonorAcceptorPair(const unsigned char atomicNumber1,const  unsigned int atom1VariantFlags,
                                                                  const unsigned char atomicNumber2, const unsigned int atom2VariantFlags) noexcept {
        return isHydrogenDonor(atomicNumber1, atom1VariantFlags) &&
               isHydrogenAcceptor(atomicNumber2, atom2VariantFlags);
    }

    CONST_FUN_ATTR inline bool hydrogenBondingPossible(const unsigned char atomicNumber1, const unsigned int atom1VariantFlags,
                                        unsigned char atomicNumber2, unsigned int atom2VariantFlags) noexcept {
        return hydrogenDonorAcceptorPair(atomicNumber1, atom1VariantFlags, atomicNumber2, atom2VariantFlags) ||
               hydrogenDonorAcceptorPair(atomicNumber2, atom2VariantFlags, atomicNumber1, atom1VariantFlags);
    }

    // ////////////////// VINA CODE END //////////////////////////////////////////////////

    CONST_FUN_ATTR inline bool isCovalentReversibleAcceptor(const unsigned char atomicNumber,const unsigned int atomVariantFlags) noexcept {
        return static_cast<bool>((atomVariantFlags & ((const unsigned int) Atom::AtomVariant::covalentReversibleAcceptor)));
    }

    CONST_FUN_ATTR inline bool isCovalentReversibleDonor(const unsigned char atomicNumber,const unsigned int atomVariantFlags) noexcept {
        return static_cast<bool>((atomVariantFlags & ((const unsigned int) Atom::AtomVariant::covalentReversibleDonor)));
    }

    CONST_FUN_ATTR inline bool covalentReversiblePair(const unsigned char atomicNumber1,const  unsigned int atom1VariantFlags,
                                                                  const unsigned char atomicNumber2, const unsigned int atom2VariantFlags) noexcept {
        return isCovalentReversibleDonor(atomicNumber1, atom1VariantFlags) &&
                isCovalentReversibleAcceptor(atomicNumber2, atom2VariantFlags);
    }


    CONST_FUN_ATTR inline bool covalentReversibleBondingPossible(const unsigned char atomicNumber1, const unsigned int atom1VariantFlags,
                                                                unsigned char atomicNumber2, unsigned int atom2VariantFlags) noexcept {
        return covalentReversiblePair(atomicNumber1, atom1VariantFlags, atomicNumber2, atom2VariantFlags) ||
                covalentReversiblePair(atomicNumber2, atom2VariantFlags, atomicNumber1, atom1VariantFlags);
    }


    // The following scoring function was implemented without knowledge of the Vina code proper, using papers such as
    // DOI:10.1371/journal.pone.0155183 or the original Vina paper, that describe the scoring function but do not
    // show code. The functions/constants above this comment (VINA CODE) were used after testing revealed small numerical
    // differences between this software and Vina.
    // The following function is thus thought to not be derived work under copyright law. However, if it actually is,
    // it could be distributed under the GPLv 3 as the original Vina code is released under Apache 2.0.


    namespace VinaClassic {
        constexpr const double coeff_gauss1 = -0.035579;
        constexpr const double coeff_gauss2 = -0.005156;
        constexpr const double coeff_repulsion = 0.840245;
        constexpr const double coeff_hydrophobic = -0.035069;
        constexpr const double coeff_hydrogen = -0.587439;
        constexpr const double coeff_entropic = 0.058459999999999998;

        constexpr const double interaction_cutoff = 8.0;
    }




    CONST_FUN_ATTR inline double distanceFromRawDistance(const double rawDistance, const double atomicRadiusLig, const double atomicRadiusProt) noexcept {
        return rawDistance - (atomicRadiusLig + atomicRadiusProt);
    }


    CONST_FUN_ATTR inline double vinaGaussComponent(const double distance, const double offset, const double multiplier) noexcept {
        return std::exp(-1 * std::pow((distance - offset)/ multiplier, 2));
    }

    CONST_FUN_ATTR inline double vinaRepulsionComponent(const double distance, const double cutoff) noexcept {
        if (distance < cutoff) {
            return std::pow(distance, 2);
        }else{
            return 0.0;
        }
    }

    CONST_FUN_ATTR inline double vinaHydrophobicComponent(const double distance, const unsigned char atomicNumberAtom1, const unsigned int variantFlagsAtom1,
                                                                   const unsigned char atomicNumberAtom2, const unsigned int variantFlagsAtom2) noexcept
    {
        if (isHydrophobic(atomicNumberAtom1, variantFlagsAtom1) &&
            isHydrophobic(atomicNumberAtom2, variantFlagsAtom2)) // "Hydrophobic" atoms
        {
            if (distance >= 1.5)
                return 0.0;
            if (distance <= 0.5)
                return 1.0;
            //then  (0.5 < distance && distance < 1.5) holds
            return (1.5 - distance);
        }
        return 0.0;
    }

    CONST_FUN_ATTR inline double vinaHydrogenComponent(const double distance, const unsigned char atomicNumberAtom1, const unsigned int variantFlagsAtom1,
                                                                const unsigned char atomicNumberAtom2, const unsigned int variantFlagsAtom2) noexcept
    {
        if (hydrogenBondingPossible(atomicNumberAtom1, variantFlagsAtom1, atomicNumberAtom2,
                                    variantFlagsAtom2)) // Hydrogen donor and acceptor
        {
            if (distance < -0.7) {
                return 1.0;
            } else if (distance < 0) //  // ==> distance between -0.7 and 0
            {
                return -distance / 0.7;
            }
        }
        return 0.0;
    }

    namespace VinaExtended {

        constexpr const double coeff_CovalentReversible = -2.5;


        CONST_FUN_ATTR inline double covalentReversibleComponent(const double distance, const unsigned char atomicNumberAtom1, const unsigned int variantFlagsAtom1,
                                                                            const unsigned char atomicNumberAtom2, const unsigned int variantFlagsAtom2) noexcept
        {
            if (covalentReversibleBondingPossible(atomicNumberAtom1, variantFlagsAtom1, atomicNumberAtom2,
                                        variantFlagsAtom2)) // Hydrogen donor and acceptor
            {
                // this corresponds to ~ 1.43 Angstrom between the two, aka classical C-O bond length
                // these value seem weird because VdW radius as understood by vina, and covalent bond length
                // are quite different (distance = 0 means there is 1.9+1.7 angstrom between the atoms, as given by the "atomic radii" table)
                // Which is too long for our purpose of modeling covalent bonds.
                // TODO: take into account the atom type (this is currently only for C-O covalent reversible)
                if (distance < -2.1 && distance > -2.2)
                {
                    return 1.0;
                }
            }
            return 0.0;
        }
    }




    CONST_FUN_ATTR inline double
    scoreForAtomCouple(const double distance,const unsigned char atom1AtomicNumber,const unsigned int atom1AtomVariant,
                       const unsigned char atom2AtomicNumber,const unsigned int atom2AtomVariant) noexcept {

        double score_intermol = 0.0;

        // exp −(d/0.5Å)^2
        const double gauss1 = vinaGaussComponent(distance, 0.0, 0.5);
        score_intermol += VinaClassic::coeff_gauss1 * gauss1;

        // exp −((d−3Å)/2Å)^2
        const double gauss2 = vinaGaussComponent(distance, 3.0, 2.0);
        score_intermol += VinaClassic::coeff_gauss2 * gauss2;

        const double repuls = vinaRepulsionComponent(distance, 0.0);
        score_intermol += VinaClassic::coeff_repulsion * repuls;

        const double hydrophobic = vinaHydrophobicComponent(distance,
                                                      atom1AtomicNumber, atom1AtomVariant,
                                                      atom2AtomicNumber, atom2AtomVariant);
        score_intermol += VinaClassic::coeff_hydrophobic * hydrophobic;


        const double hydrogen = vinaHydrogenComponent(distance,
                                                      atom1AtomicNumber, atom1AtomVariant,
                                                      atom2AtomicNumber, atom2AtomVariant);
        score_intermol += VinaClassic::coeff_hydrogen  * hydrogen;

        return score_intermol;
    }




}


#endif //SMOLDOCK_VINALIKECOMMON_H
