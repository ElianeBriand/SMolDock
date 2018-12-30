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

#include "VinaLikeScoringFunction.h"

#undef BOOST_LOG
#include <boost/log/trivial.hpp>

#include <cmath>
#include <cassert>
#include <iomanip>
#include <Structures/Atom.h>

namespace SmolDock {
    namespace Score {


        // ////////////////// VINA CODE //////////////////////////////////////////////////
        // Until the closing "VINA CODE END" comment, the following code is adapted from Autodock Vina
        // Copyright (c) 2006-2010, The Scripps Research Institute
        // Released under the Apache 2.0 licence   http://www.apache.org/licenses/LICENSE-2.0
        // See COPYING for more details on licence information

        inline bool isHydrophobic(unsigned char atomicNumber, unsigned int variantFlags) {
            return (atomicNumber == 6 && (variantFlags & ((unsigned int) Atom::AtomVariant::apolar) ))|| // C && apolar
                    atomicNumber == 9 || // F
                    atomicNumber == 17 || // Cl
                    atomicNumber == 35 || // Br
                    atomicNumber == 53; // I
        }


        inline bool isHydrogenAcceptor(unsigned char atomicNumber) {
            return atomicNumber == 7 || //N
                    atomicNumber == 8; // O
        }

        inline bool isHydrogenDonor(unsigned char atomicNumber) {
            return atomicNumber == 7 ||
                    atomicNumber == 8;
            // TODO : better differentiate donor atom and non-donor using variant
        }

        inline bool hydrogenDonorAcceptorPair(unsigned char atomicNumber1, unsigned char atomicNumber2) {
            return isHydrogenDonor(atomicNumber1) && isHydrogenAcceptor(atomicNumber2);
        }

        inline bool hydrogenBondingPossible(unsigned char atomicNumber1, unsigned char atomicNumber2) {
            return hydrogenDonorAcceptorPair(atomicNumber1, atomicNumber2) || hydrogenDonorAcceptorPair(atomicNumber2, atomicNumber1);
        }

        // ////////////////// VINA CODE END //////////////////////////////////////////////////


        // The following scoring function was implemented without knowledge of the Vina code proper, using papers such as
        // DOI:10.1371/journal.pone.0155183 or the original Vina paper, that describe the scoring function but do not
        // show code. The functions/constants above this comment (VINA CODE) were used after testing revealed small numerical
        // differences between this software and Vina.
        // The following function is thus thought to not be derived work under copyright law. However, if it actually is,
        // it could be distributed under the GPLv 3 as the original Vina code is released under Apache 2.0.

        double score_gauss1 = 0;
        double score_gauss2 = 0;
        double score_repulsion = 0;
        double score_hydrophobic = 0;
        double score_hydrogen = 0;

        inline double scoreForAtomCouple(double distance, unsigned char atom1AtomicNumber, unsigned int atom1AtomVariant,
                                                            unsigned char atom2AtomicNumber, unsigned int atom2AtomVariant)
        {

            double score = 0.0;
            /*
                    // Special rubber band term not in Vina scoring func
                    double x_diff_center = std::pow(LigPosition.x - protein.center_x, 2);
                    double y_diff_center = std::pow(LigPosition.y - protein.center_y, 2);
                    double z_diff_center = std::pow(LigPosition.z - protein.center_z, 2);
                    double rawDistCenter = std::sqrt(x_diff_center + y_diff_center + z_diff_center);

                    score += 0.00001 * rawDistCenter;
            //*/

            // exp −(d/0.5Å)^2
            double gauss1 =  std::exp(-1 * std::pow(distance / 0.5, 2));
            score += -0.035579 * gauss1;
            score_gauss1 += gauss1;

            // exp −((d−3Å)/2Å)^2
            double gauss2 =  std::exp(-1 * std::pow((distance - 3.0) / 2.0, 2));
            score += -0.005156 * gauss2;
            score_gauss2 += gauss2;

            if (distance < 0)
            {
                double repuls =  std::pow(distance, 2);
                score += 0.840245 * repuls;
                score_repulsion += repuls;
            }

            if ( isHydrophobic(atom1AtomicNumber, atom1AtomVariant) && isHydrophobic(atom2AtomicNumber, atom2AtomVariant)) // "Hydrophobic" atoms
            {
                double hydrophobic_contrib = 0;
                if(distance >= 1.5)
                    hydrophobic_contrib = 0;
                if(distance <= 0.5)
                    hydrophobic_contrib = 1;
                if(0.5 < distance &&  distance < 1.5 )
                    hydrophobic_contrib = (1.5 - distance);

                score += -0.035069 * hydrophobic_contrib;
                score_hydrophobic += hydrophobic_contrib;
            }

            if (hydrogenBondingPossible(atom1AtomicNumber,atom2AtomicNumber)) // Hydrogen donor and acceptor
            {
                double hbond_contrib = 0;
                if(distance < - 0.7)
                {
                    hbond_contrib = 1;
                }else if(distance < 0) //  // ==> distance between -0.7 and 0
                {
                    hbond_contrib = -distance/0.7;
                }
                score += -0.587439*hbond_contrib;
                score_hydrogen += hbond_contrib;

            }
            return score;
        }


        double vina_like_rigid_inter_scoring_func(const iConformer &ligand, const iTransform& transform,const iProtein &protein) {



            double score = 0;



            assert(std::abs(quaternionNorm(transform.rota) - 1) < 0.01);

            for (unsigned int idxLig = 0; idxLig < ligand.x.size(); idxLig++) {
                for (unsigned int idxProt = 0; idxProt < protein.x.size(); idxProt++) {


                    iVect LigPosition = {ligand.x[idxLig],ligand.y[idxLig],ligand.z[idxLig]};
                    applyTransformInPlace(LigPosition, transform);

                    double x_diff_sq = std::pow(LigPosition.x - protein.x[idxProt], 2);
                    double y_diff_sq = std::pow(LigPosition.y - protein.y[idxProt], 2);
                    double z_diff_sq = std::pow(LigPosition.z - protein.z[idxProt], 2);


                    double atomicRadiusLig = ligand.atomicRadius[idxLig];
                    double atomicRadiusProt = protein.atomicRadius[idxProt];
                    double rawDist = std::sqrt(x_diff_sq + y_diff_sq + z_diff_sq);

                    const double cutoff = 8.0;
                    if(rawDist >= cutoff)
                        continue;

                    double radToRemove = (atomicRadiusLig + atomicRadiusProt);
                    double distance = rawDist - radToRemove;

                    score += scoreForAtomCouple(distance,ligand.type[idxLig],ligand.variant[idxLig],protein.type[idxProt],protein.variant[idxProt]);

                } // for
            } // for
            /*
            BOOST_LOG_TRIVIAL(debug) << "Intermolecular scoring contribution without weighting :";
            BOOST_LOG_TRIVIAL(debug) << "Gauss 1      : " << std::fixed<< std::setprecision(5) << score_gauss1;
            BOOST_LOG_TRIVIAL(debug) << "Gauss 2      : " << std::fixed<< std::setprecision(5) << score_gauss2;
            BOOST_LOG_TRIVIAL(debug) << "Repulsion    : " << std::fixed<< std::setprecision(5) << score_repulsion;
            BOOST_LOG_TRIVIAL(debug) << "Hydrophobic  : " << std::fixed<< std::setprecision(5) << score_hydrophobic;
            BOOST_LOG_TRIVIAL(debug) << "Hydrogen     : " << std::fixed<< std::setprecision(5) << score_hydrogen;
            BOOST_LOG_TRIVIAL(debug) << "------------------------------------------";
            BOOST_LOG_TRIVIAL(debug) << "Raw Score    : " << score;
            BOOST_LOG_TRIVIAL(debug) << "Nrotatable   : " << ligand.num_rotatable_bond;
            */
            double final_score = score / (1 + (0.058459999999999998*ligand.num_rotatable_bond));
            /*
            BOOST_LOG_TRIVIAL(debug) << "------------------------------------------";
            BOOST_LOG_TRIVIAL(debug) << "Final Score  : " << final_score;
            */
            return final_score;
        }


    }
}