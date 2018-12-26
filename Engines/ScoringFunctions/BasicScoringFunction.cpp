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

#include "BasicScoringFunction.h"

#undef BOOST_LOG
#include <boost/log/trivial.hpp>

#include <cmath>
#include <cassert>


namespace SmolDock {
    namespace Score {


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
                    double radToRemove = (atomicRadiusLig + atomicRadiusProt);
                    double distance = rawDist - radToRemove;

/*
                    // Special rubber band term not in Vina scoring func
                    double x_diff_center = std::pow(LigPosition.x - protein.center_x, 2);
                    double y_diff_center = std::pow(LigPosition.y - protein.center_y, 2);
                    double z_diff_center = std::pow(LigPosition.z - protein.center_z, 2);
                    double rawDistCenter = std::sqrt(x_diff_center + y_diff_center + z_diff_center);

                    score += 0.00001 * rawDistCenter;
//*/

                    // exp −(d/0.5Å)^2
                    double gauss1 = -0.0356 * std::exp(-1 * std::pow(distance / 0.5, 2));
                    score += gauss1;

                    // exp −((d−3Å)/2Å)^2
                    double gauss2 = -0.00516 * std::exp(-1 * std::pow((distance - 3) / 2, 2));
                    score += gauss2;


                    if (distance < 0)
                        score += 0.840 * std::pow(distance, 2);
//*
                    unsigned char ligandAtomicNum = ligand.type[idxLig];
                    unsigned char proteinAtomicNum = protein.type[idxProt];


                    if (ligandAtomicNum == 6 && proteinAtomicNum == 6) // "Hydrophobic" atoms
                    {
                        double hydrophobic_contrib = 0;
                        if (distance < 0.5) {
                            hydrophobic_contrib = 1;
                        } else if (distance < 1.5) // ==> distance between 0.5 and 1.5
                        {
                            hydrophobic_contrib = 1.5 - distance;
                        }
                        score += -0.0351 * hydrophobic_contrib;
                    }

                    if ((ligandAtomicNum == 8 && proteinAtomicNum == 8) ||
                            (ligandAtomicNum == 7 && proteinAtomicNum == 8) ||
                            (ligandAtomicNum == 8 && proteinAtomicNum == 7) ||
                            (ligandAtomicNum == 7 && proteinAtomicNum == 8)) // Hydrogen donor and acceptor
                        // TODO : better conditional for hydrogen bond donor and acceptor
                    {
                        double hbond_contrib = 0;
                        if(distance < - 0.7)
                        {
                            hbond_contrib = 1;
                        }else if(distance < 0) //  // ==> distance between -0.7 and 0
                        {
                            hbond_contrib = -distance/0.7;
                        }
                        score += -0.597*hbond_contrib;

                    }
//*/
                } // for
            } // for

            return score;
        }


    }
}