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

#include "VinaLikeRigidScoringFunction.h"

#undef BOOST_LOG

#include <boost/log/trivial.hpp>

#include <cmath>
#include <cassert>
#include <iomanip>
#include <Structures/Atom.h>

#include <Engines/Internals/InternalsUtilityFunctions.h>

#include "VinaLikeCommon.h"

namespace SmolDock {
    namespace Score {


        double vina_like_rigid_inter_scoring_func(const iConformer &ligand, const iTransform &transform,
                                                  const iProtein &protein) {

            assert(!ligand.x.empty());
            assert(!protein.x.empty());
            assert(std::abs(transform.rota.norm() - 1) < 0.01);


            double score_raw = 0;


            Eigen::Vector3d ProtCenterPosition = {protein.center_x, protein.center_y, protein.center_z};

            for (unsigned int idxLig = 0; idxLig < ligand.x.size(); idxLig++) {
                for (unsigned int idxProt = 0; idxProt < protein.x.size(); idxProt++) {


                    Eigen::Vector3d LigPosition = {ligand.x[idxLig], ligand.y[idxLig], ligand.z[idxLig]};
                    applyRigidTransformInPlace(LigPosition, transform);

                    Eigen::Vector3d ProtPosition = {protein.x[idxProt], protein.y[idxProt], protein.z[idxProt]};

                    Eigen::Vector3d distToCenterVector = LigPosition - ProtCenterPosition;

                    double distanceToProteinCenter = distToCenterVector.norm();

                    if (distanceToProteinCenter > (protein.radius - 1)) {
                        score_raw += std::pow((distanceToProteinCenter - protein.radius), 4) + 10;
                        continue;
                    }

                    Eigen::Vector3d distVect = ProtPosition - LigPosition;

                    double atomicRadiusLig = ligand.atomicRadius[idxLig];
                    double atomicRadiusProt = protein.atomicRadius[idxProt];
                    double rawDist = distVect.norm();


                    const double cutoff = 8.0;
                    if (rawDist >= cutoff)
                        continue;

                    double radToRemove = (atomicRadiusLig + atomicRadiusProt);
                    double distance = rawDist - radToRemove;

                    score_raw += scoreForAtomCouple(distance, ligand.type[idxLig], ligand.variant[idxLig],
                                                    protein.type[idxProt], protein.variant[idxProt]);
/*
                    // Special rubber band term not in Vina scoring func
                    double x_diff_center = std::pow(LigPosition.x - protein.center_x, 2);
                    double y_diff_center = std::pow(LigPosition.y - protein.center_y, 2);
                    double z_diff_center = std::pow(LigPosition.z - protein.center_z, 2);
                    double rawDistCenter = std::sqrt(x_diff_center + y_diff_center + z_diff_center);

                    score += 0.00001 * rawDistCenter;
*/

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
            BOOST_LOG_TRIVIAL(debug) << "Raw Score    : " << score_raw;
            BOOST_LOG_TRIVIAL(debug) << "Nrotatable   : " << ligand.num_rotatable_bond;
            //*/
            double final_score = score_raw / (1 + (0.058459999999999998 * ligand.num_rotatable_bond));
            /*
            BOOST_LOG_TRIVIAL(debug) << "------------------------------------------";
            BOOST_LOG_TRIVIAL(debug) << "Final Score  : " << final_score;
            //*/

            return final_score;
        }

        double VinaLikeRigidScoringFunction::Evaluate(const arma::mat &x) {
            assert(x.n_rows == 7);

            iTransform tr = this->internalToExternalRepr(x);
            normalizeQuaternionInPlace(tr.rota);

            double score_ = vina_like_rigid_inter_scoring_func(this->startingConformation, tr, this->prot);

            return score_;
        }

        double VinaLikeRigidScoringFunction::EvaluateWithGradient(const arma::mat &x, arma::mat &grad) {

            assert(!x.has_nan());
            assert(!grad.has_nan());
            assert(x.n_rows == 7);
            assert(grad.n_rows == 7);

            iTransform tr = this->internalToExternalRepr(x);
            normalizeQuaternionInPlace(tr.rota);

            double score_ = vina_like_rigid_inter_scoring_func(this->startingConformation, tr, this->prot);

            // Translation
            {
                iTransform transform_dx = tr;
                transform_dx.transl.x() += this->differential_epsilon;
                grad[0] = vina_like_rigid_inter_scoring_func(this->startingConformation, transform_dx, this->prot) -
                          score_;
            }

            {
                iTransform transform_dy = tr;
                transform_dy.transl.y() += this->differential_epsilon;
                grad[1] = vina_like_rigid_inter_scoring_func(this->startingConformation, transform_dy, this->prot) -
                          score_;
            }

            {
                iTransform transform_dz = tr;
                transform_dz.transl.z() += this->differential_epsilon;
                grad[2] = vina_like_rigid_inter_scoring_func(this->startingConformation, transform_dz, this->prot) -
                          score_;
            }

            // Rotation

            {
                iTransform transform_ds = tr;
                transform_ds.rota.w() += this->differential_epsilon;
                normalizeQuaternionInPlace(transform_ds.rota);
                grad[3] = vina_like_rigid_inter_scoring_func(this->startingConformation, transform_ds, this->prot) -
                          score_;
            }

            {
                iTransform transform_du = tr;
                transform_du.rota.x() += this->differential_epsilon;
                normalizeQuaternionInPlace(transform_du.rota);
                grad[4] = vina_like_rigid_inter_scoring_func(this->startingConformation, transform_du, this->prot) -
                          score_;
            }

            {
                iTransform transform_dv = tr;
                transform_dv.rota.y() += this->differential_epsilon;
                normalizeQuaternionInPlace(transform_dv.rota);
                grad[5] = vina_like_rigid_inter_scoring_func(this->startingConformation, transform_dv, this->prot) -
                          score_;
            }

            {
                iTransform transform_dt = tr;
                transform_dt.rota.z() += this->differential_epsilon;
                normalizeQuaternionInPlace(transform_dt.rota);
                grad[6] = vina_like_rigid_inter_scoring_func(this->startingConformation, transform_dt, this->prot) -
                          score_;
            }
            /*
            BOOST_LOG_TRIVIAL(debug) << "Transform: " << x.t();
            BOOST_LOG_TRIVIAL(debug) << "Score: " << score_;
            BOOST_LOG_TRIVIAL(debug) << "Gradient";
            BOOST_LOG_TRIVIAL(debug) << "     ds: " << grad[3];
            BOOST_LOG_TRIVIAL(debug) << "     du: " << grad[4] << "   dx: " << grad[0];
            BOOST_LOG_TRIVIAL(debug) << "     dv: " << grad[5] << "   dy: " << grad[1];
            BOOST_LOG_TRIVIAL(debug) << "     dt: " << grad[6] << "   dx: " << grad[2];
            //*/

            assert(score_ == score_); // catches NaN
            return score_;
        }


        VinaLikeRigidScoringFunction::VinaLikeRigidScoringFunction(const iConformer &startingConformation_,
                                                                   const iProtein &p,
                                                                   const iTransform &initialTransform_,
                                                                   double differential_epsilon_) :

                startingConformation(startingConformation_),
                prot(p),
                initialTransform(initialTransform_),
                differential_epsilon(differential_epsilon_) {}


        double VinaLikeRigidScoringFunction::getDifferentialEpsilon() const {
            return this->differential_epsilon;
        }

        arma::mat VinaLikeRigidScoringFunction::getStartingConditions() const {
            return this->externalToInternalRepr(this->initialTransform);
        }

        iConformer VinaLikeRigidScoringFunction::getConformerForParamMatrix(const arma::mat &x) {
            assert(x.n_rows == 7);

            iTransform tr = this->internalToExternalRepr(x);
            normalizeQuaternionInPlace(
                    tr.rota); //!< Note that the internal representation arma::mat is not by itself normalized
            //! (because we always normalize it in the scoring function, so no constraint)

            iConformer ret = this->startingConformation;
            applyRigidTransformInPlace(ret, tr);

            return ret;
        }

        unsigned int VinaLikeRigidScoringFunction::getParamVectorDimension() const {
            return 7;
        }


    }
}