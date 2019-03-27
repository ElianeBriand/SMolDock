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

#include "VinaLikeRigid.h"

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

        const std::array<std::string, VinaLikeRigid::numCoefficients>
                VinaLikeRigid::coefficientsNames =  {"Gauss1", "Gauss2", "RepulsionExceptCovalent", "Hydrophobic","Hydrogen"};


        double vina_like_rigid_inter_scoring_func(const iConformer &ligand, iTransform &transform,
                                                  const iProtein &protein) {

            assert(!ligand.x.empty());
            assert(!protein.x.empty());

            if(std::abs(transform.rota.norm() - 1) > 0.1) {
                transform.rota.normalize();
            }

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

        double VinaLikeRigid::Evaluate(const arma::mat &x) {
            assert(x.n_rows == 7);

            iTransform tr = this->internalToExternalRepr(x);
            normalizeQuaternionInPlace(tr.rota);

            double score_ = vina_like_rigid_inter_scoring_func(this->startingConformation, tr, this->prot);

            return score_;
        }

        double VinaLikeRigid::EvaluateWithGradient(const arma::mat &x, arma::mat &grad) {

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


        VinaLikeRigid::VinaLikeRigid(const iConformer &startingConformation_,
                                                                   const iProtein &p,
                                                                   const iTransform &initialTransform_,
                                                                   double differential_epsilon_) :

                startingConformation(startingConformation_),
                prot(p),
                initialTransform(initialTransform_),
                differential_epsilon(differential_epsilon_) {}


        double VinaLikeRigid::getDifferentialEpsilon() const {
            return this->differential_epsilon;
        }

        arma::mat VinaLikeRigid::getStartingConditions() const {
            return this->externalToInternalRepr(this->initialTransform);
        }

        iConformer VinaLikeRigid::getConformerForParamMatrix(const arma::mat &x) {
            assert(x.n_rows == 7);

            iTransform tr = this->internalToExternalRepr(x);
            normalizeQuaternionInPlace(
                    tr.rota); //!< Note that the internal representation arma::mat is not by itself normalized
            //! (because we always normalize it in the scoring function, so no constraint)

            iConformer ret = this->startingConformation;
            applyRigidTransformInPlace(ret, tr);

            return ret;
        }

        unsigned int VinaLikeRigid::getParamVectorDimension() const {
            return 7;
        }


        std::vector<std::tuple<std::string, double>> VinaLikeRigid::EvaluateSubcomponents(const arma::mat &x) {
            std::vector<std::tuple<std::string, double>> ret;

            assert(x.n_rows == 7);

            iTransform tr = this->internalToExternalRepr(x);
            normalizeQuaternionInPlace(tr.rota);

            assert(!this->startingConformation.x.empty());
            assert(!this->prot.x.empty());

            if(std::abs(tr.rota.norm() - 1) > 0.1) {
                tr.rota.normalize();
            }

            assert(tr.bondRotationsAngles.size() == this->startingConformation.num_rotatable_bond);


            double gauss1_total = 0.0;
            double gauss2_total = 0.0;
            double repulsion_total = 0.0;
            double hydrogen_total = 0.0;
            double hydrophobic_total = 0.0;
            double score_raw = 0.0;

            iConformer ligand = this->startingConformation;
            applyBondRotationInPlace(ligand, tr);

            Eigen::Vector3d ProtCenterPosition = {this->prot.center_x, this->prot.center_y, this->prot.center_z};

            for (unsigned int idxLig = 0; idxLig < ligand.x.size(); idxLig++) {
                for (unsigned int idxProt = 0; idxProt < this->prot.x.size(); idxProt++) {


                    Eigen::Vector3d LigPosition = {ligand.x[idxLig], ligand.y[idxLig], ligand.z[idxLig]};
                    applyRigidTransformInPlace(LigPosition, tr);

                    Eigen::Vector3d ProtPosition = {this->prot.x[idxProt], this->prot.y[idxProt], this->prot.z[idxProt]};
                    Eigen::Vector3d distToCenterVector = LigPosition - ProtCenterPosition;

                    double distanceToProteinCenter = distToCenterVector.norm();

                    if (distanceToProteinCenter > (this->prot.radius - 1)) {
                        score_raw += std::pow((distanceToProteinCenter - this->prot.radius), 4) + 10;
                        continue;
                    }

                    Eigen::Vector3d distVect = ProtPosition - LigPosition;

                    double rawDist = distVect.norm();

                    const double cutoff = 8.0;
                    if (rawDist >= cutoff)
                        continue;

                    double distance = distanceFromRawDistance(rawDist,  ligand.atomicRadius[idxLig], this->prot.atomicRadius[idxProt]);

                    unsigned char atom1AtomicNumber = ligand.type[idxLig];
                    unsigned int atom1AtomVariant = ligand.variant[idxLig];
                    unsigned char atom2AtomicNumber = this->prot.type[idxProt];
                    unsigned int atom2AtomVariant = this->prot.variant[idxProt];

                    gauss1_total += vinaGaussComponent(distance, 0.0, 0.5);
                    gauss2_total += vinaGaussComponent(distance, 3.0, 2.0);
                    repulsion_total += vinaRepulsionComponent(distance, 0.0);
                    hydrogen_total += vinaHydrophobicComponent(distance,
                                                               atom1AtomicNumber, atom1AtomVariant,
                                                               atom2AtomicNumber, atom2AtomVariant);
                    hydrophobic_total += vinaHydrogenComponent(distance,
                                                               atom1AtomicNumber, atom1AtomVariant,
                                                               atom2AtomicNumber, atom2AtomVariant);
                    score_raw += scoreForAtomCouple(distance,
                                                    atom1AtomicNumber, atom1AtomVariant,
                                                    atom2AtomicNumber, atom2AtomVariant);

                } // for
            } // for

            double final_score = score_raw / (1 + (0.058459999999999998 * ligand.num_rotatable_bond));


            ret.emplace_back(std::make_tuple("Gauss1",gauss1_total));
            ret.emplace_back(std::make_tuple("Gauss2",gauss2_total));
            ret.emplace_back(std::make_tuple("Repulsion",repulsion_total));
            ret.emplace_back(std::make_tuple("Hydrophobic",hydrogen_total));
            ret.emplace_back(std::make_tuple("Hydrogen",hydrophobic_total));
            ret.emplace_back(std::make_tuple("ScoreRaw",score_raw));
            ret.emplace_back(std::make_tuple("Score",final_score));
            return ret;
        }

        double VinaLikeRigid::EvaluateOnlyIntermolecular(const arma::mat &x) {
            return this->Evaluate(x);
        }

        unsigned int VinaLikeRigid::getCoefficientsVectorWidth() {
            return this->numCoefficients;
        }

        std::vector<std::string> VinaLikeRigid::getCoefficientsNames() {
            return std::vector<std::string>(this->coefficientsNames.begin(),this->coefficientsNames.end());
        }

        std::vector<double> VinaLikeRigid::getCurrentCoefficients() {
            return {VinaClassic::coeff_gauss1, VinaClassic::coeff_gauss2,
                    VinaClassic::coeff_repulsion, VinaClassic::coeff_hydrophobic,
                    VinaClassic::coeff_hydrogen };
        }

        bool VinaLikeRigid::setNonDefaultCoefficients(std::vector<double> coeffs) {
            BOOST_LOG_TRIVIAL(debug) << "Non default coefficients not supported yet on Vina rigid scoring function.";
            return false;
        }


    }




}