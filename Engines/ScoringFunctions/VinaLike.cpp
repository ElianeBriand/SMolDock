//
// Created by eliane on 04/03/19.
//

#include "VinaLike.h"

#include <exception>

#include <Structures/Atom.h>

#include <Engines/Internals/InternalsUtilityFunctions.h>

#include <boost/log/trivial.hpp>

#include "VinaLikeCommon.h"

namespace SmolDock::Score {


    double VinaLikeIntermolecularScoringFunction(const iConformer &ligand_, const iTransform &transform,
                                                 const iProtein &protein) {

        assert(!ligand_.x.empty());
        assert(!protein.x.empty());
        assert(std::abs(transform.rota.norm() - 1) < 0.01);
        assert(transform.bondRotationsAngles.size() == ligand_.num_rotatable_bond);

        double score_raw = 0;

        iConformer ligand = ligand_;
        applyBondRotationInPlace(ligand, transform);

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

                const double rawDist = distVect.norm();


                if (rawDist >= VinaClassic::interaction_cutoff)
                    continue;


                const double atomicRadiusLig = ligand.atomicRadius[idxLig];
                const double atomicRadiusProt = protein.atomicRadius[idxProt];

                const double radToRemove = (atomicRadiusLig + atomicRadiusProt);

                const double distance = rawDist - radToRemove;

                const unsigned char atom1AtomicNumber = ligand.type[idxLig];
                const unsigned int atom1AtomVariant = ligand.variant[idxLig];
                const unsigned char atom2AtomicNumber = protein.type[idxProt];
                const unsigned int atom2AtomVariant = protein.variant[idxProt];


                score_raw += VinaClassic::coeff_gauss1      * vinaGaussComponent(distance, 0.0, 0.5);
                score_raw += VinaClassic::coeff_gauss2      * vinaGaussComponent(distance, 3.0, 2.0);
                score_raw += VinaClassic::coeff_repulsion   * vinaRepulsionComponent(distance, 0.0);
                score_raw += VinaClassic::coeff_hydrophobic * vinaHydrophobicComponent(distance,
                                                      atom1AtomicNumber, atom1AtomVariant,
                                                      atom2AtomicNumber, atom2AtomVariant);

                score_raw += VinaClassic::coeff_hydrogen    * vinaHydrogenComponent(distance,
                                                   atom1AtomicNumber, atom1AtomVariant,
                                                   atom2AtomicNumber, atom2AtomVariant);

            } // for
        } // for

        double final_score = score_raw / (1 + (VinaClassic::coeff_entropic * ligand.num_rotatable_bond));
        return final_score;
    }


    VinaLike::VinaLike(const iConformer &startingConformation_,
                                                     const iProtein &p,
                                                     const iTransform &initialTransform_,
                                                     double differential_epsilon_) :

            startingConformation(startingConformation_),
            prot(p),
            initialTransform(initialTransform_),
            differential_epsilon(differential_epsilon_) {
        this->numberOfRotatableBonds = this->startingConformation.num_rotatable_bond;
        this->numberOfParamInState = 7 + (this->numberOfRotatableBonds);

        if (this->initialTransform.bondRotationsAngles.size() != this->numberOfRotatableBonds) {
            BOOST_LOG_TRIVIAL(error)
                << "Discrepency between the number of rotatable bonds in the iConformer and iTransform ("
                << this->numberOfRotatableBonds << " != " << this->initialTransform.bondRotationsAngles.size() << ")";
            std::terminate();
        }

    }


    double VinaLike::Evaluate(const arma::mat &x) {
        assert(x.n_rows == this->numberOfParamInState);

        iTransform tr = this->internalToExternalRepr(x);

        normalizeQuaternionInPlace(tr.rota);

        double score_ = VinaLikeIntermolecularScoringFunction(this->startingConformation, tr, this->prot);

        return score_;
    }

    double VinaLike::EvaluateWithGradient(const arma::mat &x, arma::mat &grad) {

        assert(!x.has_nan());
        assert(!grad.has_nan());
        assert(x.n_rows == this->numberOfParamInState);
        assert(grad.n_rows == this->numberOfParamInState);

        iTransform tr = this->internalToExternalRepr(x);
        normalizeQuaternionInPlace(tr.rota);


        double score_ = VinaLikeIntermolecularScoringFunction(this->startingConformation, tr, this->prot);

        // Translation
        {
            iTransform transform_dx = tr;
            transform_dx.transl.x() += this->differential_epsilon;
            grad[0] = VinaLikeIntermolecularScoringFunction(this->startingConformation, transform_dx, this->prot) -
                      score_;
        }

        {
            iTransform transform_dy = tr;
            transform_dy.transl.y() += this->differential_epsilon;
            grad[1] = VinaLikeIntermolecularScoringFunction(this->startingConformation, transform_dy, this->prot) -
                      score_;
        }

        {
            iTransform transform_dz = tr;
            transform_dz.transl.z() += this->differential_epsilon;
            grad[2] = VinaLikeIntermolecularScoringFunction(this->startingConformation, transform_dz, this->prot) -
                      score_;
        }

        // Rotation

        {
            iTransform transform_ds = tr;
            transform_ds.rota.w() += this->differential_epsilon;
            normalizeQuaternionInPlace(transform_ds.rota);
            grad[3] = VinaLikeIntermolecularScoringFunction(this->startingConformation, transform_ds, this->prot) -
                      score_;
        }

        {
            iTransform transform_du = tr;
            transform_du.rota.x() += this->differential_epsilon;
            normalizeQuaternionInPlace(transform_du.rota);
            grad[4] = VinaLikeIntermolecularScoringFunction(this->startingConformation, transform_du, this->prot) -
                      score_;
        }

        {
            iTransform transform_dv = tr;
            transform_dv.rota.y() += this->differential_epsilon;
            normalizeQuaternionInPlace(transform_dv.rota);
            grad[5] = VinaLikeIntermolecularScoringFunction(this->startingConformation, transform_dv, this->prot) -
                      score_;
        }

        {
            iTransform transform_dt = tr;
            transform_dt.rota.z() += this->differential_epsilon;
            normalizeQuaternionInPlace(transform_dt.rota);
            grad[6] = VinaLikeIntermolecularScoringFunction(this->startingConformation, transform_dt, this->prot) -
                      score_;
        }

        for (unsigned int i = 0; i < this->numberOfRotatableBonds; i++) {
            iTransform transform_dbondrot = tr;
            transform_dbondrot.bondRotationsAngles[i] += this->differential_epsilon;
            grad[6 + i] =
                    VinaLikeIntermolecularScoringFunction(this->startingConformation, transform_dbondrot, this->prot) -
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


    double VinaLike::getDifferentialEpsilon() const {
        return this->differential_epsilon;
    }

    arma::mat VinaLike::getStartingConditions() const {
        return this->externalToInternalRepr(this->initialTransform);
    }

    iConformer VinaLike::getConformerForParamMatrix(const arma::mat &x) {
        assert(x.n_rows == this->numberOfParamInState);

        iTransform tr = this->internalToExternalRepr(x);
        tr.rota.normalize();

        iConformer ret = this->startingConformation;
        applyBondRotationInPlace(ret, tr);
        applyRigidTransformInPlace(ret, tr);

        return ret;
    }

    unsigned int VinaLike::getParamVectorDimension() const {
        return this->numberOfParamInState;
    }

    std::vector<std::tuple<std::string, double>> VinaLike::EvaluateSubcomponents(const arma::mat &x) {
        std::vector<std::tuple<std::string, double>> ret;

        assert(x.n_rows == this->numberOfParamInState);

        iTransform tr = this->internalToExternalRepr(x);
        normalizeQuaternionInPlace(tr.rota);

        assert(!this->startingConformation.x.empty());
        assert(!this->prot.x.empty());
        assert(std::abs(tr.rota.norm() - 1) < 0.01);
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

                if (rawDist >= VinaClassic::interaction_cutoff)
                    continue;

                double distance = distanceFromRawDistance(rawDist,  ligand.atomicRadius[idxLig], this->prot.atomicRadius[idxProt]);

                const unsigned char atom1AtomicNumber = ligand.type[idxLig];
                const unsigned int atom1AtomVariant = ligand.variant[idxLig];
                const unsigned char atom2AtomicNumber = this->prot.type[idxProt];
                const unsigned int atom2AtomVariant = this->prot.variant[idxProt];

                gauss1_total += vinaGaussComponent(distance, 0.0, 0.5);
                gauss2_total += vinaGaussComponent(distance, 3.0, 2.0);
                repulsion_total += vinaRepulsionComponent(distance, 0.0);
                hydrophobic_total+= vinaHydrophobicComponent(distance,
                                                           atom1AtomicNumber, atom1AtomVariant,
                                                           atom2AtomicNumber, atom2AtomVariant);
                hydrogen_total += vinaHydrogenComponent(distance,
                                                           atom1AtomicNumber, atom1AtomVariant,
                                                           atom2AtomicNumber, atom2AtomVariant);

            } // for
        } // for

        double score_sum =   VinaClassic::coeff_gauss1 * gauss1_total
                      + VinaClassic::coeff_gauss2 * gauss2_total
                      + VinaClassic::coeff_repulsion * repulsion_total
                      + VinaClassic::coeff_hydrophobic * hydrophobic_total
                      + VinaClassic::coeff_hydrogen * hydrogen_total;


        double final_score_fromSum = score_sum / (1 + (VinaClassic::coeff_entropic * ligand.num_rotatable_bond));

        ret.push_back(std::make_tuple("Gauss1", gauss1_total));
        ret.push_back(std::make_tuple("Gauss2", gauss2_total));
        ret.push_back(std::make_tuple("Repulsion", repulsion_total));
        ret.push_back(std::make_tuple("Hydrophobic", hydrophobic_total));
        ret.push_back(std::make_tuple("Hydrogen", hydrogen_total));
        ret.push_back(std::make_tuple("ScoreRawSum", score_sum));
        ret.push_back(std::make_tuple("Score", final_score_fromSum));
        return ret;
    }
}