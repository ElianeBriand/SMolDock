//
// Created by eliane on 04/03/19.
//

#include "VinaLikeScoringFunction.h"

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

                double rawDist = distVect.norm();

                const double cutoff = 8.0;
                if (rawDist >= cutoff)
                    continue;


                double atomicRadiusLig = ligand.atomicRadius[idxLig];
                double atomicRadiusProt = protein.atomicRadius[idxProt];

                double radToRemove = (atomicRadiusLig + atomicRadiusProt);

                double distance = rawDist - radToRemove;

                score_raw += scoreForAtomCouple(distance, ligand.type[idxLig], ligand.variant[idxLig],
                                                protein.type[idxProt], protein.variant[idxProt]);

            } // for
        } // for

        double final_score = score_raw / (1 + (0.058459999999999998 * ligand.num_rotatable_bond));
        return final_score;
    }


    VinaLikeScoringFunction::VinaLikeScoringFunction(const iConformer &startingConformation_,
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


    double VinaLikeScoringFunction::Evaluate(const arma::mat &x) {
        assert(x.n_rows == this->numberOfParamInState);

        iTransform tr = this->internalToExternalRepr(x);

        normalizeQuaternionInPlace(tr.rota);

        double score_ = VinaLikeIntermolecularScoringFunction(this->startingConformation, tr, this->prot);

        return score_;
    }

    double VinaLikeScoringFunction::EvaluateWithGradient(const arma::mat &x, arma::mat &grad) {

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


    double VinaLikeScoringFunction::getDifferentialEpsilon() const {
        return this->differential_epsilon;
    }

    arma::mat VinaLikeScoringFunction::getStartingConditions() const {
        return this->externalToInternalRepr(this->initialTransform);
    }

    iConformer VinaLikeScoringFunction::getConformerForParamMatrix(const arma::mat &x) {
        assert(x.n_rows == this->numberOfParamInState);

        iTransform tr = this->internalToExternalRepr(x);
        tr.rota.normalize();

        iConformer ret = this->startingConformation;
        applyBondRotationInPlace(ret, tr);
        applyRigidTransformInPlace(ret, tr);

        return ret;
    }

    unsigned int VinaLikeScoringFunction::getParamVectorDimension() const {
        return this->numberOfParamInState;
    }
}