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

    const std::array<std::string, VinaLike::numCoefficients>
            VinaLike::coefficientsNames =  {"Gauss1", "Gauss2", "RepulsionExceptCovalent", "Hydrophobic","Hydrogen"};

    template<bool OnlyIntermolecular, bool useNonDefaultCoefficients>
    double VinaLikeIntermolecularScoringFunction(const iConformer &ligand_, const iTransform &transform,
                                                 const iProtein &protein, std::array<double, VinaLike_numCoefficients> nonDefaultCoeffs) {

        assert(!ligand_.x.empty());
        assert(!protein.x.empty());
        assert(std::abs(transform.rota.norm() - 1) < 0.1);
        assert(transform.bondRotationsAngles.size() == ligand_.num_rotatable_bond);

        double score_raw = 0;

        iConformer ligand = ligand_;
        applyBondRotationInPlace(ligand, transform);

        Eigen::Vector3d ProtCenterPosition = {protein.center_x, protein.center_y, protein.center_z};

        if constexpr(!OnlyIntermolecular) // C++17
        {
            for (unsigned int idxLig = 0; idxLig < ligand.x.size(); idxLig++) {
                for (unsigned int idxLig2 = idxLig; idxLig2 < ligand.x.size(); idxLig2++) {
                    Eigen::Vector3d LigDistance = {ligand.x[idxLig] - ligand.x[idxLig2],
                                                   ligand.y[idxLig] - ligand.y[idxLig2],
                                                   ligand.z[idxLig] - ligand.z[idxLig2]};

                    double distance_raw = LigDistance.norm();
                    const double distance = distanceFromRawDistance(distance_raw, ligand.atomicRadius[idxLig],
                                                                    ligand.atomicRadius[idxLig2]);

//                score_raw += VinaClassic::coeff_gauss1      * vinaGaussComponent(distance, 0.0, 0.5);
//                score_raw += VinaClassic::coeff_gauss2      * vinaGaussComponent(distance, 3.0, 2.0);
                    score_raw += VinaClassic::coeff_repulsion   * vinaRepulsionComponent(distance, 0.0);
                }
            }
        }


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

                if constexpr(useNonDefaultCoefficients)
                {
                    score_raw += nonDefaultCoeffs[0]      * vinaGaussComponent(distance, 0.0, 0.5);
                    score_raw += nonDefaultCoeffs[1]     * vinaGaussComponent(distance, 3.0, 2.0);
                    score_raw += nonDefaultCoeffs[2]   * vinaRepulsionComponent(distance, 0.0);
                    score_raw += nonDefaultCoeffs[3] * vinaHydrophobicComponent(distance,
                                                                                           atom1AtomicNumber, atom1AtomVariant,
                                                                                           atom2AtomicNumber, atom2AtomVariant);

                    score_raw += nonDefaultCoeffs[4]    * vinaHydrogenComponent(distance,
                                                                                        atom1AtomicNumber, atom1AtomVariant,
                                                                                        atom2AtomicNumber, atom2AtomVariant);
                }else {
                    score_raw += VinaClassic::coeff_gauss1      * vinaGaussComponent(distance, 0.0, 0.5);
                    score_raw += VinaClassic::coeff_gauss2      * vinaGaussComponent(distance, 3.0, 2.0);
                    score_raw += VinaClassic::coeff_repulsion   * vinaRepulsionComponent(distance, 0.0);
                    score_raw += VinaClassic::coeff_hydrophobic * vinaHydrophobicComponent(distance,
                                                                                           atom1AtomicNumber, atom1AtomVariant,
                                                                                           atom2AtomicNumber, atom2AtomVariant);

                    score_raw += VinaClassic::coeff_hydrogen    * vinaHydrogenComponent(distance,
                                                                                        atom1AtomicNumber, atom1AtomVariant,
                                                                                        atom2AtomicNumber, atom2AtomVariant);
                }


            } // for
        } // for

        double final_score = score_raw / (1 + (VinaClassic::coeff_entropic * ligand.num_rotatable_bond));
        return final_score;
    }


    VinaLike::VinaLike(const iConformer &startingConformation_,
                                                     const iProtein &p,
                                                     const iTransform &initialTransform_,
                                                     double differential_epsilon_,
                                                     bool useNonDefaultCoefficient)  :
            useNonDefaultCoefficient(useNonDefaultCoefficient),
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

        double score_ = this->useNonDefaultCoefficient ?
                        VinaLikeIntermolecularScoringFunction<false,true>(this->startingConformation, tr, this->prot, this->nonDefaultCoefficients)
                                                       : VinaLikeIntermolecularScoringFunction<false,false>(this->startingConformation, tr, this->prot);

        return score_;
    }

    double VinaLike::EvaluateWithGradient(const arma::mat &x, arma::mat &grad) {

        assert(!x.has_nan());
        assert(!grad.has_nan());
        assert(x.n_rows == this->numberOfParamInState);
        assert(grad.n_rows == this->numberOfParamInState);

        iTransform tr = this->internalToExternalRepr(x);
        normalizeQuaternionInPlace(tr.rota);

        double score_ = this->useNonDefaultCoefficient ?
                          VinaLikeIntermolecularScoringFunction<false,true>(this->startingConformation, tr, this->prot, this->nonDefaultCoefficients)
                        : VinaLikeIntermolecularScoringFunction<false,false>(this->startingConformation, tr, this->prot);


        // Translation
        {
            iTransform transform_dx = tr;
            transform_dx.transl.x() += this->differential_epsilon;
            const double gradScore = this->useNonDefaultCoefficient ?
                              VinaLikeIntermolecularScoringFunction<false,true>(this->startingConformation, transform_dx, this->prot, this->nonDefaultCoefficients)
                                                             : VinaLikeIntermolecularScoringFunction<false,false>(this->startingConformation, transform_dx, this->prot);
            grad[0] = gradScore - score_;
        }

        {
            iTransform transform_dy = tr;
            transform_dy.transl.y() += this->differential_epsilon;
            const double gradScore = this->useNonDefaultCoefficient ?
                               VinaLikeIntermolecularScoringFunction<false,true>(this->startingConformation, transform_dy, this->prot, this->nonDefaultCoefficients)
                                                              : VinaLikeIntermolecularScoringFunction<false,false>(this->startingConformation, transform_dy, this->prot);
            grad[1] = gradScore - score_;
        }

        {
            iTransform transform_dz = tr;
            transform_dz.transl.z() += this->differential_epsilon;
            const double gradScore = this->useNonDefaultCoefficient ?
                               VinaLikeIntermolecularScoringFunction<false,true>(this->startingConformation, transform_dz, this->prot, this->nonDefaultCoefficients)
                                                              : VinaLikeIntermolecularScoringFunction<false,false>(this->startingConformation, transform_dz, this->prot);
            grad[2] = gradScore - score_;
        }

        // Rotation

        {
            iTransform transform_dqs = tr;
            transform_dqs.rota.w() += this->differential_epsilon;
            normalizeQuaternionInPlace(transform_dqs.rota);
            const double gradScore = this->useNonDefaultCoefficient ?
                               VinaLikeIntermolecularScoringFunction<false,true>(this->startingConformation, transform_dqs, this->prot, this->nonDefaultCoefficients)
                                                              : VinaLikeIntermolecularScoringFunction<false,false>(this->startingConformation, transform_dqs, this->prot);
            grad[3] = gradScore - score_;
        }

        {
            iTransform transform_dqx = tr;
            transform_dqx.rota.x() += this->differential_epsilon;
            normalizeQuaternionInPlace(transform_dqx.rota);
            const double gradScore = this->useNonDefaultCoefficient ?
                               VinaLikeIntermolecularScoringFunction<false,true>(this->startingConformation, transform_dqx, this->prot, this->nonDefaultCoefficients)
                                                              : VinaLikeIntermolecularScoringFunction<false,false>(this->startingConformation, transform_dqx, this->prot);
            grad[4] = gradScore - score_;

        }

        {
            iTransform transform_dqy = tr;
            transform_dqy.rota.x() += this->differential_epsilon;
            normalizeQuaternionInPlace(transform_dqy.rota);
            const double gradScore = this->useNonDefaultCoefficient ?
                                     VinaLikeIntermolecularScoringFunction<false,true>(this->startingConformation, transform_dqy, this->prot, this->nonDefaultCoefficients)
                                                                    : VinaLikeIntermolecularScoringFunction<false,false>(this->startingConformation, transform_dqy, this->prot);
            grad[5] = gradScore - score_;
        }

        {
            iTransform transform_dqz = tr;
            transform_dqz.rota.x() += this->differential_epsilon;
            normalizeQuaternionInPlace(transform_dqz.rota);
            const double gradScore = this->useNonDefaultCoefficient ?
                                     VinaLikeIntermolecularScoringFunction<false,true>(this->startingConformation, transform_dqz, this->prot, this->nonDefaultCoefficients)
                                                                    : VinaLikeIntermolecularScoringFunction<false,false>(this->startingConformation, transform_dqz, this->prot);
            grad[6] = gradScore - score_;
        }

        for (unsigned int i = 0; i < this->numberOfRotatableBonds; i++) {
            iTransform transform_dbondrot = tr;
            transform_dbondrot.bondRotationsAngles[i] += this->differential_epsilon;
            const double gradScore = this->useNonDefaultCoefficient ?
                                     VinaLikeIntermolecularScoringFunction<false,true>(this->startingConformation, transform_dbondrot, this->prot, this->nonDefaultCoefficients)
                                                                    : VinaLikeIntermolecularScoringFunction<false,false>(this->startingConformation, transform_dbondrot, this->prot);
            grad[7 + i] = gradScore - score_;

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

        double intramolecular_repuls_total = 0.0;
        double intramolecular_score = 0.0;

        iConformer ligand = this->startingConformation;
        applyBondRotationInPlace(ligand, tr);

        Eigen::Vector3d ProtCenterPosition = {this->prot.center_x, this->prot.center_y, this->prot.center_z};


        for (unsigned int idxLig = 0; idxLig < ligand.x.size(); idxLig++) {
            for (unsigned int idxLig2 = idxLig; idxLig2 < ligand.x.size(); idxLig2++) {

                if(idxLig == idxLig2)
                    continue;

                Eigen::Vector3d LigDistance = {ligand.x[idxLig] - ligand.x[idxLig2],
                                               ligand.y[idxLig] - ligand.y[idxLig2],
                                               ligand.z[idxLig] - ligand.z[idxLig2]};

                double distance_raw = LigDistance.norm();
                const double distance = distanceFromRawDistance(distance_raw, ligand.atomicRadius[idxLig],
                                                                ligand.atomicRadius[idxLig2]);

//                score_raw += VinaClassic::coeff_gauss1      * vinaGaussComponent(distance, 0.0, 0.5);
//                score_raw += VinaClassic::coeff_gauss2      * vinaGaussComponent(distance, 3.0, 2.0);
                intramolecular_repuls_total +=  vinaRepulsionComponent(distance, 0.0);

            }
        }

        intramolecular_score = VinaClassic::coeff_repulsion   * intramolecular_repuls_total;


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


        double score_sum = 0.0;

        if(this->useNonDefaultCoefficient)
        {
            ret.emplace_back(std::make_tuple("NonDefaultCoeffs", 1.0));
            score_sum =   this->nonDefaultCoefficients[0] * gauss1_total
                          + this->nonDefaultCoefficients[1] * gauss2_total
                          + this->nonDefaultCoefficients[2] * repulsion_total
                          + this->nonDefaultCoefficients[3] * hydrophobic_total
                          + this->nonDefaultCoefficients[4] * hydrogen_total;
        }else {
            score_sum =   VinaClassic::coeff_gauss1 * gauss1_total
                          + VinaClassic::coeff_gauss2 * gauss2_total
                          + VinaClassic::coeff_repulsion * repulsion_total
                          + VinaClassic::coeff_hydrophobic * hydrophobic_total
                          + VinaClassic::coeff_hydrogen * hydrogen_total;
        }


        double final_score_fromSum = score_sum / (1 + (VinaClassic::coeff_entropic * ligand.num_rotatable_bond));

        ret.emplace_back(std::make_tuple("Gauss1", gauss1_total));
        ret.emplace_back(std::make_tuple("Gauss2", gauss2_total));
        ret.emplace_back(std::make_tuple("Repulsion", repulsion_total));
        ret.emplace_back(std::make_tuple("Hydrophobic", hydrophobic_total));
        ret.emplace_back(std::make_tuple("Hydrogen", hydrogen_total));
        ret.emplace_back(std::make_tuple("Intra_Repuls", intramolecular_repuls_total));
        ret.emplace_back(std::make_tuple("Intra_Score", intramolecular_score));
        ret.emplace_back(std::make_tuple("ScoreRawSum", score_sum));
        ret.emplace_back(std::make_tuple("Score", final_score_fromSum));
        return ret;
    }

    double VinaLike::EvaluateOnlyIntermolecular(const arma::mat &x) {
        assert(x.n_rows == this->numberOfParamInState);

        iTransform tr = this->internalToExternalRepr(x);

        normalizeQuaternionInPlace(tr.rota);

        // Template parameter controls whether onlyIntermolecular interaction are taken into account. Here we want true
        double score_ = this->useNonDefaultCoefficient ?
                        VinaLikeIntermolecularScoringFunction<false,true>(this->startingConformation, tr, this->prot, this->nonDefaultCoefficients)
                                                       : VinaLikeIntermolecularScoringFunction<false,false>(this->startingConformation, tr, this->prot);


        return score_;
    }

    unsigned int VinaLike::getCoefficientsVectorWidth() {
        return this->numCoefficients;
    }

    std::vector<std::string> VinaLike::getCoefficientsNames() {
        return std::vector<std::string>(this->coefficientsNames.begin(),this->coefficientsNames.end());
    }

    std::vector<double> VinaLike::getCurrentCoefficients() {
        if(this->useNonDefaultCoefficient)
        {
            return std::vector<double>(this->nonDefaultCoefficients.begin(),this->nonDefaultCoefficients.end());
        }
        return {VinaClassic::coeff_gauss1, VinaClassic::coeff_gauss2,
                VinaClassic::coeff_repulsion, VinaClassic::coeff_hydrophobic,
                VinaClassic::coeff_hydrogen};
    }

    bool VinaLike::setNonDefaultCoefficients(std::vector<double> coeffs) {
        if(coeffs.size() != this->numCoefficients)
        {
            BOOST_LOG_TRIVIAL(error) << "Trying to set " << this->numCoefficients <<" coefficients with vector of " << coeffs.size() << " values.";
            return false;
        }
        if(this->useNonDefaultCoefficient == false)
        {
            BOOST_LOG_TRIVIAL(error) << "Trying to set non default coefficient, but this scoring function was constructed with default coefficients only.";
            BOOST_LOG_TRIVIAL(error) << "Check the parameters passed to the scoring function constructor.";
            return false;
        }
        for (unsigned int j = 0; j < this->nonDefaultCoefficients.size(); ++j) {
            this->nonDefaultCoefficients[j] = coeffs[j];
        }
        return true;
    }
}