//
// Created by eliane on 04/03/19.
//

#ifndef SMOLDOCK_VINALIKESCORINGFUNCTION_H
#define SMOLDOCK_VINALIKESCORINGFUNCTION_H

#include <Engines/Internals/iConformer.h>
#include <Engines/Internals/iProtein.h>
#include <Engines/Internals/iTransform.h>
#include "ScoringFunctionInterface.h"

#include <Structures/Molecule.h>
#include <Structures/Protein.h>

namespace SmolDock::Score {


    double VinaLikeIntermolecularScoringFunction(const iConformer &conformer, const iTransform &transform,
                                                 const iProtein &protein);

    class VinaLikeScoringFunction : public ScoringFunction {
    public:
        VinaLikeScoringFunction(const iConformer &startingConformation_,
                                const iProtein &p,
                                const iTransform &initialTransform_,
                                double differential_epsilon_ = 1e-3);


        double Evaluate(const arma::mat &x) final;

        double EvaluateWithGradient(const arma::mat &x, arma::mat &gradient) final;

        arma::mat getStartingConditions() const final;


        double getDifferentialEpsilon() const final;


        iConformer getConformerForParamMatrix(const arma::mat &x) final;

        unsigned int getParamVectorDimension() const final;


        ~VinaLikeScoringFunction() final = default;


    private:

        inline iTransform internalToExternalRepr(const arma::mat &x_) const {
            assert(x_.n_rows == this->numberOfParamInState);

            iTransform tr_ = iTransform();

            tr_.transl.x() = x_[0];
            tr_.transl.y() = x_[1];
            tr_.transl.z() = x_[2];

            tr_.rota.w() = x_[3];
            tr_.rota.x() = x_[4];
            tr_.rota.y() = x_[5];
            tr_.rota.z() = x_[6];

            for (unsigned int i = 0; i < this->numberOfRotatableBonds; i++) {
                tr_.bondRotationsAngles.push_back(x_[7 + i]);
            }

            return tr_;
        }

        inline arma::mat externalToInternalRepr(const iTransform &tr_) const {
            arma::mat ret(this->numberOfParamInState, 1);

            ret[0] = tr_.transl.x();
            ret[1] = tr_.transl.y();
            ret[2] = tr_.transl.z();

            ret[3] = tr_.rota.w();
            ret[4] = tr_.rota.x();
            ret[5] = tr_.rota.y();
            ret[6] = tr_.rota.z();

            for (unsigned int i = 0; i < this->numberOfRotatableBonds; i++) {
                ret[7 + i] = tr_.bondRotationsAngles[i];
            }

            return ret;
        }

        iConformer startingConformation;
        const iProtein &prot;
        const iTransform initialTransform;
        double differential_epsilon;

        unsigned int numberOfParamInState;
        unsigned int numberOfRotatableBonds;
    };
}


#endif //SMOLDOCK_VINALIKESCORINGFUNCTION_H
