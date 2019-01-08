//
// Created by eliane on 03/01/19.
//

#ifndef SMOLDOCK_INTERNALSUTILITYFUNCTIONS_H
#define SMOLDOCK_INTERNALSUTILITYFUNCTIONS_H

#include "iConformer.h"
#include "iTransform.h"

namespace SmolDock {

    inline iQuaternion iQuaternionIdentityInit() {
        iQuaternion qt;
        qt.s = 1.0;
        qt.u = 0.0;
        qt.v = 0.0;
        qt.t = 0.0;
        return qt;
    }

    //! Holds a gradient for the score function with regard to an iTransform


    //! Returns a neutral transform
    inline iTransform iTransformIdentityInit() {
        iTransform tr;
        tr.transl.x = 0.0;
        tr.transl.y = 0.0;
        tr.transl.z = 0.0;
        tr.rota = iQuaternionIdentityInit();
        return tr;
    }

    //! Apply the rotation given by a quaternion to a vector
    inline void applyRotationInPlace(iVect &vec, const iQuaternion &qt) {
        iVect v = vec;

        //TODO : Check this quaternion rotation formula
        //TODO : Refactor this formula to be clearer
        vec.x = qt.s * qt.s * v.x + 2 * qt.v * qt.s * v.z - 2 * qt.t * qt.s * v.y + qt.u * qt.u * v.x +
                2 * qt.v * qt.u * v.y + 2 * qt.t * qt.u * v.z - qt.t * qt.t * v.x - qt.v * qt.v * v.x;
        vec.y = 2 * qt.u * qt.v * v.x + qt.v * qt.v * v.y + 2 * qt.t * qt.v * v.z + 2 * qt.s * qt.t * v.x -
                qt.t * qt.t * v.y + qt.s * qt.s * v.y - 2 * qt.u * qt.s * v.z - qt.u * qt.u * v.y;
        vec.z = 2 * qt.u * qt.t * v.x + 2 * qt.v * qt.t * v.y + qt.t * qt.t * v.z - 2 * qt.s * qt.v * v.x -
                qt.v * qt.v * v.z + 2 * qt.s * qt.u * v.y - qt.u * qt.u * v.z + qt.s * qt.s * v.z;
    }

    //! Apply the rotation given by a quaternion to a vector
    inline iVect applyRotation(const iVect &v, const iQuaternion &qt) {
        iVect res(v);
        applyRotationInPlace(res, qt);
        return res;
    }

    //! Apply the provided translation to a vector
    inline void applyTranslationInPlace(iVect &v, const iTranslation &t) {
        v.x += t.x;
        v.y += t.y;
        v.z += t.z;
    }

    //! Apply the provided translation to a vector
    inline iVect applyTranslation(const iVect &v, const iTranslation &t) {
        iVect res(v);
        applyTranslationInPlace(res, t);
        return res;
    }


    //! Apply the provided iTransform to a vector
    inline void applyTransformInPlace(iVect &v, const iTransform &tr) {
        applyRotationInPlace(v, tr.rota);
        applyTranslationInPlace(v, tr.transl);
    }

    //! Apply the provided iTransform to a vector
    inline iVect applyTransform(const iVect &v, const iTransform &tr) {
        iVect res(v);
        applyTransformInPlace(res, tr);
        return res;
    }


    //! Apply the provided iTransform to a confomer, atom-by-atom
    inline void applyTransformInPlace(SmolDock::iConformer &conformer, const iTransform &tr) {
        for (unsigned int i = 0; i < conformer.x.size(); i++) {
            iVect posVect = {conformer.x[i], conformer.y[i], conformer.z[i]};
            applyTransformInPlace(posVect, tr);
            conformer.x[i] = posVect.x;
            conformer.y[i] = posVect.y;
            conformer.z[i] = posVect.z;
        }
    }

    //! Returns the norm of the given quaternion
    inline double quaternionNorm(const iQuaternion &quat) {
        return std::sqrt(
                std::pow(quat.s, 2) + std::pow(quat.u, 2) + std::pow(quat.v, 2) + std::pow(quat.t, 2)
        );
    }

    //! Returns the norm of the given vector
    inline double vectorNorm(iVect v) {
        return std::sqrt(
                std::pow(v.x, 2) + std::pow(v.y, 2) + std::pow(v.z, 2)
        );
    }

    //! Normalize a quaternion (which then has norm == 1.0)
    inline void normalizeQuaternionInPlace(iQuaternion &quat) {
        double norm = quaternionNorm(quat);
        assert(norm == norm); // This catches NaN
        if (norm == 0.0)
            return;

        quat.s /= norm;
        quat.u /= norm;
        quat.v /= norm;
        quat.t /= norm;
    }

    //! Normalize a quaternion (which then has norm == 1.0)
    inline iQuaternion normalizeQuaternion(const iQuaternion &quat) {
        iQuaternion ret = quat;
        normalizeQuaternionInPlace(ret);
        return ret;
    }


}

#endif //SMOLDOCK_INTERNALSUTILITYFUNCTIONS_H
