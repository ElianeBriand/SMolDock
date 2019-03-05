//
// Created by eliane on 03/01/19.
//

#ifndef SMOLDOCK_INTERNALSUTILITYFUNCTIONS_H
#define SMOLDOCK_INTERNALSUTILITYFUNCTIONS_H


#include <iostream>

#include <cmath>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "iConformer.h"
#include "iTransform.h"


namespace SmolDock {


    /*
     * Scratch pad : copy-pasting zone for useful expression
     *
     * unsigned int atom1Idx =  conformer.bondEnds1Index[rotBdIdx];
            unsigned int atom2Idx =  conformer.bondEnds2Index[rotBdIdx];
     *             Eigen::Vector3d rotAxis = {conformer.x[atom2Idx] - conformer.x[atom1Idx],
                             conformer.y[atom2Idx] - conformer.y[atom1Idx],
                             conformer.z[atom2Idx] - conformer.z[atom1Idx]};
     *
     */

    inline Eigen::Quaternion<double> QuaternionIdentityInit() {
        Eigen::Quaternion<double> qt(1.0,
                                     0.0,
                                     0.0,
                                     0.0);
        return qt;
    }

    //! Holds a gradient for the score function with regard to an iTransform


    //! Returns a neutral transform
    inline iTransform iTransformIdentityInit() {
        iTransform tr;
        tr.transl.x() = 0.0;
        tr.transl.y() = 0.0;
        tr.transl.z() = 0.0;
        tr.rota = QuaternionIdentityInit();
        for (double &bondAngle : tr.bondRotationsAngles) {
            bondAngle = 0.0;
        }
        return tr;
    }

    inline iTransform iTransformIdentityInit(unsigned int num_rot_bond) {
        iTransform tr;
        tr.transl.x() = 0.0;
        tr.transl.y() = 0.0;
        tr.transl.z() = 0.0;
        tr.rota = QuaternionIdentityInit();
        for (unsigned int i = 0; i < num_rot_bond; i++) {
            tr.bondRotationsAngles.push_back(0.0);
        }
        return tr;
    }

    //! Apply the rotation given by a quaternion to a three positional coordinate
    inline void applyRotationInPlace(Eigen::Vector3d &vec, const Eigen::Quaternion<double> &qt) {
        vec = qt._transformVector(vec);
    }


    //! Apply the rotation given by a quaternion to a vector
    inline Eigen::Vector3d applyRotation(const Eigen::Vector3d &vec, const Eigen::Quaternion<double> &qt) {
        Eigen::Vector3d res(vec);
        res = qt._transformVector(vec);
        return res;
    }

    //! Apply the provided translation to a vector
    inline void applyTranslationInPlace(Eigen::Vector3d &v, const Eigen::Vector3d &t) {
        v.x() += t.x();
        v.y() += t.y();
        v.z() += t.z();
    }

    //! Apply the provided translation to a vector
    inline Eigen::Vector3d applyTranslation(const Eigen::Vector3d &v, const Eigen::Translation<double, 3> &t) {
        Eigen::Vector3d res(v);
        res = t * v;
        return res;
    }


    //! Normalize a quaternion (which then has norm == 1.0)
    inline void normalizeQuaternionInPlace(Eigen::Quaternion<double> &quat) {
        quat.normalize();
    }

    inline void applyBondRotationInPlace(SmolDock::iConformer &conformer, const iTransform &tr) {
        for (unsigned int rotBondGroupIdx = 0; rotBondGroupIdx < tr.bondRotationsAngles.size(); rotBondGroupIdx++) {
            const double &rotAngleForBond = tr.bondRotationsAngles[rotBondGroupIdx];
            unsigned int atom1Idx = conformer.bondEnds1Index[rotBondGroupIdx];
            unsigned int atom2Idx = conformer.bondEnds2Index[rotBondGroupIdx];

//            std::cout << "-- Rotation --\n";
//            std::cout << "  Atom 1 : " << atom1Idx << "   " << conformer.x[atom1Idx] << "," << conformer.y[atom1Idx] << "," << conformer.z[atom1Idx] << "\n";
//            std::cout << "  Atom 2 : " << atom2Idx << "   "<< conformer.x[atom2Idx] << "," << conformer.y[atom2Idx] << "," << conformer.z[atom2Idx] << "\n";

            Eigen::Vector3d axis = {
                    (conformer.x[atom2Idx] - conformer.x[atom1Idx]),
                    (conformer.y[atom2Idx] - conformer.y[atom1Idx]),
                    (conformer.z[atom2Idx] - conformer.z[atom1Idx])
            };

            Eigen::Translation<double, 3> moveFromOrigin(conformer.x[atom1Idx], conformer.y[atom1Idx],
                                                         conformer.z[atom1Idx]);
            Eigen::Translation<double, 3> moveToOrigin = moveFromOrigin.inverse();

            // std::cout << "  Axis : " << axis.x()  << "," << axis.y() << "," << axis.z() << " norm = " << axis.norm() << "\n";

            axis.normalize();

//            std::cout << "  Axis : " << axis.x()  << "," << axis.y() << "," << axis.z() << " norm = " << axis.norm() << "\n";


            assert((axis.norm() - 1) < 0.01);

            Eigen::AngleAxis<double> rotation(rotAngleForBond, axis);

            Eigen::Matrix<double, 3, 3> rotMatrix = rotation.toRotationMatrix();

//            std::cout << "RotMatrix \n" << rotMatrix << "\n";

            for (unsigned int &rotatedAtomIdx: conformer.rotatableGroups[rotBondGroupIdx]) {
//                std::cout << " Affecting : " << rotatedAtomIdx << "\n";
                Eigen::Vector3d startPosition(conformer.x[rotatedAtomIdx],
                                              conformer.y[rotatedAtomIdx],
                                              conformer.z[rotatedAtomIdx]);

                Eigen::Vector3d endPosition = moveFromOrigin * (rotMatrix * (moveToOrigin * startPosition));

                conformer.x[rotatedAtomIdx] = endPosition.x();
                conformer.y[rotatedAtomIdx] = endPosition.y();
                conformer.z[rotatedAtomIdx] = endPosition.z();
            }

//            std::cout << "-- END Rotation --" << "\n";

        }
    }


    //! Apply the provided iTransform to a vector (only rigid component : translation + global rotation)
    inline void applyRigidTransformInPlace(Eigen::Vector3d &v, const iTransform &tr) {
        v = tr.rota * v;
        v = tr.transl * v;
    }

    //! Apply the provided iTransform to a vector( only rigid component : translation + global rotation)
    inline Eigen::Vector3d applyRigidTransform(const Eigen::Vector3d &v, const iTransform &tr) {
        Eigen::Vector3d res(v);
        applyRigidTransformInPlace(res, tr);
        return res;
    }


    //! Apply the provided iTransform to a confomer, atom-by-atom (only rigid component : translation + global rotation)
    // TODO: rename this function that operate on whole conformer (versus the other overload that operate on a vector, ie a temporary object)
    inline void applyRigidTransformInPlace(SmolDock::iConformer &conformer, const iTransform &tr) {
        for (unsigned int i = 0; i < conformer.x.size(); i++) {
            Eigen::Vector3d posVect = {conformer.x[i], conformer.y[i], conformer.z[i]};
            applyRigidTransformInPlace(posVect, tr);
            conformer.x[i] = posVect.x();
            conformer.y[i] = posVect.y();
            conformer.z[i] = posVect.z();
        }
    }


    inline void applyTransformInPlace(SmolDock::iConformer &conformer, const iTransform &tr) {
        std::terminate(); // Check for deprecated use
        applyRigidTransformInPlace(conformer, tr);
        applyBondRotationInPlace(conformer, tr);

    }


    //! Normalize a quaternion (which then has norm == 1.0)
    inline Eigen::Quaternion<double> normalizeQuaternion(const Eigen::Quaternion<double> &quat) {
        Eigen::Quaternion<double> ret = quat;
        ret.normalize();
        return ret;
    }


}

#endif //SMOLDOCK_INTERNALSUTILITYFUNCTIONS_H
