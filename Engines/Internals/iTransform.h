//
// Created by eliane on 14/12/18.
//

#ifndef SMOLDOCK_ITRANSFORM_H
#define SMOLDOCK_ITRANSFORM_H

#include <vector>
#include <memory>
#include <array>
#include <cmath>
#include <cassert>


namespace SmolDock {

    //! Represents a rotation
    /*!
     *
     *
     * \sa iTransform
     */
    struct iQuaternion {
        // Quaternion :
        double s; //! Scalar part
        double u, v, t; //! Vector part
    };

    //! General purpose 3D-vector type
    struct iVect {
        double x, y, z;
    };

    //! Represents a translation
    struct iTranslation {
        double x, y, z;
    };


    //! Represents a rotation and translation transform to be applied on ligand/ligand part
    /*!
     * \sa iTranslation, iQuaternion
    */
    struct iTransform {
        iTranslation transl;
        iQuaternion rota;
    };

    /*!
     * \sa iTranslation, iQuaternion
     */
    struct iGradient {
        double dx, dy, dz; //! Translation part

        double ds; //! Scalar quaternion part
        double du, dv, dt; //! Vector quaternion part
    };


}

#endif //SMOLDOCK_ITRANSFORM_H
