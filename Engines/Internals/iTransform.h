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

#include <Eigen/Core>
#include <Eigen/Dense>

namespace SmolDock {

    //! Represents a rotation
    /*!
     *
     *
     * \sa iTransform
     */

    //! Represents a rotation and translation transform to be applied on ligand/ligand part
    /*!
     * \sa iTranslation, Eigen::Quaternion<double>
    */
    struct iTransform {
        // Global component
        Eigen::Translation<double, 3> transl;
        Eigen::Quaternion<double> rota;

        std::vector<double> bondRotationsAngles;
    };


}

#endif //SMOLDOCK_ITRANSFORM_H
