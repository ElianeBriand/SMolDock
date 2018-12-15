//
// Created by eliane on 14/12/18.
//

#ifndef SMOLDOCK_ITRANSFORM_H
#define SMOLDOCK_ITRANSFORM_H

#include <vector>
#include <memory>
#include <array>

namespace SmolDock {

    struct iQuaternion {
        // Quaternion :
        double s;
        std::array<double,3> v;
    };

    struct iVect {
        double x,y,z;
    };

    struct iTranslation {
        double x,y,z;
    };

    struct iTransform {
        iTranslation transl;

        iQuaternion rota;
    };


    inline iQuaternion iQuaternionZeroInit()
    {
        iQuaternion qt;
        qt.s = 0.0;
        qt.v[0] = 0.0;
        qt.v[1] = 0.0;
        qt.v[2] = 0.0;
        return qt;
    }

    struct iGradient {
        double dx,dy,dz;

        double dw;
        std::array<double,3> dr;
    };
    inline iTransform iTransformZeroInit()
    {
        iTransform tr;
        tr.transl.x = 0.0;
        tr.transl.y = 0.0;
        tr.transl.z = 0.0;
        tr.rota = iQuaternionZeroInit();
        return tr;
    }

    inline std::array<double,3> scale3DArray(std::array<double,3> a, double factor)
    {
        return std::array<double, 3>{a[0]*factor,a[1]*factor,a[2]*factor};
    }

    inline std::array<double,3> crossProduct3DArray(std::array<double,3> a, std::array<double,3> b)
    {
        std::array<double,3> ret;
        ret[0] = (a[1]*b[2]) - (a[2]*b[1]);
        ret[1] = (a[2]*b[0]) - (a[0]*b[2]);
        ret[2] = (a[0]*b[1]) - (a[1]*b[0]);
        return ret;
    }

    inline std::array<double,3> elementWiseAdd(std::array<double,3> a, std::array<double,3> b)
    {
        std::array<double,3> ret;
        ret[0] = a[0] + b[0];
        ret[1] = a[1] + b[1];
        ret[2] = a[2] + b[2];
        return ret;
    }

    inline iVect rotatePosition(std::array<double,3> vec, iQuaternion qt)
    {
        // Formula for quaternion rotation :
        // Quatertion = w (scalar) + |R>
        // |V_rotated> = |V_initial> + 2*|R> x (|R> x |V_initial> + w|V_initial>)
        iVect res;
        std::array<double,3> twoScaledR = scale3DArray(qt.v,2.0);
        std::array<double,3> RxV = crossProduct3DArray(qt.v,vec);
        std::array<double,3> wScaledV = scale3DArray(vec,qt.s);
        auto a = elementWiseAdd(vec, crossProduct3DArray(twoScaledR,elementWiseAdd(RxV,wScaledV)));
        res.x = a[0];
        res.y = a[1];
        res.z = a[2];
        return res;
    }











}

#endif //SMOLDOCK_ITRANSFORM_H
