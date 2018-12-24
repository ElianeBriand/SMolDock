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


    inline iQuaternion iQuaternionIdentityInit()
    {
        iQuaternion qt;
        qt.s = 1.0;
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

    inline iTransform iTransformIdentityInit()
    {
        iTransform tr;
        tr.transl.x = 0.0;
        tr.transl.y = 0.0;
        tr.transl.z = 0.0;
        tr.rota = iQuaternionIdentityInit();
        return tr;
    }

    inline void applyRotationInPlace(iVect& vec, const iQuaternion& qt)
    {
        iVect v = vec;
        vec.x = qt.s*qt.s*v.x + 2*qt.v[1]*qt.s*v.z - 2*qt.v[2]*qt.s*v.y + qt.v[0]*qt.v[0]*v.x + 2*qt.v[1]*qt.v[0]*v.y + 2*qt.v[2]*qt.v[0]*v.z - qt.v[2]*qt.v[2]*v.x - qt.v[1]*qt.v[1]*v.x;
        vec.y = 2*qt.v[0]*qt.v[1]*v.x + qt.v[1]*qt.v[1]*v.y + 2*qt.v[2]*qt.v[1]*v.z + 2*qt.s*qt.v[2]*v.x - qt.v[2]*qt.v[2]*v.y + qt.s*qt.s*v.y - 2*qt.v[0]*qt.s*v.z - qt.v[0]*qt.v[0]*v.y;
        vec.z = 2*qt.v[0]*qt.v[2]*v.x + 2*qt.v[1]*qt.v[2]*v.y + qt.v[2]*qt.v[2]*v.z - 2*qt.s*qt.v[1]*v.x - qt.v[1]*qt.v[1]*v.z + 2*qt.s*qt.v[0]*v.y - qt.v[0]*qt.v[0]*v.z + qt.s*qt.s*v.z;
    }

    inline iVect applyRotation(const iVect& v,const iQuaternion& qt)
    {
        iVect res(v);
        applyRotationInPlace(res,qt);
        return res;
    }

    inline void applyTranslationInPlace(iVect& v,const iTranslation& t)
    {
        v.x += t.x;
        v.y += t.y;
        v.z += t.z;
    }

    inline iVect applyTranslation(const iVect& v,const iTranslation& t)
    {
        iVect res(v);
        applyTranslationInPlace(res,t);
        return res;
    }



    inline void applyTransformInPlace(iVect& v, const iTransform& tr)
    {
        applyRotationInPlace(v,tr.rota);
        applyTranslationInPlace(v,tr.transl);
    }

    inline iVect applyTransform(const iVect& v, const iTransform& tr)
    {
        iVect res(v);
        applyTransformInPlace(res,tr);
        return res;
    }


    inline void applyTransformInPlace(iConformer& conformer, const iTransform& tr)
    {
        for (unsigned int i = 0; i < conformer.x.size(); i++) {
            iVect posVect = {conformer.x[i],conformer.y[i],conformer.z[i]};
            applyTransformInPlace(posVect, tr);
            conformer.x[i] = posVect.x;
            conformer.y[i] = posVect.y;
            conformer.z[i] = posVect.z;
        }
    }










}

#endif //SMOLDOCK_ITRANSFORM_H
