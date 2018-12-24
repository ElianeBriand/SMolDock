

#include "Engines/Internals/iTransform.h"


#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_SUITE(engine_test_suite)


    BOOST_AUTO_TEST_CASE(internal_math_iTransform_test) {


        SmolDock::iVect vect1{1.0,2.0,3.0};

        SmolDock::iTransform neutralTransform = SmolDock::iTransformIdentityInit();
        SmolDock::iQuaternion neutralQt = SmolDock::iQuaternionIdentityInit();

        // Check quaternion identity
        SmolDock::iVect resVect1 = SmolDock::applyRotation(vect1, neutralQt);

        BOOST_CHECK((resVect1.x - (+1.0)) < 0.001);
        BOOST_CHECK((resVect1.y - (+2.0)) < 0.001);
        BOOST_CHECK((resVect1.z - (+3.0)) < 0.001);

        // Check neutral transform

        SmolDock::iVect resVect2 = SmolDock::applyTransform(vect1,neutralTransform);

        BOOST_CHECK((resVect2.x - (+1.0)) < 0.001);
        BOOST_CHECK((resVect2.y - (+2.0)) < 0.001);
        BOOST_CHECK((resVect2.z - (+3.0)) < 0.001);

        SmolDock::iQuaternion rotateXQt = SmolDock::iQuaternionIdentityInit();
        rotateXQt.s = 0;
        rotateXQt.v[0] = 1;

        SmolDock::iVect resVect3 = SmolDock::applyRotation(vect1,rotateXQt);

        BOOST_CHECK((resVect3.x - (+1.0)) < 0.001);
        BOOST_CHECK((resVect3.y - (-2.0)) < 0.001);
        BOOST_CHECK((resVect3.z - (-3.0)) < 0.001);

    }

BOOST_AUTO_TEST_SUITE_END();