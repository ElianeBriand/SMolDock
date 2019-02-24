

#include "Engines/Internals/iTransform.h"
#include "Engines/Internals/InternalsUtilityFunctions.h"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <iostream>

#include <Utilities/IntermediateConformerCollector.h>

SmolDock::IntermediateConformerCollector* conformerCollector;


BOOST_AUTO_TEST_SUITE(engine_test_suite)


    BOOST_AUTO_TEST_CASE(internal_math_iTransform_test) {


        SmolDock::iVect vect1{1.0, 2.0, 3.0};

        SmolDock::iTransform neutralTransform = SmolDock::iTransformIdentityInit();
        SmolDock::iQuaternion neutralQt = SmolDock::iQuaternionIdentityInit();

        // Check quaternion identity
        SmolDock::iVect resVect1 = SmolDock::applyRotation(vect1, neutralQt);

        BOOST_CHECK((resVect1.x - (+1.0)) < 0.001);
        BOOST_CHECK((resVect1.y - (+2.0)) < 0.001);
        BOOST_CHECK((resVect1.z - (+3.0)) < 0.001);

        // Check neutral transform

        SmolDock::iVect resVect2 = SmolDock::applyTransform(vect1, neutralTransform);

        BOOST_CHECK((resVect2.x - (+1.0)) < 0.001);
        BOOST_CHECK((resVect2.y - (+2.0)) < 0.001);
        BOOST_CHECK((resVect2.z - (+3.0)) < 0.001);

        SmolDock::iQuaternion rotateXQt = SmolDock::iQuaternionIdentityInit();
        rotateXQt.s = std::sqrt(1.0 - std::pow(0.02, 2));
        rotateXQt.u = 0.02;

        SmolDock::iVect resVect3 = SmolDock::applyRotation(vect1, rotateXQt);

        BOOST_CHECK((resVect3.x - (+1.0)) < 0.001);
        BOOST_CHECK((resVect3.y - (1.87842)) < 0.001);
        BOOST_CHECK((resVect3.z - (3.07758)) < 0.001);


    }

    BOOST_AUTO_TEST_CASE(internal_type_checkProperties) {
        BOOST_CHECK(std::is_polymorphic<SmolDock::iConformer>::value == false);

    }


BOOST_AUTO_TEST_SUITE_END();