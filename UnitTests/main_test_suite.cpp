//
// Created by eliane on 13/11/18.
//

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE main_test_suite

#include <boost/test/unit_test.hpp>

#include "../Structures/Molecule.h"

BOOST_AUTO_TEST_CASE(structure_tests) {
    std::shared_ptr<SmolDock::Molecule> mol = std::make_shared<SmolDock::Molecule>();
    mol->_dev_populateSampleMolecule();


    BOOST_CHECK(1 == 1);
}
