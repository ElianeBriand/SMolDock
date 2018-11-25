//
// Created by eliane on 13/11/18.
//

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE main_test_suite

#include <boost/test/unit_test.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
/*
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionRunner.h>
#include <GraphMol/ChemReactions/PreprocessRxn.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>
 */
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <DataStructs/ExplicitBitVect.h>

#define SMOLDOCK_VERBOSE_DEBUG

#include "../Structures/Molecule.h"
#include "UnitTestHelper.h"

BOOST_AUTO_TEST_CASE(molecule_from_smiles_tests) {

    SmolDock::UnitTestHelper helper;


    /****
     *
     * Summary of test :
     *
     * Test 1 : Compare manually populated mol and SMILES populated mol
     * Test 2 : Handling of implicit hydrogen : in smiles and added afterward
     *
     */

    std::shared_ptr<SmolDock::Molecule> mol1 = std::make_shared<SmolDock::Molecule>();
    /* Molecule created : CH3-O-CH2-OH */
    helper.populateMol_COCOH(mol1);

    std::shared_ptr<SmolDock::Molecule> mol2 = std::make_shared<SmolDock::Molecule>("COCO");

    auto AtomAndBondNumber_mol1 = helper.getNumberOfAtomAndBonds(mol1);
    auto AtomAndBondNumber_mol2 = helper.getNumberOfAtomAndBonds(mol2);
    BOOST_CHECK(AtomAndBondNumber_mol1 == AtomAndBondNumber_mol2);


    std::shared_ptr<SmolDock::Molecule> mol3 = std::make_shared<SmolDock::Molecule>("[CH3]O[CH2][OH]");

    auto fingerprint_mol2 = RDKit::RDKFingerprintMol(*(mol2->getInternalRWMol()));
    auto fingerprint_mol3 = RDKit::RDKFingerprintMol(*(mol3->getInternalRWMol()));

    BOOST_CHECK(*fingerprint_mol2 == *fingerprint_mol3);


}

