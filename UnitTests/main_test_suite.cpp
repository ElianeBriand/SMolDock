//
// Created by eliane on 13/11/18.
//

#include <chrono>

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

#define BOOST_TEST_MODULE main_test_module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_SUITE(main_test_suite)

    BOOST_AUTO_TEST_CASE(molecule_class_tests) {

        SmolDock::UnitTestHelper helper;


        /****
         *
         * Summary of test :
         *
         * Test 1 : Compare manually populated mol and SMILES populated mol
         * Test 2 : Handling of implicit hydrogen : in SMILES and added afterward
         * Test 3 : Checking for error with an invalid SMILES
         * Test 4 : Big but valid SMILES
         */

        SmolDock::Molecule mol1;
        /* Molecule created : CH3-O-CH2-OH */
        helper.populateMol_COCOH(&mol1);

        SmolDock::Molecule mol2;
        mol2.populateFromSMILES("COCO");

        BOOST_CHECK(mol1.numberOfAtoms() == mol2.numberOfAtoms());
        BOOST_CHECK(mol1.numberOfBonds() == mol2.numberOfBonds());


        SmolDock::Molecule mol3;
        mol3.populateFromSMILES("[CH3]O[CH2][OH]");

        auto fingerprint_mol2 = RDKit::RDKFingerprintMol(*(mol2.getInternalRWMol()));
        auto fingerprint_mol3 = RDKit::RDKFingerprintMol(*(mol3.getInternalRWMol()));

        BOOST_CHECK(*fingerprint_mol2 == *fingerprint_mol3);


        SmolDock::Molecule failed_mol;
        bool ret_failed_smiles = failed_mol.populateFromSMILES("HcndDH"); // Garbage SMILES

        BOOST_CHECK(!ret_failed_smiles);


        SmolDock::Molecule mol4;
        // This is ciclosporin
        auto start_largesmile = std::chrono::system_clock::now();
        mol4.populateFromSMILES("CC[C@H]1C(=O)N(CC(=O)N([C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@@H]"
                                "(C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N1)[C@@H]([C@H](C)C/C"
                                "=C/C)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C");
        auto end_largesmile = std::chrono::system_clock::now();
        std::cout << "[ ] Ciclosporin SMILE parsing took "
                  << static_cast< std::chrono::duration<double> >(end_largesmile - start_largesmile).count()
                  << std::endl;
        BOOST_CHECK(mol4.numberOfAtoms() == 196);

    }

    BOOST_AUTO_TEST_CASE(conformer_generation_tests) {

        SmolDock::Molecule mol1;
        // This is ibuprofen
        mol1.populateFromSMILES("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O");


        SmolDock::iConformer conformer1;
        bool res1 = mol1.generateConformer(conformer1, 234);

        BOOST_CHECK(res1);


        auto start_conformersgen = std::chrono::system_clock::now();

        std::vector<SmolDock::iConformer> vec;
        unsigned int desired_num_conformer = 20;
        unsigned int res2 = mol1.generateConformers(&vec, desired_num_conformer, 234);

        auto end_conformersgen = std::chrono::system_clock::now();


        BOOST_CHECK(res2 == desired_num_conformer);

        auto duration_second = static_cast< std::chrono::duration<double> >(end_conformersgen - start_conformersgen).count();
        std::cout << "[ ] Conformer generation took "
                  << duration_second
                  << "s for "
                  << res2 << " of " << desired_num_conformer << " generated." << std::endl;

        BOOST_CHECK(duration_second < 1); // we hope to be faster,  but 50ms/conformer is our alert level.

    }

BOOST_AUTO_TEST_SUITE_END();