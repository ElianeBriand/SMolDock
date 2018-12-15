/*
 * Copyright (c) 2018 Eliane Briand
 *
 * This file is part of SmolDock.
 *
 * SmolDock is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Foobar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

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

#ifndef SMOLDOCK_VERBOSE_DEBUG
#define SMOLDOCK_VERBOSE_DEBUG
#endif

#include "../Structures/Molecule.h"
#include "UnitTestHelper.h"
#include "Engines/Internals/iTransform.h"

#define BOOST_TEST_MODULE main_test_module
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <Structures/Protein.h>


BOOST_AUTO_TEST_SUITE(main_test_suite)

    BOOST_AUTO_TEST_CASE(molecule_class_tests) {



/*  We dont care about the hydrogen anyore
        SmolDock::Molecule mol1;
        // Molecule created : CH3-O-CH2-OH
        helper.populateMol_COCOH(&mol1);
*/

        SmolDock::Molecule mol2;
        mol2.populateFromSMILES("COCO");

        /*
        BOOST_CHECK(mol1.numberOfAtoms() == mol2.numberOfAtoms());
        BOOST_CHECK(mol1.numberOfBonds() == mol2.numberOfBonds());

    */

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
        BOOST_CHECK(mol4.numberOfAtoms() == 85);

    }

    BOOST_AUTO_TEST_CASE(conformer_generation_tests) {

        SmolDock::Molecule mol1;
        // This is ibuprofen
        mol1.populateFromSMILES("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O");


        SmolDock::iConformer conformer1;
        bool res1 = mol1.generateConformer(conformer1, 234);

        BOOST_CHECK(res1);


        auto start_conformersgen = std::chrono::system_clock::now();

        std::vector<SmolDock::iConformer> vecConformer;
        unsigned int desired_num_conformer = 20;
        unsigned int res2 = mol1.generateConformers(vecConformer, desired_num_conformer, 234);

        auto end_conformersgen = std::chrono::system_clock::now();


        BOOST_CHECK(res2 == desired_num_conformer);

        auto duration_second = static_cast< std::chrono::duration<double> >(end_conformersgen - start_conformersgen).count();
        std::cout << "[ ] Conformer generation took "
                  << duration_second
                  << "s for "
                  << res2 << " of " << desired_num_conformer << " generated." << std::endl;

        BOOST_CHECK(duration_second < 1); // we hope to be faster (for ibuprofen, consistently < 2ms) but 50ms/conformer is our "something is broken" alert level



        // We check that the conformers are (at least superficially) differents
        bool coordinate_test_tripped = false;
        for(int i = 0; i < vecConformer.size()-1; i++)
        {
            SmolDock::iConformer& conformer1 = vecConformer.at(i);
            SmolDock::iConformer& conformer2 = vecConformer.at(i+1);

            for(int j = 0; j < conformer1.x.size(); j++)
            {
                if(conformer1.x[j] == conformer2.x[j])
                    coordinate_test_tripped = true;
                if(conformer1.y[j] == conformer2.y[j])
                    coordinate_test_tripped = true;
                if(conformer1.z[j] == conformer2.z[j])
                    coordinate_test_tripped = true;
            }
        }

        BOOST_CHECK(coordinate_test_tripped == false);

        std::vector<SmolDock::iConformer> vecConformer2;
        mol1.generateConformers(vecConformer2, desired_num_conformer, 234); // Same seed, same everything --> are they the same conformers ?

        // We check that the conformers are the same
        bool sameconformer_test_tripped = false;
        for(int i = 0; i < vecConformer.size(); i++)
        {
            SmolDock::iConformer& conformer1 = vecConformer.at(i);
            SmolDock::iConformer& conformer2 = vecConformer2.at(i);

            for(int j = 0; j < conformer1.x.size(); j++)
            {
                if(conformer1.x[j] != conformer2.x[j])
                    sameconformer_test_tripped = true;
                if(conformer1.y[j] != conformer2.y[j])
                    sameconformer_test_tripped = true;
                if(conformer1.z[j] != conformer2.z[j])
                    sameconformer_test_tripped = true;
            }
        }

        BOOST_CHECK(sameconformer_test_tripped == false);


    }

    BOOST_AUTO_TEST_CASE(internal_representation_generation_test) {

        SmolDock::Molecule mol1;
        mol1.populateFromSMILES("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"); // Ibuprofen

        SmolDock::Protein prot;
        // prot.populateFromPDB("1dpx.pdb"); // Lysozyme
        bool res = prot.populateFromPDB("../DockingTests/COX2_Ibuprofen/3LN1_NoHeme_NoLigand.pdb"); // COX-2

        BOOST_CHECK(res);

        SmolDock::iProtein iprot = prot.getiProtein();

        BOOST_CHECK(iprot.x.size() == iprot.y.size());
        BOOST_CHECK(iprot.x.size() == iprot.z.size());

        BOOST_CHECK(iprot.x.size() == iprot.type.size());
        BOOST_CHECK(iprot.type.size() == iprot.variant.size());
        BOOST_CHECK(iprot.x.size() == iprot.atomicRadius.size());

        std::vector<SmolDock::iConformer> vecConformer;
        unsigned int desired_num_conformer = 3;
        mol1.generateConformers(vecConformer, desired_num_conformer, 234);

        for(auto& conformer: vecConformer)
        {
            BOOST_CHECK(conformer.x.size() == conformer.y.size());
            BOOST_CHECK(conformer.x.size() == conformer.z.size());

            BOOST_CHECK(conformer.x.size() == conformer.type.size());
            BOOST_CHECK(conformer.type.size() == conformer.variant.size());
            BOOST_CHECK(conformer.x.size() == conformer.atomicRadius.size());
        }



    }

    BOOST_AUTO_TEST_CASE(internal_math_iTransform_test) {
        // Some example vectors
        std::array<double,3> vec1 = {1.0,2.0,0.0};
        std::array<double,3> vec2 = {2.0,1.0,3.0};
        std::array<double,3> vec3 = {120.0,-2.0,3.0};

        // Testing the cross product
        std::array<double,3> vecRes1 = SmolDock::crossProduct3DArray(vec1, vec2);
        std::array<double,3> vecRes2 = SmolDock::crossProduct3DArray(vec3, vec2);

        BOOST_CHECK((vecRes1[0] - (+6.0)) < 0.001);
        BOOST_CHECK((vecRes1[1] - (-3.0)) < 0.001);
        BOOST_CHECK((vecRes1[2] - (-3.0)) < 0.001);

        BOOST_CHECK((vecRes2[0] - (-9.0)) < 0.001);
        BOOST_CHECK((vecRes2[1] - (-354.0)) < 0.001);
        BOOST_CHECK((vecRes2[2] - (+124.0)) < 0.001);

        // Testing the scalar scaling
        std::array<double,3> vecRes3 = SmolDock::scale3DArray(vec3, 3.0);
        std::array<double,3> vecRes4 = SmolDock::scale3DArray(vec3, -1.0);

        BOOST_CHECK((vecRes3[0] - (+360.0)) < 0.001);
        BOOST_CHECK((vecRes3[1] - (-6.0)) < 0.001);
        BOOST_CHECK((vecRes3[2] - (+9.0)) < 0.001);

        BOOST_CHECK((vecRes4[0] - (-120.0)) < 0.001);
        BOOST_CHECK((vecRes4[1] - (+2.0)) < 0.001);
        BOOST_CHECK((vecRes4[2] - (-3.0)) < 0.001);

    }

BOOST_AUTO_TEST_SUITE_END();