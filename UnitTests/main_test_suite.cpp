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
 * SmolDock is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SmolDock.  If not, see <https://www.gnu.org/licenses/>.
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
#include "Engines/Internals/iTransform.h"

#define BOOST_TEST_MODULE main_test_module
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <Structures/Protein.h>


BOOST_AUTO_TEST_SUITE(main_test_suite)

    BOOST_AUTO_TEST_CASE(molecule_class_tests) {


        SmolDock::Molecule mol2;
        mol2.populateFromSMILES("COCO");

        /*
        BOOST_CHECK(mol1.numberOfAtoms() == mol2.numberOfAtoms());
        BOOST_CHECK(mol1.numberOfBonds() == mol2.numberOfBonds());

    */


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
        bool res1 = mol1.generateConformer(conformer1, true, 234);

        BOOST_CHECK(res1);


        auto start_conformersgen = std::chrono::system_clock::now();

        std::vector<SmolDock::iConformer> vecConformer;
        unsigned int desired_num_conformer = 20;
        unsigned int res2 = mol1.generateConformers(vecConformer, desired_num_conformer, true, 234);

        auto end_conformersgen = std::chrono::system_clock::now();


        BOOST_CHECK(res2 == desired_num_conformer);

        auto duration_second = static_cast< std::chrono::duration<double> >(end_conformersgen -
                                                                            start_conformersgen).count();
        std::cout << "[ ] Conformer generation took "
                  << duration_second
                  << "s for "
                  << res2 << " of " << desired_num_conformer << " generated." << std::endl;

        BOOST_CHECK(duration_second <
                    1); // we hope to be faster (for ibuprofen, consistently < 2ms) but 50ms/conformer is our "something is broken" alert level



        // We check that the conformers are (at least superficially) differents
        bool coordinate_test_tripped = false;
        for (unsigned int i = 0; i < vecConformer.size() - 1; i++) {
            SmolDock::iConformer &conformer1 = vecConformer.at(i);
            SmolDock::iConformer &conformer2 = vecConformer.at(i + 1);

            for (int j = 0; j < conformer1.x.size(); j++) {
                if (conformer1.x[j] == conformer2.x[j])
                    coordinate_test_tripped = true;
                if (conformer1.y[j] == conformer2.y[j])
                    coordinate_test_tripped = true;
                if (conformer1.z[j] == conformer2.z[j])
                    coordinate_test_tripped = true;
            }
        }

        BOOST_CHECK(coordinate_test_tripped == false);

        std::vector<SmolDock::iConformer> vecConformer2;
        mol1.generateConformers(vecConformer2, desired_num_conformer,
                                ,true,234); // Same seed, same everything --> are they the same conformers ?

        // We check that the conformers are the same
        bool sameconformer_test_tripped = false;
        for (int i = 0; i < vecConformer.size(); i++) {
            SmolDock::iConformer &conformer1 = vecConformer.at(i);
            SmolDock::iConformer &conformer2 = vecConformer2.at(i);

            for (int j = 0; j < conformer1.x.size(); j++) {
                if (conformer1.x[j] != conformer2.x[j])
                    sameconformer_test_tripped = true;
                if (conformer1.y[j] != conformer2.y[j])
                    sameconformer_test_tripped = true;
                if (conformer1.z[j] != conformer2.z[j])
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
        bool res = prot.populateFromPDB("../DockingTests/COX2_Ibuprofen/4PH9_COX2_without_Ibuprofen.pdb"); // COX-2

        BOOST_CHECK(res);

        SmolDock::iProtein iprot = prot.getiProtein();

        BOOST_CHECK(iprot.x.size() == iprot.y.size());
        BOOST_CHECK(iprot.x.size() == iprot.z.size());

        BOOST_CHECK(iprot.x.size() == iprot.type.size());
        BOOST_CHECK(iprot.type.size() == iprot.variant.size());
        BOOST_CHECK(iprot.x.size() == iprot.atomicRadius.size());

        std::vector<SmolDock::iConformer> vecConformer;
        unsigned int desired_num_conformer = 3;
        mol1.generateConformers(vecConformer, desired_num_conformer, true, 234);

        for (auto &conformer: vecConformer) {
            BOOST_CHECK(conformer.x.size() == conformer.y.size());
            BOOST_CHECK(conformer.x.size() == conformer.z.size());

            BOOST_CHECK(conformer.x.size() == conformer.type.size());
            BOOST_CHECK(conformer.type.size() == conformer.variant.size());
            BOOST_CHECK(conformer.x.size() == conformer.atomicRadius.size());
        }


    }


BOOST_AUTO_TEST_SUITE_END();