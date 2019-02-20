//
// Created by eliane on 23/12/18.
//

#include "PDBLigandUtils.h"
#include <boost/log/trivial.hpp>

namespace SmolDock {


    void AssignBondOrderFromTemplateSMILES(std::shared_ptr<RDKit::RWMol> mol, const std::string &smiles) {
        RDKit::RWMol* attempt_ptr = nullptr;
        try {
            /* This can throw */
            attempt_ptr = RDKit::SmilesToMol(smiles);
        }
        catch (...) {
            attempt_ptr = nullptr;
        }
        if (attempt_ptr == nullptr) {
            BOOST_LOG_TRIVIAL(error)
                << "Error in AssignBondOrderFromTemplateSMILES(std::shared_ptr<RDKit::RWMol> mol,const std::string& smiles)";
            BOOST_LOG_TRIVIAL(error) << "for teplate smiles = " << smiles;
            BOOST_LOG_TRIVIAL(error) << "The SMILES string was not correctly parsed by RDKit.";
            BOOST_LOG_TRIVIAL(error)
                << "This often indicated a malformed SMILES. (but not always, RDKit has parsing bugs)";
            BOOST_LOG_TRIVIAL(error) << "Non-modified mol will be used (likely with incorrect bond order)";
        }

        std::shared_ptr<RDKit::RWMol> templatemol(attempt_ptr);


        RDKit::MolOps::removeHs(*(mol));

        // The following code, until the END RDKIT CODE comment, is a C++ translation
        // of python code from RDKit, in particular, the function AssignBondOrdersFromTemplate
        //
        // As a reminder, the copyright for this code is
        //  Copyright (C) 2001-2017 Greg Landrum and Rational Discovery LLC
        //
        //   @@ All Rights Reserved @@
        //  The contents are covered by the terms of the BSD license
        //  which is included in the file COPYING
        //

        // RDKit::MatchVectType matches = std::vector<std::pair<int, int>>
        //  The format is (queryAtomIdx, molAtomIdx)
        RDKit::MatchVectType matches;
        bool matching = RDKit::SubstructMatch((RDKit::ROMol) *mol, (RDKit::ROMol) *templatemol, matches);

        if (!matching) {
            // We make a copy on which we will remove all double bond/charges/... to compare
            std::shared_ptr<RDKit::RWMol> templatemol_plainstruct(new RDKit::RWMol(*templatemol));
            // Set all bonds to single, non-aromatic
            for (unsigned int i = 0; i < mol->getNumBonds(); i++) {
                auto bond = mol->getBondWithIdx(i);
                bond->setIsAromatic(false);
                bond->setBondType(RDKit::Bond::BondType::SINGLE);
            }

            // Also in templatemol_plainstruct
            for (unsigned int i = 0; i < templatemol_plainstruct->getNumBonds(); i++) {
                auto bond = templatemol_plainstruct->getBondWithIdx(i);
                bond->setIsAromatic(false);
                bond->setBondType(RDKit::Bond::BondType::SINGLE);
            }

            // Remove all charges
            for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
                auto atom = mol->getAtomWithIdx(i);
                atom->setFormalCharge(0);
            }

            for (unsigned int i = 0; i < templatemol_plainstruct->getNumAtoms(); i++) {
                auto atom = templatemol_plainstruct->getAtomWithIdx(i);
                atom->setFormalCharge(0);
            }


            bool matching = RDKit::SubstructMatch((RDKit::ROMol) *mol, (RDKit::ROMol) *templatemol_plainstruct,
                                                  matches);
            if (!matching) // Cant reconcile even without charges or double bond
            {
                BOOST_LOG_TRIVIAL(error)
                    << "Cannot match molecule with SMILES template, thus cannot adjust bond order ";
                BOOST_LOG_TRIVIAL(error) << "Plain mol will be used (all bond order single, no charges, no nothing)";
                return;
            }

        }

        // After this point, we are matching (or we returned)
        // And we have matches in RDKit::MatchVectType matches

        // For each bond of the template, we find the corresponding one in the molecule to edit
        // Then we set it to the same characteristics
        for (int i = 0; i < templatemol->getNumBonds(); i++) {
            auto template_bond = templatemol->getBondWithIdx(i);

            // RDKit::MatchVectType matches = std::vector<std::pair<int, int>>
            //  The format is (queryAtomIdx, molAtomIdx)
            // So we find_if for matching queryAtomIdx to the beginning and end of the bond
            int idx_atom1 = std::get<1>(*(std::find_if(matches.begin(), matches.end(),
                                                       [&](const std::pair<int, int> &e) {
                                                           return (std::get<0>(e) == template_bond->getBeginAtomIdx());
                                                       })));
            int idx_atom2 = std::get<1>(*(std::find_if(matches.begin(), matches.end(),
                                                       [&](const std::pair<int, int> &e) {
                                                           return (std::get<0>(e) == template_bond->getEndAtomIdx());
                                                       })));

            auto bond_to_edit = mol->getBondBetweenAtoms(idx_atom1, idx_atom2);
            bond_to_edit->setBondType(template_bond->getBondType());
            bond_to_edit->setIsAromatic(template_bond->getIsAromatic());
        }

        for (unsigned int i = 0; i < templatemol->getNumAtoms(); i++) {
            auto template_atom = templatemol->getAtomWithIdx(i);

            int idx_atom_to_edit = std::get<1>(*(std::find_if(matches.begin(), matches.end(),
                                                              [&](const std::pair<int, int> &e) {
                                                                  return (std::get<0>(e) == i);
                                                              })));

            auto atom_to_edit = mol->getAtomWithIdx(idx_atom_to_edit);

            atom_to_edit->setHybridization(template_atom->getHybridization());
            atom_to_edit->setIsAromatic(template_atom->getIsAromatic());
            atom_to_edit->setNumExplicitHs(template_atom->getNumExplicitHs());
            atom_to_edit->setFormalCharge(template_atom->getFormalCharge());
            atom_to_edit->setChiralTag(template_atom->getChiralTag());
        }

        RDKit::MolOps::sanitizeMol(*mol);

        // //////////////// END RDKIT CODE  //////////////////////////////////


    }

}