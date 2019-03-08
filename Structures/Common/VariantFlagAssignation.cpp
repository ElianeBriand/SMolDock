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

#include <cassert>


#include "VariantFlagAssignation.h"
#include "Structures/Bond.h"

#include <boost/log/trivial.hpp>

#include "PDBResiduePropertiesTable.h"

namespace SmolDock {


    void assignVariantFlagsForResidueAtom(AminoAcid &residue, PDBResidueVariantAssignationType assignation_type) {

        std::set<std::tuple<AminoAcid::AAType, std::vector<std::tuple<std::string, Atom::AtomVariant> > > > *ResidueAtomPropertiesLookupTable;
        if (assignation_type == PDBResidueVariantAssignationType::GeneralPurpose) {
            ResidueAtomPropertiesLookupTable = &ResidueAtomPropertiesLookupTable_General;
        } else {
            ResidueAtomPropertiesLookupTable = &ResidueAtomPropertiesLookupTable_General;
        }

        AminoAcid::AAType type = residue.getType();
        if (type == AminoAcid::AAType::unknown) {
            BOOST_LOG_TRIVIAL(error)
                << "Encountered unknown residue type while assigning atom variant. Check PDB file. Undefined behaviour may occur.";

        } else {
            auto record_it = std::find_if(std::begin(*ResidueAtomPropertiesLookupTable),
                                          std::end(*ResidueAtomPropertiesLookupTable),
                                          [&](const std::tuple<AminoAcid::AAType, std::vector<std::tuple<std::string, Atom::AtomVariant> > > &e) {
                                              return std::get<0>(e) == type;
                                          });
            if (record_it == std::end(*ResidueAtomPropertiesLookupTable)) {
                // Not found
                BOOST_LOG_TRIVIAL(error)
                    << "Amino acid not found in ResidueAtomPropertiesLookupTable. Variant assignation not done for this amino acid.";
                return;
            }


            std::vector<std::tuple<std::string, Atom::AtomVariant> > record_for_this_AA = std::get<1>(*record_it);


            for (auto &atom : residue.atoms) {
                std::string complete_atom_name = atom->rawPDBAtomName;
                auto nameTypeTuple_it = std::find_if(std::begin(record_for_this_AA), std::end(record_for_this_AA),
                                                     [&](const std::tuple<std::string, Atom::AtomVariant> &e) {
                                                         return std::get<0>(e) == complete_atom_name;
                                                     });
                if (nameTypeTuple_it != std::end(record_for_this_AA)) {
                    // We have an association between complete atom name and AtomVariant
                    Atom::AtomVariant newVariant = atom->getAtomVariant() | std::get<1>(*nameTypeTuple_it);
                    atom->setAtomVariant(newVariant);
                } else {
                    BOOST_LOG_TRIVIAL(error)
                        << "Missing AtomVariant record for residue " << resTypeToString(type) << ", atom "
                        << complete_atom_name;
                }

            }

        }

    }

    void assignApolarCarbonFlag(std::vector<std::shared_ptr<Atom> > &atomVect) {
        for (auto &atom: atomVect) {
            if (atom->getAtomType() == Atom::AtomType::carbon) {
                bool isLinkedToHeteroatom = false;
                for (const auto &bond_weakptr: atom->bonds) {
                    std::shared_ptr<Bond> bond = bond_weakptr.lock();
                    if (bond == nullptr) {
                        BOOST_LOG_TRIVIAL(error)
                            << "Encountered nullptr weak_ptr while traversing molecule. This is definitely a bug: please report.";
                        BOOST_LOG_TRIVIAL(error) << "Exiting unsucessfully...";
                        std::exit(6);
                    }
                    std::shared_ptr<Atom> endA = bond->getEndA();
                    std::shared_ptr<Atom> endB = bond->getEndB();

                    assert((endA.get() == atom.get()) !=
                           (endB.get() ==
                            atom.get())); // Check that the current atom is exactly one of the end (!= is logical XOR for bool)

                    std::shared_ptr<Atom> opposite_end;
                    if (endA.get() == atom.get())
                        // The current atom is EndA
                        opposite_end = endB;
                    else if (endB.get() == atom.get())
                        // The current atom is EndB
                        opposite_end = endA;

                    Atom::AtomType linked_atom_type = opposite_end->getAtomType();
                    if ((linked_atom_type != Atom::AtomType::carbon) &&
                        (linked_atom_type != Atom::AtomType::hydrogen)) {
                        isLinkedToHeteroatom = true;
                    }
                }
                if (!isLinkedToHeteroatom) {
                    atom->variant = atom->variant | Atom::AtomVariant::apolar;
                }
            }
        }
    }

}