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

#include <tuple>

#include "Molecule.h"
#include "Atom.h"
#include "Utilities/PDBLigandUtils.h"

#include <ESBTL/default.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/weighted_atom_iterator.h>
#include <ESBTL/compressed_ifstream.h>


#include <boost/algorithm/string.hpp>
#include <boost/log/trivial.hpp>

namespace SmolDock {

    Molecule::Molecule() = default;


    bool Molecule::populateFromSMILES(const std::string &smiles, unsigned int seed) {
        RDKit::RWMol *a = nullptr;
        try {
            /* This can throw */
            a = RDKit::SmilesToMol(smiles);
        }
        catch (...) {
            a = nullptr;
        }
        if (a == nullptr) {
            BOOST_LOG_TRIVIAL(error) << "Error in constructor Molecule::Molecule(const std::string &smiles)";
            BOOST_LOG_TRIVIAL(error) << "for smiles = " << smiles;
            BOOST_LOG_TRIVIAL(error) << "The SMILES string was not correctly parsed by RDKit.";
            BOOST_LOG_TRIVIAL(error)
                << "This often indicated a malformed SMILES. (but not always, RDKit has parsing bugs)";
            return false;
        }

        this->rwmol.reset(a);

        if (this->populateInternalAtomAndBondFromRWMol(seed) == false)
            return false;

        this->smiles = smiles;

        return true;
    }


    bool Molecule::generateConformer(iConformer &conformer, int seed) {

        // Hydrogen added for conformer gen
        RDKit::MolOps::addHs(*rwmol);
        int conformer_id = RDKit::DGeomHelpers::EmbedMolecule(*rwmol, 10, seed, false);
        RDKit::MolOps::removeHs(*rwmol);

        if (conformer_id == -1) // Failed to generate
            return false;


        conformer.x.reserve(static_cast<size_t>(this->atoms.size()));
        conformer.y.reserve(static_cast<size_t>(this->atoms.size()));
        conformer.z.reserve(static_cast<size_t>(this->atoms.size()));
        conformer.type.reserve(static_cast<size_t>(this->atoms.size()));

        RDKit::Conformer &rdkit_conformer = rwmol->getConformer(conformer_id);


        for (auto atom_it = rwmol->beginAtoms(); atom_it != rwmol->endAtoms(); ++atom_it) {
            assert((*atom_it)->getAtomicNum() != 1); // Conformer should not contain Hs (which were removed earlier)

            conformer.type.push_back(static_cast<unsigned char>((*atom_it)->getAtomicNum()));

            const RDGeom::Point3D &position = rdkit_conformer.getAtomPos((*atom_it)->getIdx());
            conformer.x.push_back(position.x);
            conformer.y.push_back(position.y);
            conformer.z.push_back(position.z);

            // FIXME : this means there is an encapsulation problem, likely due to the RWMol* being in Molecule but not linked to Atom
            conformer.atomicRadius.push_back(atomTypeToAtomicRadius(
                    stringToAtomType(boost::to_upper_copy<std::string>((*atom_it)->getSymbol()))));
        }

        return true;
    }

    unsigned int Molecule::generateConformers(std::vector<iConformer> &viConformers, unsigned int num, int seed) {


        std::vector<int> conformer_ids = RDKit::DGeomHelpers::EmbedMultipleConfs(*rwmol, // Molecule
                                                                                 num, // Num of conformer
                                                                                 30, // Max attempts
                                                                                 seed, // RNG seed
                                                                                 false); // Erase existing conformer
        if (conformer_ids.size() == 0)
            return 0; // Early failure case

        viConformers.reserve(viConformers.capacity() + conformer_ids.size() + 1);

        // Add the initial conformer
        conformer_ids.push_back(this->initial_conformer_id);

        unsigned int num_atoms = this->atoms.size();

        for (int i : conformer_ids) {
            iConformer conformer;
            RDKit::Conformer rdkit_conformer = this->rwmol->getConformer(i);


            conformer.x.reserve(static_cast<size_t>(num_atoms));
            conformer.y.reserve(static_cast<size_t>(num_atoms));
            conformer.z.reserve(static_cast<size_t>(num_atoms));
            conformer.type.reserve(static_cast<size_t>(num_atoms));

            for (auto atom_it = rwmol->beginAtoms(); atom_it != rwmol->endAtoms(); ++atom_it) {
                conformer.type.push_back(static_cast<unsigned char>((*atom_it)->getAtomicNum()));
                conformer.variant.push_back(0); // Not implemented

                const RDGeom::Point3D &position = rdkit_conformer.getAtomPos((*atom_it)->getIdx());
                conformer.x.push_back(position.x);
                conformer.y.push_back(position.y);
                conformer.z.push_back(position.z);

                // FIXME : this means there is an encapsulation problem, likely due to the RWMol* being in Molecule but not linked to Atom
                conformer.atomicRadius.push_back(atomTypeToAtomicRadius(
                        stringToAtomType(boost::to_upper_copy<std::string>((*atom_it)->getSymbol()))));
            }

            viConformers.push_back(std::move(conformer));
        }

        return conformer_ids.size();
    }

    unsigned int Molecule::numberOfAtoms() {
        return atoms.size();
    }

    unsigned int Molecule::numberOfBonds() {
        return bonds.size();
    }

    bool Molecule::populateFromPDB(const std::string &filename, const std::string &smiles_hint, unsigned int seed) {


        // Flavor = 1 --> ignore alternate location
        RDKit::RWMol *mol = RDKit::PDBFileToMol(filename, false, false, 1, true);

        if (mol == nullptr) {
            return false;
        }

        this->rwmol.reset(mol);

        if (smiles_hint != "") {
            AssignBondOrderFromTemplateSMILES(this->rwmol, smiles_hint);
        }

        populateInternalAtomAndBondFromRWMol(seed);

        return true;
    }

    std::string Molecule::getResidueName() const {
        return this->residue_name;
    }

    void Molecule::setResidueName(const std::string &res_name) {
        this->residue_name = res_name.substr(0, 3);

    }

    bool Molecule::populateInternalAtomAndBondFromRWMol(unsigned int seed) {

        // We care about Hs only during conformer generation
        // TODO: find out if we need to care at other moments...
        RDKit::MolOps::addHs(*(this->rwmol));

        this->initial_conformer_id = 0;
        if (this->rwmol->getNumConformers() == 0) // We need to generate the first conformer
        {
            this->initial_conformer_id = RDKit::DGeomHelpers::EmbedMolecule(*rwmol, 1000, seed, true);
        } else { // Else we will just use the first
            this->initial_conformer_id = (*(this->rwmol->beginConformers()))->getId();
        }

        RDKit::Conformer &starting_conformer = rwmol->getConformer(this->initial_conformer_id);

        // We remove them afterward
        RDKit::MolOps::removeHs((*(this->rwmol)));

        for (auto atom_it = rwmol->beginAtoms(); atom_it != rwmol->endAtoms(); ++atom_it) {
            std::shared_ptr<Atom> current_atom;
            // TODO : refactor this using the std::set of type-name-properties
            switch ((*atom_it)->getAtomicNum()) {
                case 1:
                    current_atom = std::make_shared<Atom>(Atom::AtomType::hydrogen, (*atom_it)->getIdx());
                    break;
                case 6:
                    current_atom = std::make_shared<Atom>(Atom::AtomType::carbon, (*atom_it)->getIdx());
                    break;
                case 7:
                    current_atom = std::make_shared<Atom>(Atom::AtomType::nitrogen, (*atom_it)->getIdx());
                    break;
                case 8:
                    current_atom = std::make_shared<Atom>(Atom::AtomType::oxygen, (*atom_it)->getIdx());
                    break;
                default:
                    std::cout << "[!] Unsupported atom type" << std::endl;
                    std::cout << "[!] Atomic num = " << (*atom_it)->getAtomicNum() << std::endl;
                    return false;
                    // break;
            }
            const RDGeom::Point3D &position = starting_conformer.getAtomPos((*atom_it)->getIdx());
            current_atom->setAtomPosition(std::make_tuple(position.x, position.y, position.z));
            atoms.push_back(current_atom);
        }

        for (auto bond_it = rwmol->beginBonds(); bond_it != rwmol->endBonds(); ++bond_it) {
            auto beginAtom = (*bond_it)->getBeginAtom();
            auto endAtom = (*bond_it)->getEndAtom();

            unsigned int beginAtomID = beginAtom->getIdx();
            unsigned int endAtomID = endAtom->getIdx();

            auto resultBegin = std::find_if(std::begin(atoms), std::end(atoms),
                                            [beginAtomID](const std::shared_ptr<Atom> &atom) -> bool {
                                                return atom->getAtomID() == beginAtomID;
                                            }
            );

            auto resultEnd = std::find_if(std::begin(atoms), std::end(atoms),
                                          [endAtomID](const std::shared_ptr<Atom> &atom) -> bool {
                                              return atom->getAtomID() == endAtomID;
                                          }
            );

            if (resultBegin == std::end(atoms) || resultEnd == std::end(atoms)) {
                std::cout << "[!] Unable to find atom with ID = " << beginAtomID << " or " << endAtomID << std::endl;
                return false;
            }

            auto new_bond = std::make_shared<Bond>(*resultBegin, *resultEnd);


            /*
              For reference, from RDKit's bond.h. We do not intent to support the exotic/dative bond for the
              foreseeable future.

              typedef enum {
                 UNSPECIFIED = 0,
                 SINGLE,
                 DOUBLE,
                 TRIPLE,
                 QUADRUPLE,
                 QUINTUPLE,
                 HEXTUPLE,
                 ONEANDAHALF,
                 TWOANDAHALF,
                 THREEANDAHALF,
                 FOURANDAHALF,
                 FIVEANDAHALF,
                 AROMATIC,
                 IONIC,
                 HYDROGEN,
                 THREECENTER,
                 DATIVEONE,  //!< one-electron dative (e.g. from a C in a Cp ring to a metal)
                 DATIVE,     //!< standard two-electron dative
                 DATIVEL,    //!< standard two-electron dative
                 DATIVER,    //!< standard two-electron dative
                 OTHER,
                 ZERO  //!< Zero-order bond (from
                 // http://pubs.acs.org/doi/abs/10.1021/ci200488k)
              } BondType;

             */

            switch ((*bond_it)->getBondType()) {
                case RDKit::Bond::SINGLE:
                    new_bond->setBondType(Bond::BondType::singlebond);
                    break;
                case RDKit::Bond::DOUBLE:
                    new_bond->setBondType(Bond::BondType::doublebond);
                    break;
                case RDKit::Bond::TRIPLE:
                    new_bond->setBondType(Bond::BondType::triplebond);
                    break;
                case RDKit::Bond::AROMATIC: // As expected aromatic is set for all bond in the cycle
                    // , not just the kekule-style "double bond"
                    new_bond->setBondType(Bond::BondType::aromatic);
                    break;
                default:
                    std::cout << "[!] Unsupported bond type" << std::endl;
                    std::cout << "[!] Enum bond type = " << (*bond_it)->getBondType() << std::endl;
                    return false;
                    // break;
            }
            new_bond->publicizeToAtom();
            bonds.push_back(new_bond);
        }
        return true;
    }

    bool Molecule::updateAtomPositionsFromiConformer(const iConformer &conformer) {

        for (int i = 0; i < atoms.size(); ++i) {
            atoms[i]->setAtomPosition(
                    std::make_tuple(conformer.x[i], conformer.y[i], conformer.z[i])
            );

        }

        return true;
    }

    Molecule Molecule::deepcopy() {
        Molecule newmol(*this);


        std::vector< std::tuple<std::shared_ptr<Atom>,Atom*> > new_old_lookup_table;
        newmol.atoms.clear();
        for(auto& ptr: this->atoms)
        {
            Atom* newatom_ptr = new Atom(*ptr);
            newatom_ptr->bonds.clear(); // We will restore it after rebuilding the bonds
            std::shared_ptr<Atom>& inserted = newmol.atoms.emplace_back(newatom_ptr);
            new_old_lookup_table.push_back(std::make_tuple(inserted,ptr.get())); // Register new/old pair to rebuild bonds

        }

        newmol.bonds.clear();
        for(auto& ptr: this->bonds)
        {
            Atom* endA_ptr = ptr->getEndA().get();
            Atom* endB_ptr = ptr->getEndB().get();

            // Find matches in lookup table :
            std::shared_ptr<Atom> endA = std::get<0>(*(std::find_if(std::begin(new_old_lookup_table), std::end(new_old_lookup_table),
                                   [&](const std::tuple<std::shared_ptr<Atom>,Atom*> &e) {
                                       return (std::get<1>(e) == endA_ptr);
                                   })));

            std::shared_ptr<Atom> endB = std::get<0>(*(std::find_if(std::begin(new_old_lookup_table), std::end(new_old_lookup_table),
                                                                    [&](const std::tuple<std::shared_ptr<Atom>,Atom*> &e) {
                                                                        return (std::get<1>(e) == endB_ptr);
                                                                    })));

            Bond* newbond = new Bond(endA,endB,ptr->getBondID());
            newbond->setBondType(ptr->getBondType());
            std::shared_ptr<Bond>& inserted_bond = newmol.bonds.emplace_back(newbond);
            // PublicizeToAtom uses enable_shared_from_this, which means the previous statement, creating a shared_ptr
            // must occur before the following statement, lest throw std::bad_weak_ptr
            inserted_bond->publicizeToAtom();
        }

        newmol.rwmol.reset(new RDKit::RWMol(*rwmol));

        return newmol;
    }

    iConformer Molecule::getInitialConformer() {
        iConformer conformer;
        auto rdkit_conformer = (*(this->rwmol->beginConformers()));
        for (auto atom_it = this->rwmol->beginAtoms(); atom_it != this->rwmol->endAtoms(); ++atom_it) {
            assert((*atom_it)->getAtomicNum() != 1); // Conformer should not contain Hs (which were removed earlier)

            conformer.type.push_back(static_cast<unsigned char>((*atom_it)->getAtomicNum()));

            const RDGeom::Point3D &position = rdkit_conformer->getAtomPos((*atom_it)->getIdx());
            conformer.x.push_back(position.x);
            conformer.y.push_back(position.y);
            conformer.z.push_back(position.z);

            conformer.atomicRadius.push_back(atomTypeToAtomicRadius(
                    stringToAtomType(boost::to_upper_copy<std::string>((*atom_it)->getSymbol()))));
        }
        return conformer;
    }


}