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

#include <tuple>
#include <algorithm>
#include <utility>


#include "Molecule.h"
#include "Atom.h"
#include "Utilities/PDBLigandUtils.h"

#include <GraphMol/Descriptors/Lipinski.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>

#include <ESBTL/default.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/weighted_atom_iterator.h>
#include <ESBTL/compressed_ifstream.h>


#include <boost/algorithm/string.hpp>
#include <boost/log/trivial.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>


#include <boost/graph/visitors.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/connected_components.hpp> // connected_components
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/copy.hpp>
#include <boost/bimap.hpp>


namespace SmolDock {

    Molecule::Molecule() {
        Molecule::LastMolID++;
        this->molID = LastMolID;
    };


    bool Molecule::populateFromSMILES(const std::string &smiles, unsigned int seed,
                                      std::vector<std::shared_ptr<InputModifier::InputModifier> > modifiers) {
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

        if (this->populateInternalAtomAndBondFromRWMol(seed, std::move(modifiers)) == false)
            return false;

        this->smiles = smiles;

        return true;
    }


    bool Molecule::generateConformer(iConformer &conformer, bool centroidNormalization, int seed) {

        // Hydrogen added for conformer gen
        RDKit::MolOps::addHs(*rwmol);
        int conformer_id = RDKit::DGeomHelpers::EmbedMolecule(*rwmol, 10, seed, false);
        RDKit::MMFF::MMFFOptimizeMolecule(*rwmol, 1000, "MMFF94s", 10.0, conformer_id);
        RDKit::MolOps::removeHs(*rwmol);

        if (conformer_id == -1) // Failed to generate
            return false;

        conformer = generateIConformerForGivenRDKitConformerID(conformer_id);

        return true;
    }

    unsigned int
    Molecule::generateConformers(std::vector<iConformer> &viConformers, unsigned int num, bool centroidNormalization,
                                 int seed) {

        RDKit::MolOps::addHs(*this->rwmol);
        std::vector<int> conformer_ids = RDKit::DGeomHelpers::EmbedMultipleConfs(*rwmol, // Molecule
                                                                                 num, // Num of conformer
                                                                                 30, // Max attempts
                                                                                 seed, // RNG seed
                                                                                 false); // Erase existing conformer
        RDKit::MolOps::removeHs(*rwmol);

        if (conformer_ids.size() == 0)
            return 0; // Early failure case


        viConformers.reserve(viConformers.capacity() + conformer_ids.size() + 1);



        //* // FIXME : uncomment this after fixing scoring function, and remove the initial conformer
        for (int i : conformer_ids) {
            RDKit::MMFF::MMFFOptimizeMolecule(*rwmol, 1000, "MMFF94s", 10.0, i);
            iConformer conformer = this->generateIConformerForGivenRDKitConformerID(i, centroidNormalization);
            viConformers.push_back(conformer);
        }
        //*/


        return conformer_ids.size();
    }

    unsigned int Molecule::numberOfAtoms() {
        return atoms.size();
    }

    unsigned int Molecule::numberOfBonds() {
        return bonds.size();
    }

    bool Molecule::populateFromPDBFile(const std::string &filename, const std::string &smiles_hint, unsigned int seed,
                                       std::vector<std::shared_ptr<InputModifier::InputModifier> > modifiers) {


        // Flavor = 1 --> ignore alternate location
        RDKit::RWMol *mol = RDKit::PDBFileToMol(filename, false, false, 1, true);

        if (mol == nullptr) {
            return false;
        }

        this->rwmol.reset(mol);

        if (smiles_hint != "") {
            AssignBondOrderFromTemplateSMILES(this->rwmol, smiles_hint);
        }

        if (this->populateInternalAtomAndBondFromRWMol(seed, modifiers) == false)
            return false;

        return true;
    }

    bool Molecule::populateFromMolFile(const std::string &filename, unsigned int seed,
                                       std::vector<std::shared_ptr<InputModifier::InputModifier> > modifiers) {


        RDKit::RWMol *mol = RDKit::MolFileToMol(filename, true, true);

        if (mol == nullptr) {
            return false;
        }

        this->rwmol.reset(mol);

        if (this->populateInternalAtomAndBondFromRWMol(seed, modifiers) == false)
            return false;

        return true;
    }

    bool Molecule::populateFromMol2File(const std::string &filename, unsigned int seed,
                                        std::vector<std::shared_ptr<InputModifier::InputModifier> > modifiers) {


        RDKit::RWMol *mol = RDKit::Mol2FileToMol(filename, true, true);

        if (mol == nullptr) {
            return false;
        }

        this->rwmol.reset(mol);

        if (this->populateInternalAtomAndBondFromRWMol(seed, modifiers) == false)
            return false;

        return true;
    }

    std::string Molecule::getResidueName() const {
        return this->residue_name;
    }

    void Molecule::setResidueName(const std::string &res_name) {
        this->residue_name = res_name.substr(0, 3);

    }

    bool Molecule::populateInternalAtomAndBondFromRWMol(unsigned int seed,
                                                        std::vector<std::shared_ptr<InputModifier::InputModifier> > modifiers) {

        // We care about Hs only during conformer generation
        // TODO: find out if we need to care at other moments...
        RDKit::MolOps::addHs(*(this->rwmol));

        this->initial_conformer_id = 0;
        if (this->rwmol->getNumConformers() == 0) // We need to generate the first conformer
        {
            this->initial_conformer_id = RDKit::DGeomHelpers::EmbedMolecule(*rwmol, 1000, seed, true);
            RDKit::MMFF::MMFFOptimizeMolecule(*rwmol, 1000, "MMFF94s", 10.0, this->initial_conformer_id);
        } else { // Else we will just use the first
            this->initial_conformer_id = (*(this->rwmol->beginConformers()))->getId();
        }

        RDKit::Conformer &starting_conformer = rwmol->getConformer(this->initial_conformer_id);

        // We remove them afterward
        RDKit::MolOps::removeHs((*(this->rwmol)));

        if (!this->rwmol->hasProp("_Name")) {
            this->rwmol->setProp("_Name", "ligand");
        }

        unsigned int currentAtomIndexInVector = 0;
        for (auto atom_it = rwmol->beginAtoms(); atom_it != rwmol->endAtoms(); ++atom_it) {
            std::shared_ptr<Atom> current_atom;
            // TODO : refactor this using the std::set of type-name-properties

            Atom::AtomType atomicNumToFind = static_cast<Atom::AtomType>((*atom_it)->getAtomicNum());

            auto it = std::find_if(Atom::AtomTypeLabel.begin(), Atom::AtomTypeLabel.end(),
                                   [atomicNumToFind](
                                           const std::tuple<Atom::AtomType, std::string, std::string, double> &e) {
                                       return std::get<0>(e) == atomicNumToFind;
                                   });

            if (it == Atom::AtomTypeLabel.end()) {
                BOOST_LOG_TRIVIAL(error) << "Unsupported atom type";
                BOOST_LOG_TRIVIAL(error) << "Atomic num = " << (*atom_it)->getAtomicNum();
                return false;
            }

            current_atom = std::make_shared<Atom>(std::get<0>(*it), (*atom_it)->getIdx());

            const RDGeom::Point3D &position = starting_conformer.getAtomPos((*atom_it)->getIdx());
            current_atom->setAtomPosition(std::make_tuple(position.x, position.y, position.z));

            atoms.push_back(current_atom);
            this->RDKitAtomIdxToAtomsPosInVector[(*atom_it)->getIdx()] = currentAtomIndexInVector;
            this->atomsPosInVectorToRDKitAtomIdx[currentAtomIndexInVector] = (*atom_it)->getIdx();
            currentAtomIndexInVector++;
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
                BOOST_LOG_TRIVIAL(error) << "Unable to find atom with ID = " << beginAtomID << " or " << endAtomID;
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
                    BOOST_LOG_TRIVIAL(error) << "Unsupported bond type";
                    BOOST_LOG_TRIVIAL(error) << "Enum bond type = " << (*bond_it)->getBondType();
                    return false;
                    // break;
            }
            new_bond->publicizeToAtom();
            bonds.push_back(new_bond);
        }

        this->numberOfRotatableBonds = RDKit::Descriptors::calcNumRotatableBonds((RDKit::ROMol) (*this->rwmol));


        unsigned int numAliphaticRing = RDKit::Descriptors::calcNumAliphaticRings((RDKit::ROMol) (*this->rwmol));
        if (numAliphaticRing > 0) {
            BOOST_LOG_TRIVIAL(error) << "Non-rigid (ie, non-aromatic) rings are unsupported yet";
            BOOST_LOG_TRIVIAL(error) << "  Found " << numAliphaticRing << " aliphatic ring(s)";
            BOOST_LOG_TRIVIAL(error) << "  Processing will continue but correctness is not guaranteed";
        }


        std::vector<std::tuple<int,int>> bondsListToAvoid;
        for (const auto &modifier : modifiers) {
                std::vector<std::tuple<int,int>> result = modifier->deselectRotatableBonds(this->rwmol);
                for(auto& item : result)
                {
                    bondsListToAvoid.push_back(item);
                }
        }


        // rotatable bond SMARTS from RdKit source (Lipinski.cpp)
        std::string strict_pattern =
                "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])("
                "[CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]="
                "[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-!@[!$(*#"
                "*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])"
                "[CH3])]";

        std::shared_ptr<RDKit::RWMol> rotatableBondPatt(RDKit::SmartsToMol(strict_pattern));

        std::vector<RDKit::MatchVectType> matchsRotBonds;

        if (RDKit::SubstructMatch(static_cast<const RDKit::ROMol &>(*this->rwmol), *rotatableBondPatt,
                                  matchsRotBonds)) {
            BOOST_LOG_TRIVIAL(debug) << "Rotatable bond(s) found";
        }

        assert(matchsRotBonds.size() == this->numberOfRotatableBonds);

        for (unsigned int i = 0; i < matchsRotBonds.size(); ++i) {
            if (matchsRotBonds[i].size() < 2) {
                BOOST_LOG_TRIVIAL(error) << "Rotatable bond pattern matching yielded 1 end only ??";
                continue;
            }
            unsigned int atom1Index = matchsRotBonds[i][0].second;
            unsigned int atom2Index = matchsRotBonds[i][1].second;
            const RDKit::Atom *atomMatched1 = this->rwmol->getAtomWithIdx(atom1Index);
            const RDKit::Atom *atomMatched2 = this->rwmol->getAtomWithIdx(atom2Index);

            assert(atomMatched1 != nullptr);
            assert(atomMatched2 != nullptr);

            BOOST_LOG_TRIVIAL(debug) << "Rotatable bond found : " << atomMatched1->getSymbol()
                                     << "-" << atomMatched2->getSymbol() << "(" << atomMatched1->getIdx() << "-"
                                     << atomMatched2->getIdx() << ")";


            auto avoidListMatchit = std::find_if(std::begin(bondsListToAvoid), std::end(bondsListToAvoid),
                                                 [atom1Index, atom2Index](const std::tuple<int,int> &e) {
                                                     // RotatableBondRemover and variant put both orientation in the bondsListToAvoid
                                                     // so we need only check one
                                                     return (atom1Index == std::get<0>(e)) && (atom2Index == std::get<1>(e));
                                                 });
            if(avoidListMatchit != std::end(bondsListToAvoid))
            {
                // We matched a bond from the list to avoid
                BOOST_LOG_TRIVIAL(info) << "  Rotatable bond " << atomMatched1->getSymbol()
                                        << "-" << atomMatched2->getSymbol() << "(" << atomMatched1->getIdx() << "-"
                                        << atomMatched2->getIdx() << ") not taken into account due to InputModifier flagging.";
                continue;
            }

            auto itFoundBond = std::find_if(std::begin(this->bonds), std::end(this->bonds),
                                            [atom1Index, atom2Index](const std::shared_ptr<Bond> &e) {
                                                // We might get endA and endB in either order
                                                return ((e->getEndA()->getAtomID() == atom1Index) &&
                                                        (e->getEndB()->getAtomID() == atom2Index))
                                                       || ((e->getEndB()->getAtomID() == atom1Index) &&
                                                           (e->getEndA()->getAtomID() == atom2Index));
                                            });
            if (itFoundBond == std::end(this->bonds)) {
                BOOST_LOG_TRIVIAL(error)
                    << "Rotatable bond not found in Molecule::bonds (but present in RDKit::RWMol) ?";
                return false;
            }

            std::vector<unsigned int> atomIdxOfRotatedAtoms;


            // TODO: while either side of the bond can be designated as "the part that rotate"
            // it would be useful to designate the smallest one, so as to do less math
            // or do some kind of tree traversal to rationalize all this "subsection to rotate" business
            // (this will help very much for quaternion, as they compose)

            std::function<void(const RDKit::Atom *, const RDKit::Atom *, const RDKit::Atom *)> exploreAdj;
            exploreAdj = ([&atomIdxOfRotatedAtoms, this, &exploreAdj]
                    (const RDKit::Atom *atom,
                     const RDKit::Atom *avoidThisAtom,
                     const RDKit::Atom *avoidThisOtherAtom) {
                RDKit::ROMol::ADJ_ITER nbr, end_nbr;
                boost::tie(nbr, end_nbr) = this->rwmol->getAtomNeighbors(atom);
                while (nbr != end_nbr) {
                    const RDKit::Atom *nbr_atom = (*this->rwmol)[*nbr];
                    if (nbr_atom == avoidThisAtom ||
                        nbr_atom == avoidThisOtherAtom ||
                        nbr_atom == atom ||
                        std::find(atomIdxOfRotatedAtoms.begin(), atomIdxOfRotatedAtoms.end(), nbr_atom->getIdx()) !=
                        atomIdxOfRotatedAtoms.end()) {
                        ++nbr;
                        continue;
                    }
                    atomIdxOfRotatedAtoms.push_back(nbr_atom->getIdx());
                     //BOOST_LOG_TRIVIAL(debug) << "AdjAtom from " << atom->getIdx() <<"(" <<atom->getSymbol() << ")" <<
                    //" : " << nbr_atom->getIdx() << " : " << nbr_atom->getSymbol();
                    exploreAdj(nbr_atom, avoidThisAtom, avoidThisOtherAtom);
                    ++nbr;
                }
            });

            exploreAdj(atomMatched1, atomMatched1, atomMatched2);

            this->rotatableBonds.emplace_back(
                    std::make_tuple(*itFoundBond, atom1Index, atom2Index, atomIdxOfRotatedAtoms));


        }





        // Post processing to assign variant

        // Processing to assign the apolar flag to carbon
        assignApolarCarbonFlag(this->atoms);

        // Final post processing : calling
        for (const auto &modifier : modifiers) {
            for (auto &atom: this->atoms) {
                modifier->postProcessAtomFromLigand(*atom);
            }
        }

        return true;
    }

    bool Molecule::updateAtomPositionsFromiConformer(const iConformer &conformer) {
        auto first_conformer = rwmol->beginConformers();
        RDKit::Conformer *newconformer = new RDKit::Conformer(**first_conformer);
        for (unsigned int i = 0; i < atoms.size(); ++i) {
            atoms[i]->setAtomPosition(
                    std::make_tuple(conformer.x[i], conformer.y[i], conformer.z[i])
            );
            newconformer->setAtomPos(i, {conformer.x[i], conformer.y[i], conformer.z[i]});

        }
        this->updatedConformerID = rwmol->addConformer(newconformer, true);

        return true;
    }

    Molecule Molecule::deepcopy() {
        Molecule newmol(*this);


        std::vector<std::tuple<std::shared_ptr<Atom>, Atom *> > new_old_lookup_table;
        newmol.atoms.clear();
        for (auto &ptr: this->atoms) {
            Atom *newatom_ptr = new Atom(*ptr);
            newatom_ptr->bonds.clear(); // We will restore it after rebuilding the bonds
            std::shared_ptr<Atom> &inserted = newmol.atoms.emplace_back(newatom_ptr);
            new_old_lookup_table.push_back(
                    std::make_tuple(inserted, ptr.get())); // Register new/old pair to rebuild bonds

        }

        newmol.bonds.clear();
        for (auto &ptr: this->bonds) {
            Atom *endA_ptr = ptr->getEndA().get();
            Atom *endB_ptr = ptr->getEndB().get();

            // Find matches in lookup table :
            std::shared_ptr<Atom> endA = std::get<0>(
                    *(std::find_if(std::begin(new_old_lookup_table), std::end(new_old_lookup_table),
                                   [&](const std::tuple<std::shared_ptr<Atom>, Atom *> &e) {
                                       return (std::get<1>(e) == endA_ptr);
                                   })));

            std::shared_ptr<Atom> endB = std::get<0>(
                    *(std::find_if(std::begin(new_old_lookup_table), std::end(new_old_lookup_table),
                                   [&](const std::tuple<std::shared_ptr<Atom>, Atom *> &e) {
                                       return (std::get<1>(e) == endB_ptr);
                                   })));

            Bond *newbond = new Bond(endA, endB, ptr->getBondID());
            newbond->setBondType(ptr->getBondType());
            std::shared_ptr<Bond> &inserted_bond = newmol.bonds.emplace_back(newbond);
            // PublicizeToAtom uses enable_shared_from_this, which means the previous statement, creating a shared_ptr
            // must occur before the following statement, lest throw std::bad_weak_ptr
            inserted_bond->publicizeToAtom();
        }

        newmol.rwmol.reset(new RDKit::RWMol(*rwmol));

        return newmol;
    }

    iConformer Molecule::getInitialConformer(bool centroidNormalization) const {
        unsigned int rdkit_first_conformer_id = (*(this->rwmol->beginConformers()))->getId();
        iConformer conformer = this->generateIConformerForGivenRDKitConformerID(rdkit_first_conformer_id,
                                                                                centroidNormalization);
        return conformer;
    }

    iConformer Molecule::generateIConformerForGivenRDKitConformerID(unsigned int id,
                                                                    bool centroidNormalization) const {
        iConformer conformer;
        RDKit::Conformer rdkit_conformer = this->rwmol->getConformer(id);

        unsigned long num_atoms = this->atoms.size();
        conformer.x.insert(conformer.x.begin(), num_atoms, 0.0);
        conformer.y.insert(conformer.y.begin(), num_atoms, 0.0);
        conformer.z.insert(conformer.z.begin(), num_atoms, 0.0);
        conformer.variant.insert(conformer.variant.begin(), num_atoms, 0);
        conformer.type.insert(conformer.type.begin(), num_atoms, 0);
        conformer.atomicRadius.insert(conformer.atomicRadius.begin(), num_atoms, 0.0);

        conformer.num_rotatable_bond = this->numberOfRotatableBonds;

        conformer.bondEnds1Index.insert(conformer.bondEnds1Index.begin(),
                                        conformer.num_rotatable_bond,
                                        0);
        conformer.bondEnds2Index.insert(conformer.bondEnds2Index.begin(),
                                        conformer.num_rotatable_bond,
                                        0);
        conformer.rotatableGroups.insert(conformer.rotatableGroups.begin(),
                                         conformer.num_rotatable_bond,
                                         std::vector<unsigned int>());


        for (unsigned int j = 0; j < conformer.num_rotatable_bond; ++j) {
            conformer.bondEnds1Index[j] = std::get<1>(rotatableBonds[j]);
            conformer.bondEnds2Index[j] = std::get<2>(rotatableBonds[j]);
            conformer.rotatableGroups[j] = std::get<3>(rotatableBonds[j]);
        }




        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::mean, tag::count, tag::min, tag::max> > acc_x, acc_y, acc_z;

        if (centroidNormalization) {
            for (auto atom_it = rwmol->beginAtoms(); atom_it != rwmol->endAtoms(); ++atom_it) {
                const RDGeom::Point3D &position = rdkit_conformer.getAtomPos((*atom_it)->getIdx());
                acc_x(position.x);
                acc_y(position.y);
                acc_z(position.z);
            }

            conformer.centroidNormalizingTransform = Eigen::Translation<double, 3>(mean(acc_x), mean(acc_y), mean(acc_z));

        } else {
            conformer.centroidNormalizingTransform = Eigen::Translation<double, 3>(0.0, 0.0, 0.0);
        }

        for (auto atom_it = rwmol->beginAtoms(); atom_it != rwmol->endAtoms(); ++atom_it) {

            // Find the corresponding Atom in this->atoms
            int atomIdx = (*atom_it)->getIdx();
            unsigned int atomsVectorIdx = this->RDKitAtomIdxToAtomsPosInVector.at(atomIdx);
            std::shared_ptr<Atom> current_Atom = this->atoms[atomsVectorIdx];

            assert(current_Atom != nullptr);

            conformer.type[(*atom_it)->getIdx()] = static_cast<unsigned char>((*atom_it)->getAtomicNum());
            conformer.variant[(*atom_it)->getIdx()] = static_cast<unsigned int>(current_Atom->variant);

            const RDGeom::Point3D &position = rdkit_conformer.getAtomPos((*atom_it)->getIdx());

            if (centroidNormalization) {
                conformer.x[(*atom_it)->getIdx()] = (position.x - conformer.centroidNormalizingTransform.x());
                conformer.y[(*atom_it)->getIdx()] = (position.y - conformer.centroidNormalizingTransform.y());
                conformer.z[(*atom_it)->getIdx()] = (position.z - conformer.centroidNormalizingTransform.z());
            } else {

                conformer.x[(*atom_it)->getIdx()] = position.x;
                conformer.y[(*atom_it)->getIdx()] = position.y;
                conformer.z[(*atom_it)->getIdx()] = position.z;
            }

            conformer.atomicRadius[(*atom_it)->getIdx()] = atomTypeToAtomicRadius(
                    stringToAtomType(boost::to_upper_copy<std::string>((*atom_it)->getSymbol()))
            );
        }

        assert(conformer.x.size() == num_atoms);
        assert(conformer.x.size() == conformer.y.size());
        assert(conformer.x.size() == conformer.z.size());
        assert(conformer.x.size() == conformer.type.size());
        assert(conformer.x.size() == conformer.atomicRadius.size());


        assert(conformer.bondEnds1Index.size() == conformer.num_rotatable_bond);
        assert(conformer.bondEnds1Index.size() == conformer.bondEnds2Index.size());
        assert(conformer.bondEnds1Index.size() == conformer.rotatableGroups.size());

        return conformer;
    }

    unsigned int Molecule::getNumRotatableBond() {
        return this->numberOfRotatableBonds;
    }


    unsigned int Molecule::LastMolID = 0;

    bool Molecule::operator==(const Molecule &rhs) const {
        return this->molID == rhs.molID;
    }

    bool Molecule::operator!=(const Molecule &rhs) const {
        return !(*this == rhs);
    }

    std::string Molecule::writeToMolBlock() {
        std::string molBlock = RDKit::MolToMolBlock(*(this->rwmol));
        return molBlock;
    }

    bool Molecule::writeToMolFile(const std::string &filename, bool overwrite) {
        std::ofstream molFile;

        std::string filecontent = this->writeToMolBlock();

        molFile.open(filename);
        molFile << filecontent << std::endl;
        molFile.close();

        BOOST_LOG_TRIVIAL(info) << "Wrote ligand to mol file : " << filename;
        return false;
    }

    unsigned int Molecule::applyAtomVariant(std::string smarts_pattern, Atom::AtomVariant variant) {
        std::shared_ptr<RDKit::RWMol> flaggingPattern(RDKit::SmartsToMol(smarts_pattern));

        std::vector<RDKit::MatchVectType> matchFlagged;

        RDKit::SubstructMatch(static_cast<const RDKit::ROMol &>(*this->rwmol), *flaggingPattern,
                              matchFlagged);

        if(matchFlagged.size() == 0)
            return 0;

        int queryAtomToBeFlaggedIdx = -1;
        RDKit::ROMol::VERTEX_ITER it , end;
        boost::tie( it , end ) = flaggingPattern->getVertices();
        while( it != end ) {
            const RDKit::Atom* atom = (*flaggingPattern)[*it];
            int atomMapNum = atom->getAtomMapNum();
            if(atomMapNum == 1)
            {
                queryAtomToBeFlaggedIdx = *it;
                break; // We found the atom which has SMARTS atom map 1, eg [C:1]ON, then C is matched
            }
            ++it;
        }

        unsigned int num_matched_atoms = 0;
        for (unsigned int i = 0; i < matchFlagged.size(); ++i) {
            for(unsigned int j = 0; j< matchFlagged[i].size(); ++j) {

                int queryIdx = std::get<0>(matchFlagged[i][j]);
                int atomIdx = std::get<1>(matchFlagged[i][j]);
                int vectorPos = this->RDKitAtomIdxToAtomsPosInVector[atomIdx];


                if (queryIdx == queryAtomToBeFlaggedIdx) // We flag
                {
                    Atom::AtomVariant var = this->atoms[vectorPos]->getAtomVariant();
                    var = var | variant;
                    this->atoms[vectorPos]->setAtomVariant(var);

                    num_matched_atoms++;
                }
            }

        }
        BOOST_LOG_TRIVIAL(info) << num_matched_atoms << " atoms were flagged with [" << atomVariantToString(variant) << "]";
        BOOST_LOG_TRIVIAL(debug) << "   Matching SMARTS : " << smarts_pattern;
        return num_matched_atoms;
    }


}