//
// Created by eliane on 11/11/18.
//

#include "Molecule.h"
#include "Atom.h"

namespace SmolDock {

    Molecule::Molecule() = default;


    bool Molecule::populateFromSMILES(const std::string &smiles, unsigned int seed) {
        RDKit::RWMol *a = nullptr;
        try {
            /* This can throw */
            a = RDKit::SmilesToMol(smiles);
        }
        catch (...) {
            std::cerr << "[!] Error in constructor Molecule::Molecule(const std::string &smiles)" << std::endl;
            std::cerr << "[!] for smiles = " << smiles << std::endl;
            std::cerr << "[!] The SMILES string was not correctly parsed by RDKit." << std::endl;
            std::cerr << "[!] This often indicated a malformed SMILES. (but not always, RDKit has parsing bugs)"
                      << std::endl;
            return false;
        }
        if (a == nullptr) {
            std::cerr << "[!] Error in constructor Molecule::Molecule(const std::string &smiles)" << std::endl;
            std::cerr << "[!] for smiles = " << smiles << std::endl;
            std::cerr << "[!] The SMILES string was not correctly parsed by RDKit." << std::endl;
            std::cerr << "[!] This often indicated a malformed SMILES. (but not always, RDKit has parsing bugs)"
                      << std::endl;
            return false;
        }

        // If this does not throw : we have a mol in *a
        RDKit::MolOps::addHs(*a);

        rwmol.reset(a);

        int conformer_id = RDKit::DGeomHelpers::EmbedMolecule(*rwmol, 10, seed, true);
        RDKit::Conformer &starting_conformer = rwmol->getConformer(conformer_id);

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
                default:
                    std::cout << "[!] Unsupported bond type" << std::endl;
                    std::cout << "[!] Enum bond type = " << (*bond_it)->getBondType() << std::endl;
                    return false;
                    // break;
            }
            new_bond->publicizeToAtom();
            bonds.push_back(new_bond);
        }

        this->smiles = smiles;

        return true;
    }


    std::shared_ptr<RDKit::RWMol> Molecule::getInternalRWMol() const {
        return rwmol;
    }

    bool Molecule::generateConformer(iConformer &conformer, int seed) {

        int conformer_id = RDKit::DGeomHelpers::EmbedMolecule(*rwmol, 10, seed, false);

        if (conformer_id == -1) // Failed to generate
            return false;

        conformer.atoms_vect = std::make_unique<std::vector<iAtom> >();
        conformer.atoms_vect->reserve(this->atoms.size());

        RDKit::Conformer &rdkit_conformer = rwmol->getConformer(conformer_id);


        for (auto atom_it = rwmol->beginAtoms(); atom_it != rwmol->endAtoms(); ++atom_it) {
            iAtom atom;
            atom.atomicNum = (*atom_it)->getAtomicNum();

            const RDGeom::Point3D &position = rdkit_conformer.getAtomPos((*atom_it)->getIdx());
            atom.x = position.x;
            atom.y = position.y;
            atom.z = position.z;
            conformer.atoms_vect->push_back(std::move(atom));
        }

        return true;
    }

    unsigned int Molecule::generateConformers(std::vector<iConformer> *viConformers, unsigned int num, int seed) {


        std::vector<int> conformer_ids = RDKit::DGeomHelpers::EmbedMultipleConfs(*rwmol, // Molecule
                                                                                 num, // Num of conformer
                                                                                 30, // Max attempts
                                                                                 seed, // RNG seed
                                                                                 false); // Erase existing conformer
        if (conformer_ids.size() == 0)
            return 0; // Early failure case

        viConformers->reserve(viConformers->capacity() + conformer_ids.size());


        //std::vector<RDKit::Conformer> rdkit_conformers;
        //std::vector<iConformer> conformers;
        //rdkit_conformers.reserve(this->conformer_num);
        //conformers.reserve(this->conformer_num);

        unsigned int num_atoms = this->atoms.size();

        for (int i : conformer_ids) {
            iConformer conformer;
            RDKit::Conformer rdkit_conformer = this->rwmol->getConformer(i);

            conformer.atoms_vect = std::make_unique<std::vector<iAtom> >(static_cast<size_t>(num_atoms));

            for (auto atom_it = rwmol->beginAtoms(); atom_it != rwmol->endAtoms(); ++atom_it) {
                iAtom atom;
                atom.atomicNum = static_cast<unsigned char>((*atom_it)->getAtomicNum());

                const RDGeom::Point3D &position = rdkit_conformer.getAtomPos((*atom_it)->getIdx());
                atom.x = position.x;
                atom.y = position.y;
                atom.z = position.z;
                conformer.atoms_vect->push_back(std::move(atom));
            }

            viConformers->push_back(std::move(conformer));

        }

        return conformer_ids.size();
    }

    unsigned int Molecule::numberOfAtoms() {
        return atoms.size();
    }

    unsigned int Molecule::numberOfBonds() {
        return bonds.size();
    }


}