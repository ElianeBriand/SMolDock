//
// Created by eliane on 11/11/18.
//

#include "Molecule.h"


namespace SmolDock {

    Molecule::Molecule() {

    }

    void Molecule::_dev_populateSampleMolecule() {

        /* Molecule created : C-O-C-OH */
        auto atom1 = std::make_shared<Atom>(Atom::carbon);
        this->atoms.push_back(atom1);

        auto atom1_h1 = std::make_shared<Atom>(Atom::hydrogen);
        auto bond_atom1_h1 = std::make_shared<Bond>(atom1, atom1_h1);
        this->atoms.push_back(atom1_h1);
        this->bonds.push_back(bond_atom1_h1);


        auto atom1_h2 = std::make_shared<Atom>(Atom::hydrogen);
        auto bond_atom1_h2 = std::make_shared<Bond>(atom1, atom1_h2);
        this->atoms.push_back(atom1_h2);
        this->bonds.push_back(bond_atom1_h2);

        auto atom1_h3 = std::make_shared<Atom>(Atom::hydrogen);
        auto bond_atom1_h3 = std::make_shared<Bond>(atom1, atom1_h3);
        this->atoms.push_back(atom1_h3);
        this->bonds.push_back(bond_atom1_h3);

        auto atom2 = std::make_shared<Atom>(Atom::oxygen);
        auto bond_atom1_atom2 = std::make_shared<Bond>(atom1, atom2);
        this->atoms.push_back(atom2);
        this->bonds.push_back(bond_atom1_atom2);

        auto atom3 = std::make_shared<Atom>(Atom::carbon);
        auto bond_atom2_atom3 = std::make_shared<Bond>(atom2, atom3);
        this->atoms.push_back(atom3);
        this->bonds.push_back(bond_atom2_atom3);

        auto atom3_h1 = std::make_shared<Atom>(Atom::hydrogen);
        auto bond_atom3_h1 = std::make_shared<Bond>(atom3, atom3_h1);
        this->atoms.push_back(atom3_h1);
        this->bonds.push_back(bond_atom3_h1);

        auto atom3_h2 = std::make_shared<Atom>(Atom::hydrogen);
        auto bond_atom3_h2 = std::make_shared<Bond>(atom3, atom3_h2);
        this->atoms.push_back(atom3_h2);
        this->bonds.push_back(bond_atom3_h2);

        auto atom4 = std::make_shared<Atom>(Atom::oxygen);
        auto bond_atom3_atom4 = std::make_shared<Bond>(atom3, atom4);
        this->atoms.push_back(atom4);
        this->bonds.push_back(bond_atom3_atom4);

        auto atom5 = std::make_shared<Atom>(Atom::hydrogen);
        auto bond_atom4_atom5 = std::make_shared<Bond>(atom4, atom5);
        this->atoms.push_back(atom5);
        this->bonds.push_back(bond_atom4_atom5);


        for (auto bond: this->bonds) {
            bond->publicizeToAtom();
        }


    }
}