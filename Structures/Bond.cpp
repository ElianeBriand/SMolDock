//
// Created by eliane on 11/11/18.
//

#include "Bond.h"

namespace SmolDock {

    unsigned int Bond::nextBondID = 0;


    Bond::Bond(std::shared_ptr<Atom> atom_a, std::shared_ptr<Atom> atom_b) {
        this->bond_end_a = std::move(atom_a);
        this->bond_end_b = std::move(atom_b);
    }

    Bond::Bond(std::shared_ptr<Atom> atom_a, std::shared_ptr<Atom> atom_b, unsigned int id) {
        this->bond_end_a = std::move(atom_a);
        this->bond_end_b = std::move(atom_b);

        BondID = id;
        bondtype = defaultbond;
    }

    std::shared_ptr<Atom> Bond::getEndA() {
        return bond_end_a;
    }

    std::shared_ptr<Atom> Bond::getEndB() {
        return bond_end_b;
    }

    void Bond::publicizeToAtom() {
        // Publicize the bond to the atoms involved
        // (allow for easy atom based molecule traversal)
        bond_end_a->bonds.push_back(shared_from_this());
        bond_end_b->bonds.push_back(shared_from_this());
    }

    void Bond::setBondType(const Bond::BondType t) {
        bondtype = t;
    }

    Bond::BondType Bond::getBondType() const {
        return bondtype;
    }


}